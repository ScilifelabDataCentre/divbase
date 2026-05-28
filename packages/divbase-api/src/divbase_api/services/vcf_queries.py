import contextlib
import datetime
import os
import shlex
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import psutil
import structlog

from divbase_api.services.bcftools_helpers import (
    BCFTOOLS_CONTAINER_NAME,
    _raise_task_user_error_from_bcftools_stderr,
    ensure_csi_index,
    get_container_id,
    is_in_kubernetes,
    run_bcftools,
)
from divbase_lib.divbase_constants import QUERY_RESULTS_FILE_PREFIX
from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    BcftoolsPipeEmptyCommandError,
    BcftoolsPipeUnsupportedCommandError,
    TaskUserError,
)
from divbase_lib.utils import format_file_size, split_semicolon_bcftools_command_segments

logger = structlog.get_logger(__name__)

###
### Data structures for the query managers and their helper functions
###


@dataclass
class BCFToolsInput:
    """
    Contains the inputs required to run a bcftools query.
    """

    sample_and_filename_subset: list["SampleFileMapping"]
    sampleIDs: list[str]
    filenames: list[str]
    auto_sample_injection: bool = True


@dataclass
class SampleFileMapping:
    sample_id: str
    filename: str


@dataclass
class BCFToolsCommandConfig:
    command: str
    counter: int
    input_files: list[str]
    sample_subset: list[SampleFileMapping]
    output_temp_files: list[str]
    pipe_has_sample_placeholder: bool
    command_has_sample_placeholder: bool
    auto_sample_injection: bool


@dataclass
class ParsedCommandSegment:
    """
    Parsed representation of one semicolon-separated bcftools command segment.
    """

    position: int
    cmd: str
    args: list[str]  # includes the bcftools command and its options. E.g. ["view", "-s", "-r", "chr1:1-1000"]
    cmd_name: str
    normalized_segment: str


@dataclass
class ParsedViewRegionsOption:
    """
    Parsed representation of a single bcftools view regions argument.
    """

    is_regions_option: bool
    regions_value: str | None
    violation_reason: str | None
    next_index: int


###
### Bcftools command validation and parsing helper functions
###


def _parse_command_segments(command: str) -> list[ParsedCommandSegment]:
    """
    Parse a semicolon-separated bcftools command string into normalized segments.

    Raises TaskUserError on malformed quoting/syntax so callers can surface a user-friendly error.
    """
    segments: list[ParsedCommandSegment] = []

    for position, raw_cmd in enumerate(split_semicolon_bcftools_command_segments(command), start=1):
        cmd = raw_cmd.strip()
        if not cmd:
            # Defensive skip only. Empty segments are rejected by schema validation upstream.
            continue

        parse_error_message = f"Could not parse --command segment at position {position}: {cmd}"
        try:
            args = shlex.split(cmd)
        except ValueError:
            # shlex.split can raise ValueError for malformed quoting. We catch that and raise a TaskUserError with the position of the malformed segment.
            raise TaskUserError(parse_error_message) from None

        if not args:
            continue

        segments.append(
            ParsedCommandSegment(
                position=position,
                cmd=cmd,
                args=args,
                cmd_name=args[0],
                normalized_segment=shlex.join(args),
            )
        )

    return segments


def _parse_view_regions_option_argument(args: list[str], idx: int) -> ParsedViewRegionsOption:
    """
    Parse -r/--regions option forms at a given argument index in a bcftools view command.
    """
    arg = args[idx]

    def _build_result(
        *,
        is_regions_option: bool,
        regions_value: str | None = None,
        violation_reason: str | None = None,
        advance: int = 1,
    ) -> ParsedViewRegionsOption:
        return ParsedViewRegionsOption(
            is_regions_option=is_regions_option,
            regions_value=regions_value,
            violation_reason=violation_reason,
            next_index=idx + advance,
        )

    if arg in ("-r", "--regions"):
        next_arg = args[idx + 1] if idx + 1 < len(args) else None
        if next_arg is None or next_arg.startswith("-"):
            return _build_result(
                is_regions_option=True,
                violation_reason=(
                    "Option '-r/--regions' requires a non-empty region selector. "
                    "Use a value like '-r chr1:1-1000' or '--regions chr1:1-1000'."
                ),
            )
        return _build_result(
            is_regions_option=True,
            regions_value=next_arg,
            advance=2,
        )

    if arg.startswith("--regions="):
        regions_value = arg.split("=", 1)[1].strip()
        if not regions_value:
            return _build_result(
                is_regions_option=True,
                violation_reason=(
                    "Option '--regions=' requires a non-empty value. Use a value like '--regions=chr1:1-1000'."
                ),
            )
        return _build_result(
            is_regions_option=True,
            regions_value=regions_value,
        )

    if arg.startswith("-r") and len(arg) > 2 and not arg.startswith("--"):
        regions_value = arg[2:].strip()
        if not regions_value or regions_value.startswith("-"):
            return _build_result(
                is_regions_option=True,
                violation_reason=(
                    "Option '-r<REGIONS>' requires a non-empty region selector. "
                    "Use a value like '-rchr1:1-1000' or '-r chr1:1-1000'."
                ),
            )
        return _build_result(
            is_regions_option=True,
            regions_value=regions_value,
        )

    return _build_result(is_regions_option=False)


def _check_if_arg_matches_short_or_long_option(arg: str, short_opt: str, long_opt: str) -> bool:
    """Helper function to check if a given argument matches either the short or long version of a bcftools option."""
    if arg == short_opt:
        return True
    if arg.startswith(short_opt) and len(arg) > len(short_opt) and not arg.startswith("--"):
        return True
    if arg == long_opt:
        return True
    return bool(arg.startswith(f"{long_opt}="))


def _check_if_view_option_is_supported(arg: str) -> str | None:
    """Helper function to return the reason why a given bcftools view option is not supported by DivBase."""

    if _check_if_arg_matches_short_or_long_option(arg, "-h", "--header-only"):
        return (
            "Option '-h/--header-only' is not supported in DivBase queries. "
            "Instead use: `divbase-cli files stream <file.vcf.gz> | zcat | awk '/^##/ || /^#CHROM/ {print} !/^#/ {exit}'`"
        )
    if _check_if_arg_matches_short_or_long_option(arg, "-l", "--compression-level"):
        return "Option '-l/--compression-level' is handled by the DivBase server."
    if _check_if_arg_matches_short_or_long_option(arg, "-O", "--output-type"):
        return "Option '-O/--output-type' is handled by the DivBase server."
    if _check_if_arg_matches_short_or_long_option(arg, "-o", "--output"):
        return "Option '-o/--output' is handled by the DivBase server."
    if _check_if_arg_matches_short_or_long_option(arg, "-R", "--regions-file"):
        return "Option '-R/--regions-file' is not supported in DivBase queries because external filter files are not supported."
    if _check_if_arg_matches_short_or_long_option(arg, "-T", "--targets-file"):
        return "Option '-T/--targets-file' is not supported in DivBase queries because external filter files are not supported."
    if arg == "--threads" or arg.startswith("--threads="):
        return "Option '--threads' is handled by the DivBase server."
    if arg == "--verbosity" or arg.startswith("--verbosity="):
        return "Option '--verbosity' is handled by the DivBase server."
    if _check_if_arg_matches_short_or_long_option(arg, "-W", "--write-index"):
        return "Option '-W/--write-index' is handled by the DivBase server."
    if _check_if_arg_matches_short_or_long_option(arg, "-S", "--samples-file"):
        return (
            "Option '-S/--samples-file' is not supported in '--command'. "
            "Use `divbase-cli query vcf --samples-file` instead."
        )
    if _check_if_arg_matches_short_or_long_option(arg, "-f", "--apply-filters"):
        return "Option '-f/--apply-filters' is not supported in DivBase queries."

    return None


def _check_if_view_samples_option_violates_constraints(arg: str, all_samples: bool) -> str | None:
    """
    Return violation reason for -s/--samples usage in --command, or None if no violation.
    """
    if all_samples:
        return (
            "Option '-s/--samples' is not supported together with '--all-samples'. "
            "Use CLI '--samples' or '--samples-file' mode instead, or remove '-s/--samples' from '--command'."
        )

    if arg.startswith("--samples="):
        return (
            "Option '--samples=<LIST>' is not supported in '--command'. "
            "Set samples via DivBase CLI '--samples' or '--samples-file'. "
            "If you only want sample-based subsetting, use '--command \"view -s\"'."
        )

    if arg.startswith("-s") and arg != "-s" and not arg.startswith("--"):
        return (
            "Option '-s<LIST>' is not supported in '--command'. "
            "Set samples via DivBase CLI '--samples' or '--samples-file'. "
            "If you only want sample-based subsetting, use '--command \"view -s\"'."
        )

    return None


def _validate_view_segment(
    segment: ParsedCommandSegment,
    all_samples: bool,
) -> tuple[list[str], bool]:
    """
    Validate one parsed `bcftools view` command segment and return collected violations.
    """
    position = segment.position
    args = segment.args
    unsupported_view_options = []
    has_non_sample_filter_flag = False

    def _check_if_arg_looks_like_vcf_or_bcf_filename(arg: str) -> bool:
        """Reusable helper function to check if a given argument looks like a VCF or BCF filename, based on its extension."""
        normalized_arg = arg.lower()
        if normalized_arg.endswith(".vcf.gz"):
            return True
        if normalized_arg.endswith(".vcf"):
            return True
        return bool(normalized_arg.endswith(".bcf"))

    def _format_segment_argument_violation(argument: str, reason: str) -> str:
        """Reusable helper function to format a violation message for a given argument in a command segment."""
        return f"Pipe segment {position}, argument '{argument}': {reason}"

    filename_in_command_error_message = (
        "Do not provide VCF/BCF input filenames in '--command'; DivBase resolves input files from project dimensions."
    )
    stdin_in_command_error_message = (
        "Do not use stdin '-' in '--command'; DivBase resolves input files from project dimensions."
    )

    # Users must explicitly provide at least one view flag.
    # This guards against autoinject of -s for `--command "view"` when sample selection is used, which could lead to bugs.
    if not any(argument.startswith("-") for argument in args[1:]):
        for argument in args[1:]:
            if _check_if_arg_looks_like_vcf_or_bcf_filename(argument):
                unsupported_view_options.append(
                    _format_segment_argument_violation(
                        argument=argument,
                        reason=filename_in_command_error_message,
                    )
                )
        if unsupported_view_options:
            return unsupported_view_options, False
        raise TaskUserError(
            f"Pipe segment {position}: 'view' must include at least one bcftools view option flag "
            "(short or long form). "
            'If you only want sample-based subsetting, use --command "view -s".'
        )

    # Check for unsupported view options and for unsupported usage of -s/--samples in the command string.
    idx = 1
    while idx < len(args):
        arg = args[idx]
        next_idx = idx + 1

        if arg == "-":
            # Do not support bcftools stdin operator '-' for piping files into a bcftools pipe (e.g. 'cat file1.vcf | bcftools view -' or 'cat file1.vcf | bcftools view -r 1:1-10000 -')
            unsupported_view_options.append(
                _format_segment_argument_violation(
                    argument=arg,
                    reason=stdin_in_command_error_message,
                )
            )
            idx = next_idx
            continue

        reason_for_not_supported = _check_if_view_option_is_supported(arg=arg)
        if reason_for_not_supported is not None:
            unsupported_view_options.append(
                _format_segment_argument_violation(
                    argument=arg,
                    reason=reason_for_not_supported,
                )
            )
            idx = next_idx
            continue

        # Check for unsupported usage of -s/--samples in the command string.
        if _check_if_arg_matches_short_or_long_option(arg, "-s", "--samples"):
            next_arg_lookahead = args[idx + 1] if idx + 1 < len(args) else None
            sample_option_violation_reason = _check_if_view_samples_option_violates_constraints(
                arg=arg, all_samples=all_samples
            )
            if sample_option_violation_reason is not None:
                raise TaskUserError(f"Pipe segment {position}, argument '{arg}': {sample_option_violation_reason}")

            # Look ahead to the next argument after -s/--samples to check if it is an option
            # (beginning with - or --). Anything else is treated as user-provided sample names.
            # E.g. "view -s S1,S2" is invalid.
            if arg in ("-s", "--samples") and next_arg_lookahead is not None and not next_arg_lookahead.startswith("-"):
                raise TaskUserError(
                    f"Pipe segment {position}, argument '{arg}': "
                    "Do not provide sample names in '--command' via '-s/--samples'. "
                    "DivBase resolves sample IDs from '--tsv-filter', '--samples', or '--samples-file'. "
                    "If you only want sample-based subsetting, use '--command \"view -s\"'."
                )
            idx = next_idx
            continue

        parsed_region_option = _parse_view_regions_option_argument(args=args, idx=idx)
        next_idx = parsed_region_option.next_index
        if parsed_region_option.is_regions_option:
            if parsed_region_option.violation_reason:
                unsupported_view_options.append(
                    _format_segment_argument_violation(
                        argument=arg,
                        reason=parsed_region_option.violation_reason,
                    )
                )
            elif parsed_region_option.regions_value is not None and _check_if_arg_looks_like_vcf_or_bcf_filename(
                parsed_region_option.regions_value
            ):
                # For '-r/--regions' forms, reject filename-like region values so commands like
                # `view -r file.vcf.gz` do not look like user-provided input files.
                unsupported_view_options.append(
                    _format_segment_argument_violation(
                        argument=parsed_region_option.regions_value,
                        reason=filename_in_command_error_message,
                    )
                )
                idx = next_idx
                continue

        if _check_if_arg_looks_like_vcf_or_bcf_filename(arg):
            unsupported_view_options.append(
                _format_segment_argument_violation(
                    argument=arg,
                    reason=filename_in_command_error_message,
                )
            )
            idx = next_idx
            continue

        if all_samples and arg.startswith("-"):
            has_non_sample_filter_flag = True

        idx = next_idx

    return unsupported_view_options, has_non_sample_filter_flag


def validate_user_submitted_bcftools_command(command: str, all_samples: bool = False) -> str:
    """
    Validate that user-submitted bcftools command(s) are valid for DivBase.
    Intended to be run in the API layer to make early exits when needed.

    Supports two ways that users can submit "complex" bcftools commands:
    - Using semicolons to separate multiple command segments. E.g. "view -s -r 1:1-100; view -G" contains two semicolon-separated segments. bcftools manual reccommends this for more control over data processing/data ingrity.
    - Using a single command string with multiple bcftools options. E.g. "view -s -r 1:1-100 -G" contains multiple bcftools options. Will result is fewer bcftools subprocesses, but bcftools manual mentions that this **might** lead to unexpected results.

    Validation is done by iterating over each command segment separately. E.g. "view -s -r 1:1-100; view -G" contains two semicolon-separated segments.
    Some bcftools options have both a short and long version (e.g. -h and --header-only).
    The validation checks for both versions of the option.

    Returns the validated command string.
    """
    valid_commands = BcftoolsQueryManager.VALID_BCFTOOLS_COMMANDS
    segment_violations = []
    has_all_samples_and_at_least_one_view_option_that_is_not_samples = False
    normalized_segments_with_positions: dict[str, list[int]] = {}

    for segment in _parse_command_segments(command):
        # evaluate each command segment separately. E.g. "view -s -r 1:1-100; view -G" contains two semicolon-separated segments.
        position = segment.position
        cmd = segment.cmd
        cmd_name = segment.cmd_name

        if not cmd:
            # Empty command/empty ';' segments are already rejected by upstream Pydantic command field validation.
            continue

        if cmd_name not in valid_commands:
            raise TaskUserError(
                f"Unsupported bcftools command '{cmd_name}' at position {position}. "
                f"Only the following commands are supported: {', '.join(valid_commands)}"
            )

        # Duplication tracking.
        normalized_segment = segment.normalized_segment  # quoting/whitespace normalized by shlex.join
        normalized_segments_with_positions.setdefault(normalized_segment, []).append(position)

        segment_violations_for_segment, has_non_sample_filter_flag = _validate_view_segment(
            segment=segment,
            all_samples=all_samples,
        )
        segment_violations.extend(segment_violations_for_segment)
        has_all_samples_and_at_least_one_view_option_that_is_not_samples |= has_non_sample_filter_flag

    if segment_violations:
        details = "\n".join(f"  • {violation}" for violation in segment_violations)
        raise TaskUserError(f"Unsupported bcftools view option(s) found in '--command':\n{details}")

    duplicate_segments = {
        segment: positions for segment, positions in normalized_segments_with_positions.items() if len(positions) > 1
    }
    if duplicate_segments:
        duplicate_details = []
        for segment, positions in duplicate_segments.items():
            position_string = ", ".join(str(position) for position in positions)
            duplicate_details.append(f"  - Duplicate command segment '{segment}' found at positions: {position_string}")
        raise TaskUserError(
            "Duplicate bcftools command segment(s) found in '--command'. "
            "Please remove duplicate segments:\n" + "\n".join(duplicate_details)
        )

    if all_samples and not has_all_samples_and_at_least_one_view_option_that_is_not_samples:
        raise TaskUserError(
            "When using all-samples mode, --command must include at least one supported bcftools view option "
            "other than '-s/--samples' and none of the options blacklisted by DivBase. "
            "Examples: -r/--regions, -t/--targets, -i/--include, -e/--exclude, -q/--min-af, -Q/--max-af, -v/--types, -V/--exclude-types, -m/--min-alleles, -M/--max-alleles, -c/--min-ac, -C/--max-ac, -g/--genotype. "
            "Otherwise the query could potentially return all VCF data as is, and for that case it would be more efficient to download the dataset from the project instead with: divbase-cli files download-all. "
            "Please revise your command and try again."
        )

    return command


def extract_region_scaffolds_from_command(command: str) -> list[str]:
    """
    Extract scaffold/chromosome names from -r/--regions selectors in bcftools view command segments.

    Used by bcftools command validator, and by VCF query task to determin which VCF files to download if the user has specified region-based subsetting in the command.

    Supported forms:
    - -r chr1:1-1000
    - -rchr1:1-1000
    - --regions chr1:1-1000
    - --regions=chr1:1-1000
    """
    scaffolds: list[str] = []

    for segment in _parse_command_segments(command):
        if segment.cmd_name != "view":
            continue

        args = segment.args
        idx = 1
        while idx < len(args):
            parsed_region_option = _parse_view_regions_option_argument(args=args, idx=idx)
            idx = parsed_region_option.next_index

            if parsed_region_option.regions_value:
                for region in parsed_region_option.regions_value.split(","):
                    region = region.strip()
                    if not region:
                        continue
                    scaffold_name = region.split(":")[0].strip()
                    if scaffold_name:
                        scaffolds.append(scaffold_name)

    return scaffolds


###
### VCF query manager
###


class BcftoolsQueryManager:
    """
    A class that manages the execution of querys that require bcftools.

    Intended for use with a Celery architechture to run the queries as synchronous or asynchronous jobs.
    The bottom layer of the class - run_bcftools() - is designed for either being run inside a Celery worker container
    upon receiving a task from the queue (async job), or to be run inside the same Docker container as the Celery worker with
    `docker exec` instead of the queue (synchronous job). Either way, the class expects that the worker container with the
    name defined in BCFTOOLS_CONTAINER_NAME is running.

    NOTE! The current implementation does not handle starting the container.

    Users are expected to interact with this class through the CLI. The CLI layer handles Celery task management for this class.

    However, it is possible to run a query without Celery using the class directly, e.g.:
    `output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs)`

    The main entry point is the `execute_pipe` method, which takes a semicolon-separated string of
    bcftools commands and a dictionary of inputs. It builds a configuration structure for the commands,
    ensure that the relevant Docker container is available (for synchronous runs), and processes the commands in
    a structured manner.

    Current implementation is a "merge-last" strategy. In short, all subsetting and filtering commands
    are executed on each input files, temp files are created, and the results are merged at the end.
    The benefit of this approach is that merge operations can be memory intensive, and by subsetting first,
    the row-columns of the input files for the eventual merge operation are smaller.
    Downside is that the logic is more complex than in a "merge-first" strategy.

    The "merge-first" logic in based on nesting of two loops:
    - Outer loop: iterates over the input bcftools commands in the pipe, and passes each command to `run_current_command()`.
    - Inner loop: iterates over the input files for each command , and runs the command on each file using `run_bcftools()`.
    The results of the inner loop are temporary files that are passed to the next command in the outer loop.
    At the end of the outer loop, the temporary files are merged into a single output file using `merge_or_concat_bcftools_temp_files()`.
    The temporary files are cleaned up after the processing is done using `cleanup_temp_files()`.
    The class also provides a context manager for temporary file management to ensure that temporary files are cleaned
    up even if the processing fails or exits unexpectedly.

    """

    VALID_BCFTOOLS_COMMANDS = ["view"]  # white-list of valid bcftools commands to run in the pipe.

    def __init__(self, enable_subprocess_monitoring: bool = False):
        # this can be controlled by the worker config
        self.enable_subprocess_monitoring = enable_subprocess_monitoring

    def execute_pipe(self, command: str, bcftools_inputs: BCFToolsInput, job_id: int) -> tuple[Path, dict[str, float]]:
        """
        Main entrypoint for executing executing divbase queries that require bcftools.
        First calls on a method to build a structure of input parameters for bcftools, and then
        passes that on to another method that process the commands according to the "merge-last" strategy.

        Returns a tuple of (output_file, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses.
        """

        walltime_start = time.time()

        in_docker = os.path.exists("/.dockerenv")
        in_k8s = is_in_kubernetes()
        if not in_docker and not in_k8s:
            logger.info("Running outside Docker container, ensuring Docker container is available")
            try:
                container_id = get_container_id(BCFTOOLS_CONTAINER_NAME)
                logger.info(f"Found the required {BCFTOOLS_CONTAINER_NAME} container running with ID: {container_id}")
            except BcftoolsEnvironmentError:
                raise

        identifier = job_id if job_id else datetime.datetime.now().timestamp()

        commands_config_structure = self.build_commands_config(command, bcftools_inputs, identifier)
        output_file, metrics = self.process_bcftools_commands(commands_config_structure, identifier)

        walltime_end = time.time()
        metrics["walltime_seconds"] = walltime_end - walltime_start
        logger.info(f"Total bcftools pipeline walltime: {metrics['walltime_seconds']:.2f}s")

        return output_file, metrics

    @contextlib.contextmanager
    def temp_file_management(self):
        """Context manager to handle temporary file cleanup, even if processing fails/exits."""
        self.temp_files = []
        try:
            yield self
        finally:
            if self.temp_files:
                logger.info(f"Cleaning up {len(self.temp_files)} temporary files")
                self.cleanup_temp_files(self.temp_files)

    def build_commands_config(
        self, command: str, bcftools_inputs: BCFToolsInput, identifier: str = None
    ) -> list[BCFToolsCommandConfig]:
        """
        Method that builds a configuration structure for the bcftools commands based on the provided command string and inputs.
        The command string is expected to be a semicolon-separated list of bcftools commands.
        Each command is processed to create a list of dictionaries containing the command details,
        input files, sample subsets, and output temporary files. Returns a list of dictionaries
        where each dictionary represents a command configuration with the bcftools command, input files and temporary
        output files.
        """
        filenames = bcftools_inputs.filenames
        sample_and_filename_subset = bcftools_inputs.sample_and_filename_subset or []
        auto_sample_injection = bcftools_inputs.auto_sample_injection

        if not command or command.strip() == ";" or command.strip() == "":
            raise BcftoolsPipeEmptyCommandError()
        command_list = split_semicolon_bcftools_command_segments(command)
        pipe_has_sample_placeholder = any(
            self._command_has_sample_placeholder(cmd.strip()) for cmd in command_list if cmd.strip()
        )
        commands_config_structure: list[BCFToolsCommandConfig] = []
        current_inputs = filenames

        for c_counter, cmd in enumerate(command_list):
            cmd = cmd.strip()

            if not cmd:
                logger.warning(f"Skipping empty command at position {c_counter + 1} in command pipeline")
                continue

            try:
                cmd_name = shlex.split(cmd)[0]
            except ValueError:
                # shlex.split can raise ValueError for malformed quoting. If so, fallback on the raw cmd, which will raise BcftoolsPipeUnsupportedCommandError for the malformed cmd instead.
                cmd_name = cmd
            if cmd_name not in self.VALID_BCFTOOLS_COMMANDS:
                raise BcftoolsPipeUnsupportedCommandError(
                    command=cmd_name, position=c_counter + 1, valid_commands=self.VALID_BCFTOOLS_COMMANDS
                )

            output_temp_files = [
                f"temp_subset_{identifier}_{c_counter}_{f_counter}.bcf" for f_counter, _ in enumerate(current_inputs)
            ]

            command_details = BCFToolsCommandConfig(
                command=cmd,
                counter=c_counter,
                input_files=current_inputs,
                sample_subset=sample_and_filename_subset,
                output_temp_files=output_temp_files,
                pipe_has_sample_placeholder=pipe_has_sample_placeholder,
                command_has_sample_placeholder=self._command_has_sample_placeholder(cmd),
                auto_sample_injection=auto_sample_injection,
            )

            commands_config_structure.append(command_details)

            current_inputs = output_temp_files

        if not commands_config_structure:
            raise BcftoolsPipeEmptyCommandError()

        return commands_config_structure

    def process_bcftools_commands(
        self, commands_config: list[BCFToolsCommandConfig], identifier: str
    ) -> tuple[Path, dict[str, float]]:
        """
        Method that handles the outer loop of the merge-last strategy: it loops over each of commands in
        the input and passes them to the command runner run_current_command() (which in turn handles the inner loop:
        looping over VCF files). Once the outer loop is done, it calls a method to merge all temporary results files
        into a single output VCF file.

        Returns a tuple of (output_file, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses.
        """
        with self.temp_file_management() as temp_file_manager:
            logger.info(f"Loaded configuration with {len(commands_config)} commands in the pipe")

            final_output_temp_files: list[str] = []

            # Accumulate metrics across all commands
            total_cpu_seconds = 0.0
            max_peak_memory = 0
            all_memory_samples = []

            for cmd_config in commands_config:
                logger.info(f"Processing command #{cmd_config.counter + 1}: {cmd_config.command}")
                output_temp_files, cmd_metrics = self.run_current_command(cmd_config)
                temp_file_manager.temp_files.extend(output_temp_files)
                final_output_temp_files = output_temp_files

                total_cpu_seconds += cmd_metrics.get("cpu_seconds", 0)
                max_peak_memory = max(max_peak_memory, cmd_metrics.get("peak_memory_bytes", 0))
                if cmd_metrics.get("avg_memory_bytes", 0) > 0:
                    all_memory_samples.append(cmd_metrics["avg_memory_bytes"])

            output_file = self.merge_or_concat_bcftools_temp_files(final_output_temp_files, identifier)

            logger.info("bcftools processing completed successfully")

            avg_memory = sum(all_memory_samples) / len(all_memory_samples) if all_memory_samples else 0

            metrics = {
                "cpu_seconds": total_cpu_seconds,
                "peak_memory_bytes": max_peak_memory,
                "avg_memory_bytes": avg_memory,
            }

            return output_file, metrics

    def run_current_command(self, cmd_config: BCFToolsCommandConfig) -> tuple[list[str], dict[str, float]]:
        """
        Method that handles the inner loop of the merge-last strategy: for each command pass from the outer loop,
        it processes all given input VCF files individually by running the command on each file using run_bcftools().
        For each processed file, a temporary output file is created, which is then used as input for the next command
        in the outer loop. Each temporary output file is indexed with a .csi index file using ensure_csi_index().

        Automatically inserts resolved samples based on one of two policies:
        1) If the pipe contains a sample-placeholder (view -s / view --samples without sample values),
           inject at the placeholder position for those command segments.
        2) Otherwise, inject at the beginning of the first command segment.

        For instance, if the user submitted the command: "view -r chr1:1-1000" and the samples S1 and S2 are in the file,
        the command will be automatically converted to "view -s S1,S2 -r chr1:1-1000".

        Returns a tuple of (output_temp_files, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses executed in this command.
        """
        command = cmd_config.command
        input_files = cmd_config.input_files
        output_temp_files = cmd_config.output_temp_files
        sample_subset = cmd_config.sample_subset

        # Accumulate metrics across all bcftools subprocess calls
        total_cpu_seconds = 0.0
        peak_memory = 0
        memory_samples = []

        for f_counter, file in enumerate(input_files):
            temp_file = output_temp_files[f_counter]

            samples_in_file = []
            for record in sample_subset:
                if record.filename == file:
                    samples_in_file.append(record.sample_id)

            cmd_with_samples = command.strip()  # Use user-submitted command as starting point

            if cmd_config.auto_sample_injection and samples_in_file:
                samples_in_file_bcftools_formatted = ",".join(samples_in_file)
                if cmd_config.pipe_has_sample_placeholder:
                    if cmd_config.command_has_sample_placeholder:
                        cmd_with_samples = self._inject_samples_at_placeholder(
                            command=cmd_with_samples,
                            resolved_samples=samples_in_file_bcftools_formatted,
                        )
                # Explicit sample names inside --command are blocked by validate_user_submitted_bcftools_command().
                # So when no placeholder is present, inject samples at the first segment by default.
                elif cmd_config.counter == 0:
                    if cmd_with_samples == "view":
                        cmd_with_samples = f"view -s {samples_in_file_bcftools_formatted}"
                    elif cmd_with_samples.startswith("view "):
                        cmd_with_samples = f"view -s {samples_in_file_bcftools_formatted} {cmd_with_samples[5:]}"

            # TODO perhaps the default place for `view -s` should be the last place of the command? since it is faster on shorter files?
            # TODO: ensure backend strips `-s LIST_OF_SAMPLES` to just `-s`

            # Ensure source VCFs are indexed *before* command execution.
            ensure_csi_index(Path(file))

            # Output goes to temp_file (-Ou -o), not stdout, so the stdout pipe stays empty — no deadlock risk.
            formatted_cmd = f"{cmd_with_samples} {file} -Ou -o {temp_file}"

            # Run bcftools and optionally monitor the subprocess
            proc = run_bcftools(command=formatted_cmd, capture_stderr=True)

            if self.enable_subprocess_monitoring:
                current_cpu = 0.0  # Initialize to track CPU usage
                loop_iterations = 0

                try:
                    # Monitor the bcftools subprocess
                    bcftools_proc = psutil.Process(proc.pid)
                    cpu_start = bcftools_proc.cpu_times()
                    logger.info(
                        f"Started monitoring bcftools PID {proc.pid}, initial CPU: user={cpu_start.user}, system={cpu_start.system}"
                    )

                    # Take initial memory sample immediately
                    try:
                        mem = bcftools_proc.memory_info().rss
                        memory_samples.append(mem)
                        peak_memory = max(peak_memory, mem)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass

                    # Sample memory and CPU while process is running
                    while proc.poll() is None:
                        loop_iterations += 1
                        try:
                            mem = bcftools_proc.memory_info().rss
                            memory_samples.append(mem)
                            peak_memory = max(peak_memory, mem)

                            cpu_current = bcftools_proc.cpu_times()
                            current_cpu = (cpu_current.user + cpu_current.system) - (cpu_start.user + cpu_start.system)

                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            break
                        # Set interval between measurements
                        time.sleep(0.01)
                        # TODO sleep time could be longer to reduce overhead, but then fast processes might be missed

                    logger.info(
                        f"Monitoring loop completed after {loop_iterations} iterations, current_cpu={current_cpu:.6f}s"
                    )

                    self._wait_proc_and_check_return_code(proc=proc, command=formatted_cmd)

                    # Add CPU to total (will be 0.0 if process was too fast to measure)
                    total_cpu_seconds += current_cpu
                    logger.info(
                        f"Bcftools subprocess finished: CPU={current_cpu:.6f}s, Peak memory={peak_memory} bytes, Samples={len(memory_samples)}"
                    )

                except psutil.NoSuchProcess:
                    self._wait_proc_and_check_return_code(proc=proc, command=formatted_cmd)
                    logger.info("Bcftools process finished too quickly to monitor - caught NoSuchProcess")
            else:
                # Monitoring disabled, just wait for process to complete
                self._wait_proc_and_check_return_code(proc=proc, command=formatted_cmd)

                logger.info("Bcftools subprocess finished (monitoring disabled)")

            # Ensure temporary output VCFs are indexed *after* command execution.
            ensure_csi_index(Path(temp_file))
            self._log_file_size(temp_file)

        # Calculate average memory
        avg_memory = sum(memory_samples) / len(memory_samples) if memory_samples else 0

        metrics = {
            "cpu_seconds": total_cpu_seconds,
            "peak_memory_bytes": peak_memory,
            "avg_memory_bytes": avg_memory,
        }

        return output_temp_files, metrics

    def _command_has_sample_placeholder(self, command: str) -> bool:
        """
        Detect placeholder sample-option forms in a command segment.

        Assumes command has already passed validate_user_submitted_bcftools_command(), meaning that
        explicit sample-name forms in --command (for example '-s S1,S2', '-sS1,S2', '--samples S1,S2',
        '--samples=S1,S2') are have already been rejected upstream.

        Placeholder forms are:
        - view -s
        - view --samples
        - view -s <another-flag>
        - view --samples <another-flag>
        """
        args = shlex.split(command)
        for index, arg in enumerate(args[1:], start=1):
            if arg in ("-s", "--samples"):
                next_arg = args[index + 1] if index + 1 < len(args) else None
                if next_arg is None or next_arg.startswith("-"):
                    return True
        return False

    def _inject_samples_at_placeholder(self, command: str, resolved_samples: str) -> str:
        """
        Replace placeholder sample options (-s/--samples without values) with
        '-s <resolved_samples>' while preserving the rest of the command segment.
        """
        args = shlex.split(command)
        injected_args = []
        index = 0
        while index < len(args):
            arg = args[index]
            if arg in ("-s", "--samples"):
                next_arg = args[index + 1] if index + 1 < len(args) else None
                if next_arg is None or next_arg.startswith("-"):
                    injected_args.append("-s")
                    injected_args.append(resolved_samples)
                    index += 1
                    continue
            injected_args.append(arg)
            index += 1

        return shlex.join(injected_args)

    def merge_or_concat_bcftools_temp_files(self, output_temp_files: list[str], identifier: str) -> Path:
        """
        Helper method that merges the final temporary files produced by pipe_query_command into a single output file.

        This method adds some temp files list used by temp_file_management() for cleanup after processing. Note that
        the results file is not included in this list. Clean up of the results file is handled by tasks.py: only
        after successfull upload to S3 will the results file be deleted.
        """
        # TODO handle naming of output file better, e.g. by using a timestamp or a unique identifier

        unsorted_output_file = f"merged_unsorted_{identifier}.bcf"
        annotated_unsorted_output_file = f"merged_annotated_unsorted_{identifier}.bcf"
        divbase_header_for_vcf = "divbase_header.txt"  # Stored as temp txt on the worker for compatibility with bcftools annotate -h option, which requires a text file input for the header.
        self.temp_files.append(unsorted_output_file)
        self.temp_files.append(annotated_unsorted_output_file)
        self.temp_files.append(divbase_header_for_vcf)

        output_file = Path(f"{QUERY_RESULTS_FILE_PREFIX}{identifier}.vcf.gz")
        logger.info("Trying to determine if sample names overlap between temp files...")

        sample_names_per_VCF = self._get_all_sample_names_from_vcf_files(output_temp_files)
        sample_set_to_files = self._group_vcfs_by_sample_set(sample_names_per_VCF)
        non_overlapping_sample_names = self._check_non_overlapping_sample_names(sample_set_to_files)

        if len(output_temp_files) > 1:
            if non_overlapping_sample_names:
                logger.info("Sample names do not overlap between temp files, will continue with 'bcftools merge'")
                merge_command = f"merge -Ou -o {unsorted_output_file} {' '.join(output_temp_files)}"
                proc = run_bcftools(command=merge_command, capture_stderr=True)
                self._wait_proc_and_check_return_code(proc=proc, command=merge_command)
                logger.info(f"Merged all temporary files into '{unsorted_output_file}'.")
                self._log_file_size(unsorted_output_file)
            else:
                logger.info(
                    "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible."
                )
                temp_concat_files = []
                for sample_set, files in sample_set_to_files.items():
                    logger.debug(f"Processing sample set with {len(files)} files ({files}) and samples: {sample_set}")
                    if len(files) > 1:
                        logger.debug("Sample set occurs in multiple files, will concat these files.")
                        concat_temp = f"concat_{identifier}_{hash(sample_set)}.bcf"
                        concat_command = f"concat -Ou -o {concat_temp} {' '.join(files)}"
                        proc = run_bcftools(command=concat_command, capture_stderr=True)
                        self._wait_proc_and_check_return_code(proc=proc, command=concat_command)
                        temp_concat_files.append(concat_temp)
                        self.temp_files.append(concat_temp)
                        ensure_csi_index(Path(concat_temp))
                        self._log_file_size(concat_temp)
                    elif len(files) == 1:
                        logger.debug(
                            "Sample set only occurs in a single file, will use this file as is for merging in a downstream step."
                        )
                        temp_concat_files.append(files[0])
                if len(temp_concat_files) > 1:
                    merge_command = f"merge -Ou -o {unsorted_output_file} {' '.join(temp_concat_files)}"
                    proc = run_bcftools(command=merge_command, capture_stderr=True)
                    self._wait_proc_and_check_return_code(proc=proc, command=merge_command)
                    logger.info(f"Merged all files (including concatenated files) into '{unsorted_output_file}'.")
                    self._log_file_size(unsorted_output_file)
                elif len(temp_concat_files) == 1:
                    self._log_file_size(temp_concat_files[0])
                    os.rename(temp_concat_files[0], unsorted_output_file)
                    logger.info(
                        f"Only one file remained after concatenation, renamed this file to '{unsorted_output_file}'."
                    )
                    self._log_file_size(unsorted_output_file)
        elif len(output_temp_files) == 1:
            logger.info(f"Only one file was produced by the query, renamed this file to '{unsorted_output_file}'.")
            os.rename(output_temp_files[0], unsorted_output_file)
            self._log_file_size(unsorted_output_file)

        self._prepare_txt_with_divbase_header_for_vcf(header_filename=divbase_header_for_vcf)
        annotate_command = (
            f"annotate -h {divbase_header_for_vcf} -Ou -o {annotated_unsorted_output_file} {unsorted_output_file}"
        )
        proc = run_bcftools(command=annotate_command, capture_stderr=True)
        self._wait_proc_and_check_return_code(proc=proc, command=annotate_command)
        self._log_file_size(annotated_unsorted_output_file)

        sort_command = f"sort -Oz -o {output_file} {annotated_unsorted_output_file}"
        proc = run_bcftools(command=sort_command, capture_stderr=True)
        self._wait_proc_and_check_return_code(proc=proc, command=sort_command)
        self._log_file_size(str(output_file))
        logger.info(
            f"Sorting the results file to ensure proper order of variants. Final results are in '{output_file}'."
        )

        return output_file

    def cleanup_temp_files(self, output_temp_files: list[str]) -> None:
        """
        Helper method that handles deletion of all temporary files and their associated .csi index files
        that were generated during the pipe_query_command execution.
        """
        for temp_file in output_temp_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                index_file = f"{temp_file}.csi"
                if os.path.exists(index_file):
                    os.remove(index_file)
            except Exception as e:
                logger.error(f"Failed to delete temporary file or index {temp_file}: {e}")

    def _get_all_sample_names_from_vcf_files(self, output_temp_files: list[str]) -> dict[str, list[str]]:
        """
        Helper method that is used to determine if there are any sample names that recur across the temp files.
        If they do, bcftools concat is needed instead of bcftools merge.

        Use bcftools query to extract sample names from BCF/VCF files since it is designed for such operations.
        """
        sample_names_per_VCF = {}
        for vcf_file in output_temp_files:
            sample_extraction_command = f"query -l {vcf_file}"
            proc = run_bcftools(command=sample_extraction_command, capture_output=True)
            stdout, stderr = proc.communicate()

            if proc.returncode != 0:
                stderr = (stderr or "").strip()
                _raise_task_user_error_from_bcftools_stderr(
                    stderr=stderr,
                    operation="extracting sample names from a temporary BCF file",
                    target=f"temporary file '{vcf_file}'",
                )
                raise BcftoolsCommandError(
                    command=sample_extraction_command,
                    error_details=stderr or f"Process exited with code {proc.returncode}",
                ) from None

            sample_names = (stdout or "").strip().split("\n")
            sample_names_per_VCF[vcf_file] = sample_names if sample_names != [""] else []

        return sample_names_per_VCF

    def _group_vcfs_by_sample_set(self, sample_names_per_VCF: dict[str, list[str]]) -> dict[tuple[str, ...], list[str]]:
        """
        Helper method that groups VCF files by their sample sets. VCF files that contain the same sample set
        (=completely overlapping samples) need to be combined using bcftools concat instead of bcftools merge.
        Each sample set is represented as a tuple of sample names (preserving order) that acts as a hashable dict key.
        sample_set_to_files stores all files that share the same ordered sample set.
        """
        sample_set_to_files = {}
        for vcf_file, sample_list in sample_names_per_VCF.items():
            sample_set = tuple(sample_list)
            sample_set_to_files.setdefault(sample_set, []).append(vcf_file)
        return sample_set_to_files

    def _check_non_overlapping_sample_names(self, sample_set_to_files: dict[tuple[str, ...], list[str]]) -> bool:
        """
        Helper method that looks at a mapping of sample set to VCF files and checks for non-overlapping sample names.
        Simply put, if any sample set in the input dict has more than one file, samples overlap between files
        """
        return not any(len(files) > 1 for files in sample_set_to_files.values())

    def _prepare_txt_with_divbase_header_for_vcf(self, header_filename: str) -> None:
        """
        The command 'bcftools annotate -h' can be used to append a custom line to the header of a VCF file.
        The command requires a text file. This method checks if such a file exists, and if not, creates it.
        The 'annotate' command expects the pattern '##key=value' in the text file.
        """
        try:
            with open(header_filename, "w", encoding="utf-8") as fh:
                timestamp_in_format_used_by_bcftools = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
                fh.write(
                    f'##DivBase_created="This is a results file created by a DivBase query; Date={timestamp_in_format_used_by_bcftools}"\n'
                )
                logger.info("Added header to VCF saying that this results file was created with DivBase.")
        except Exception as e:
            logger.warning(f"Could not write {header_filename}: {e}")

    def _log_file_size(self, file_path: str):
        """
        Log the size of the a given file in a human-readable format.
        # TODO consider changing to logger debug later in the dev process
        """
        try:
            size_bytes = os.path.getsize(file_path)
            formatted_size = format_file_size(size_bytes)
            logger.info(f"File '{file_path}' size: {formatted_size}")
        except Exception as e:
            logger.warning(f"Could not determine size of file '{file_path}': {e}")

    def _wait_proc_and_check_return_code(self, proc: subprocess.Popen, command: str) -> None:
        """
        Wait for a bcftools subprocess and raise if it failed.
        """
        proc_stdout = getattr(proc, "stdout", None)
        proc_stderr = getattr(proc, "stderr", None)
        if proc_stdout is not None or proc_stderr is not None:
            _, stderr = proc.communicate()
            returncode = proc.returncode
        else:
            returncode = proc.wait()
            stderr = ""

        if returncode != 0:
            stderr = (stderr or "").strip()
            _raise_task_user_error_from_bcftools_stderr(
                stderr=stderr,
                operation=f"running command 'bcftools {command}'",
                target="the current query pipeline step",
            )
            raise BcftoolsCommandError(
                command=command,
                error_details=stderr or f"Process exited with code {returncode}",
            ) from None
