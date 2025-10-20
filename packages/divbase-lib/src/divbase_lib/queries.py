"""
TODO - consider split metadata query and bcftools query into two separate modules.
"""

import contextlib
import datetime
import gzip
import logging
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    BcftoolsPipeEmptyCommandError,
    BcftoolsPipeUnsupportedCommandError,
    SidecarColumnNotFoundError,
    SidecarInvalidFilterError,
    SidecarNoDataLoadedError,
    VCFDimensionsFileMissingOrEmptyError,
)
from divbase_lib.s3_client import S3FileManager
from divbase_lib.vcf_dimension_indexing import VCFDimensionIndexManager

logger = logging.getLogger(__name__)


@dataclass
class SidecarQueryResult:
    """
    Hold the results of a query run on a sidecar metadata TSV file.
    """

    sample_and_filename_subset: List[Dict[str, str]]
    unique_sample_ids: List[str]
    unique_filenames: List[str]
    query_message: str


def run_sidecar_metadata_query(
    file: Path, filter_string: str = None, bucket_name: str = None, s3_file_manager: S3FileManager = None
) -> SidecarQueryResult:
    """
    Run a query on a sidecar metadata TSV file.
    Call SidecarQueryManager to filter on the metadata, then call VCFDimensionIndexManager to get the
    sample-file name mapping. Combine the results and return a SidecarQueryResult object.
    """

    sidecar_manager = SidecarQueryManager(file=file).run_query(filter_string=filter_string)
    query_message = sidecar_manager.query_message
    unique_sample_ids = sidecar_manager.get_unique_values("Sample_ID")

    dimensions_manager = VCFDimensionIndexManager(bucket_name=bucket_name, s3_file_manager=s3_file_manager)
    dimensions_info = dimensions_manager.get_dimensions_info()

    if not dimensions_manager.dimensions_info or not dimensions_manager.dimensions_info.get("dimensions"):
        raise VCFDimensionsFileMissingOrEmptyError(dimensions_manager.bucket_name)

    sample_and_filename_subset = []
    unique_filenames = set()
    for entry in dimensions_info.get("dimensions", []):
        filename = entry["filename"]
        for sample_id in entry["dimensions"]["sample_names"]:
            if sample_id in unique_sample_ids:
                sample_and_filename_subset.append({"Sample_ID": sample_id, "Filename": filename})
                unique_filenames.add(filename)

    return SidecarQueryResult(
        sample_and_filename_subset=sample_and_filename_subset,
        unique_sample_ids=list(unique_sample_ids),
        unique_filenames=list(unique_filenames),
        query_message=query_message,
    )


@dataclass
class BCFToolsInput:
    """
    Contains the inputs required to run a bcftools query.
    """

    sample_and_filename_subset: List[Dict[str, str]]
    sampleIDs: List[str]
    filenames: List[str]


class BcftoolsQueryManager:
    """
    A class that manages the execution of querys that require bcftools.

    # TODO - support different file paths for input files.

    Intended for use with a Celery architechture to run the queries as synchronous or asynchronous jobs.
    The bottom layer of the class - self.run_bcftools() - is designed for either being run inside a Celery worker container
    upon receiving a task from the queue (async job), or to be run inside the same Docker container as the Celery worker with
    `docker exec` instead of the queue (synchronous job). Either way, the class expects that the worker container with the
    name defined in CONTAINER_NAME is running.

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
    CONTAINER_NAME = "divbase-worker-quick-1"  # for synchronous tasks, use this container name to find the container ID

    def execute_pipe(self, command: str, bcftools_inputs: dict, task_id: str = None) -> str:
        """
        Main entrypoint for executing executing divbase queries that require bcftools.
        First calls on a method to build a structure of input parameters for bcftools, and then
        passes that on to another method that process the commands according to the "merge-last" strategy.
        """

        in_docker = os.path.exists("/.dockerenv")
        in_k8s = self._is_in_kubernetes()
        if not in_docker and not in_k8s:
            logger.info("Running outside Docker container, ensuring Docker container is available")
            try:
                get_container_id = self.get_container_id(self.CONTAINER_NAME)
                logger.info(f"Found the required {self.CONTAINER_NAME} container running with ID: {get_container_id}")
            except BcftoolsEnvironmentError:
                raise

        identifier = task_id if task_id else datetime.datetime.now().timestamp()

        commands_config_structure = self.build_commands_config(command, bcftools_inputs, identifier)
        output_file = self.process_bcftools_commands(commands_config_structure, identifier)
        return output_file

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
        self, command: str, bcftools_inputs: Dict[str, Any], identifier: str = None
    ) -> List[Dict[str, Any]]:
        """
        Method that builds a configuration structure for the bcftools commands based on the provided command string and inputs.
        The command string is expected to be a semicolon-separated list of bcftools commands.
        Each command is processed to create a list of dictionaries containing the command details,
        input files, sample subsets, and output temporary files. Returns a list of dictionaries
        where each dictionary represents a command configuration with the bcftools command, input files and temporary
        output files.
        """
        filenames = bcftools_inputs.get("filenames")
        sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")

        if not command or command.strip() == ";" or command.strip() == "":
            raise BcftoolsPipeEmptyCommandError()
        command_list = command.split(";")
        commands_config_structure = []
        current_inputs = filenames

        for c_counter, cmd in enumerate(command_list):
            cmd = cmd.strip()

            if not cmd:
                logger.warning(f"Skipping empty command at position {c_counter + 1} in command pipeline")
                continue

            cmd_name = cmd.split()[0] if cmd and " " in cmd else cmd
            if cmd_name not in self.VALID_BCFTOOLS_COMMANDS:
                raise BcftoolsPipeUnsupportedCommandError(
                    command=cmd_name, position=c_counter + 1, valid_commands=self.VALID_BCFTOOLS_COMMANDS
                )

            output_temp_files = [
                f"temp_subset_{identifier}_{c_counter}_{f_counter}.vcf.gz" for f_counter, _ in enumerate(current_inputs)
            ]

            command_details = {
                "command": cmd,
                "counter": c_counter,
                "input_files": current_inputs,
                "sample_subset": sample_and_filename_subset,
                "output_temp_files": output_temp_files,
            }

            commands_config_structure.append(command_details)

            current_inputs = output_temp_files

        if not commands_config_structure:
            logger.error("No valid commands provided in input string")

        return commands_config_structure

    def process_bcftools_commands(self, commands_config: List[Dict[str, Any]], identifier: str) -> str:
        """
        Method that handles the outer loop of the merge-last strategy: it loops over each of commands in
        the input and passes them to the command runner run_current_command() (which in turn handles the inner loop:
        looping over VCF files). Once the outer loop is done, it calls a method to merge all temporary results files
        into a single output VCF file.
        """
        with self.temp_file_management() as temp_file_manager:
            logger.info(f"Loaded configuration with {len(commands_config)} commands in the pipe")

            final_output_temp_files = None

            for cmd_config in commands_config:
                logger.info(f"Processing command #{cmd_config['counter'] + 1}: {cmd_config['command']}")
                output_temp_files = self.run_current_command(cmd_config)
                temp_file_manager.temp_files.extend(output_temp_files)
                final_output_temp_files = output_temp_files

            output_file = self.merge_or_concat_bcftools_temp_files(final_output_temp_files, identifier)

            logger.info("bcftools processing completed successfully")

            return output_file

    def run_current_command(self, cmd_config: Dict[str, Any]) -> List[str]:
        """
        Method that handles the inner loop of the merge-last strategy: for each command pass from the outer loop,
        it processes all given input VCF files individually by running the command on each file using run_bcftools().
        For each processed file, a temporary output file is created, which is then used as input for the next command
        in the outer loop. Each temporary output file is indexed with a .csi index file using ensure_csi_index().
        The method returns a list of output temporary files created by running the command on each input file.
        """
        command = cmd_config["command"]
        input_files = cmd_config["input_files"]
        output_temp_files = cmd_config["output_temp_files"]
        sample_subset = cmd_config["sample_subset"]

        for f_counter, file in enumerate(input_files):
            temp_file = output_temp_files[f_counter]

            samples_in_file = []
            for record in sample_subset:
                if record["Filename"] == file:
                    samples_in_file.append(record["Sample_ID"])

            samples_in_file_bcftools_formatted = ",".join(samples_in_file)

            cmd_with_samples = command.strip().replace("SAMPLES", samples_in_file_bcftools_formatted)
            formatted_cmd = f"{cmd_with_samples} {file} -Oz -o {temp_file}"
            self.run_bcftools(command=formatted_cmd)
            self.ensure_csi_index(temp_file)
        return output_temp_files

    def run_bcftools(self, command: str) -> None:
        """
        Methid to run a bcftools command inside a Docker container that has bcftools installed.
        This method is specifically designed to with a Celery manager in an upper layer of the architechture.
        In short, the CLI layer handles the task management and submission to Celery, which can be either synchronous or asynchronous.

        Synchronous tasks (submitted with .apply()) skips the queue and tries to find the container id of the running container
        and runs the command inside it using `docker exec` and `subprocess.run`. I.e. the host machine executes the command
        inside the container.

        Asynchronous tasks (submitted with .apply_async()) are queued, picked up by the celery worker container and then execeuted
        inside the containter with `subprocess.run`. I.e. the container itself executes the command inside itself.

        To identify if the job is running async, the method evaluates if the current process is running inside a docker container by checking
        for the existence of the /.dockerenv file. If the file exists, it assumes that the command is run asynchronously inside the Celery worker
        container. If the file does not exist, it assumes that the command is run synchronously and tries to find the container ID of the running
        bcftools container using the get_container_id() method.
        """
        logger.info(f"Running: bcftools {command}")

        in_docker = os.path.exists("/.dockerenv")
        in_k8s = self._is_in_kubernetes()

        if in_docker or in_k8s:
            logger.debug("Running inside Celery worker Docker container, executing bcftools directly")
            try:
                subprocess.run(["bcftools"] + command.split(), check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run bcftools directly: {e}")
                raise BcftoolsCommandError(command=command, error_details=e) from e
        else:
            try:
                container_id = self.get_container_id(self.CONTAINER_NAME)
                logger.debug(f"Executing command in container with ID: {container_id}")
                docker_cmd = ["docker", "exec", container_id, "bcftools"] + command.split()
                subprocess.run(docker_cmd, check=True)
            except BcftoolsEnvironmentError:
                raise
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run bcftools in container: {e}")
                raise BcftoolsCommandError(command=command, error_details=e) from e

    def ensure_csi_index(self, file: str) -> None:
        """
        Helper method that ensures that the given VCF file has a .csi index. If not, create it using bcftools.

        (bcftools can sometimes handle VCF files that lack an index file, but for consistency
        it is better to create an index file for all VCF files that are created.)
        """
        index_file = f"{file}.csi"
        if not os.path.exists(index_file):
            index_command = f"index -f {file}"
            self.run_bcftools(command=index_command)

    def merge_or_concat_bcftools_temp_files(self, output_temp_files: List[str], identifier: str) -> str:
        """
        Helper method that merges the final temporary files produced by pipe_query_command into a single output file.


        # for all sets that have len(files) > 1, perform concat, save the temp filename to a new list
        # for all sets that have len(files) == 1, save the temp filename to a new list
        # for all the temp filenames in the new list, perform merge

        """
        # TODO handle naming of output file better, e.g. by using a timestamp or a unique identifier

        unsorted_output_file = f"merged_unsorted_{identifier}.vcf.gz"
        annotated_unsorted_output_file = f"merged_annotated_unsorted_{identifier}.vcf.gz"
        output_file = f"merged_{identifier}.vcf.gz"
        logger.info("Trying to determine if sample names overlap between temp files...")

        sample_names_per_VCF = self._get_all_sample_names_from_vcf_files(output_temp_files)
        sample_set_to_files = self._group_vcfs_by_sample_set(sample_names_per_VCF)
        non_overlapping_sample_names = self._check_non_overlapping_sample_names(sample_set_to_files)

        if len(output_temp_files) > 1:
            if non_overlapping_sample_names:
                logger.info("Sample names do not overlap between temp files, will continue with 'bcftools merge'")
                merge_command = f"merge --force-samples -Oz -o {unsorted_output_file} {' '.join(output_temp_files)}"
                # TODO double check if this should use output_temp_files or if that is an old remnant. the code below uses sample_set_to_files but that is perhaps to decide between concat and merge
                self.run_bcftools(command=merge_command)
                logger.info(f"Merged all temporary files into '{unsorted_output_file}'.")
            else:
                logger.info(
                    "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible."
                )
                temp_concat_files = []
                for sample_set, files in sample_set_to_files.items():
                    logger.debug(f"Processing sample set with {len(files)} files ({files}) and samples: {sample_set}")
                    if len(files) > 1:
                        logger.debug("Sample set occurs in multiple files, will concat these files.")
                        concat_temp = f"concat_{identifier}_{hash(sample_set)}.vcf.gz"
                        concat_command = f"concat -Oz -o {concat_temp} {' '.join(files)}"
                        self.run_bcftools(command=concat_command)
                        temp_concat_files.append(concat_temp)
                        self.temp_files.append(concat_temp)
                        self.ensure_csi_index(concat_temp)
                    elif len(files) == 1:
                        logger.debug(
                            "Sample set only occurs in a single file, will use this file as is for merging in a downstream step."
                        )
                        temp_concat_files.append(files[0])
                if len(temp_concat_files) > 1:
                    merge_command = f"merge --force-samples -Oz -o {unsorted_output_file} {' '.join(temp_concat_files)}"
                    self.run_bcftools(command=merge_command)
                    logger.info(f"Merged all files (including concatenated files) into '{unsorted_output_file}'.")
                elif len(temp_concat_files) == 1:
                    os.rename(temp_concat_files[0], unsorted_output_file)
                    logger.info(
                        f"Only one file remained after concatenation, renamed this file to '{unsorted_output_file}'."
                    )
        elif len(output_temp_files) == 1:
            logger.info(f"Only one file was produced by the query, renamed this file to '{unsorted_output_file}'.")
            os.rename(output_temp_files[0], unsorted_output_file)

        divbase_header_for_vcf = "divbase_header.txt"
        self._prepare_txt_with_divbase_header_for_vcf(header_filename=divbase_header_for_vcf)
        annotate_command = (
            f"annotate -h {divbase_header_for_vcf} -Oz -o {annotated_unsorted_output_file} {unsorted_output_file}"
        )
        self.run_bcftools(command=annotate_command)
        self.temp_files.append(unsorted_output_file)
        self.temp_files.append(divbase_header_for_vcf)

        sort_command = f"sort -Oz -o {output_file} {annotated_unsorted_output_file}"
        self.run_bcftools(command=sort_command)
        self.temp_files.append(annotated_unsorted_output_file)
        logger.info(
            f"Sorting the results file to ensure proper order of variants. Final results are in '{output_file}'."
        )

        return output_file

    def cleanup_temp_files(self, output_temp_files: List[str]) -> None:
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

    def get_container_id(self, container_name: str) -> str:
        """
        Helper method to get the container ID of the running bcftools Docker container.
        Raises BcftoolsEnvironmentError if the container is not found or if the command fails.
        The error can the be re-raised be the methods that call this method, e.g. run_bcftools() and execute_pipe().
        """
        try:
            result = subprocess.run(
                ["docker", "ps", "--filter", f"name={container_name}", "--format", "{{.ID}}"],
                capture_output=True,
                text=True,
                check=True,
            )
            return result.stdout.strip()
        except subprocess.SubprocessError as e:
            logger.error(f"Docker command failed: {e}")
            raise BcftoolsEnvironmentError(container_name) from e

    def _get_all_sample_names_from_vcf_files(self, output_temp_files: List[str]) -> dict[str, list[str]]:
        """
        Helper method that is used to determine if there are any sample names that recur across the temp files.
        If they do, bcftools concat is needed instead of bcftools merge.
        """
        sample_names_per_VCF = {}
        for vcf_file in output_temp_files:
            with gzip.open(vcf_file, "rt") as file:
                for line in file:
                    if line.startswith("#CHROM"):
                        header = line.strip().split("\t")
                        sample_names_per_VCF[vcf_file] = header[9:]
                        break

        return sample_names_per_VCF

    def _group_vcfs_by_sample_set(self, sample_names_per_VCF: dict[str, list[str]]) -> dict[tuple, list[str]]:
        """
        Helper method that groups VCF files by their sample sets. VCF files that contain the same sample set
        (=completely overlapping samples) need to be combined using bcftools concat instead of bcftools merge.
        Here, frozenset is used to create an immutable set of sample names for each VCF file.
        sample_set_to_files then stores all files that contain the same sample set.
        """
        sample_set_to_files = {}
        for vcf_file, sample_list in sample_names_per_VCF.items():
            sample_set = tuple(sample_list)
            sample_set_to_files.setdefault(sample_set, []).append(vcf_file)
        return sample_set_to_files

    def _check_non_overlapping_sample_names(self, sample_set_to_files: dict[frozenset, list[str]]) -> bool:
        """
        Helper method that looks at a mapping of sample set to VCF files and checks for non-overlapping sample names.
        Simply put, if any sample set in the input dict has more than one file, samples overlap between files
        """
        return not any(len(files) > 1 for files in sample_set_to_files.values())

    def _is_in_kubernetes(self) -> bool:
        """
        Check if k8s environment variable is set in containter.
        """
        return "KUBERNETES_SERVICE_HOST" in os.environ

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


class SidecarQueryManager:
    """
    A class that manages the execution queries on sidecar metadata files.

    Expects a TSV file with a header row and tab-separated values.
    The class provides methods to load the TSV file into a pandas DataFrame, run queries against the data,
    and retrieve unique values from specific columns.

    TODO - consider seperation of concerns: this class currently handles both loading the TSV file and running queries on it.
    TODO - some of the __init__ params are perhaps better as properties?
    """

    def __init__(self, file: Path):
        self.file = file
        self.filter_string = None
        self.df = None
        self.query_result = None
        self.query_message: str = ""
        self.load_file()

    def load_file(self) -> "SidecarQueryManager":
        """
        Method that loads the TSV file into a pandas DataFrame. Assumes that the first row is a header row, and that the file is tab-separated.
        Also removes any leading '#' characters from the column names

        If a sample exists in multiple files, it can be entered in the form: 'file1.vcf.gz,file2.vcf.gz' and all files will be extracted using pandas explode.
        Strip empty filenames if there e.g. are typos with trailing commas
        """
        # TODO: pandas will likely read all plain files to df, so perhaps there should be a check that the file is a TSV file? or at least has properly formatted tabular columns and rows?
        try:
            logger.info(f"Loading sidecar metadata file: {self.file}")
            self.df = pd.read_csv(self.file, sep="\t")
            self.df.columns = self.df.columns.str.lstrip("#")
            if "Sample_ID" not in self.df.columns:
                raise SidecarColumnNotFoundError("The 'Sample_ID' column is required in the metadata file.")

        except Exception as e:
            raise SidecarNoDataLoadedError(file_path=self.file, submethod="load_file") from e
        return self

    def run_query(self, filter_string: str = None) -> "SidecarQueryManager":
        """
        Method to run a query against the loaded data. The filter_string should be a semicolon-separated list of key:value pairs,
        where key is a column name and value is a comma-separated list of values to filter by.
        For example: "key1:value1,value2;key2:value3,value4".

        Summary of how different input filter values are handled:
        - If the filter_string is empty, all records are returned.
        - If the filter_string is None, an error is raised.
        - If the filter_string is not empty, the method filters the DataFrame based on the provided filter_string.
        - If any of the keys in the filter_string are not found in the DataFrame columns, a warning is logged and those conditions are skipped.
        - If none of the values in the filter_string are found in the DataFrame, a warning is logged and all records are returned.
        - If the filter_string is invalid, a SidecarInvalidFilterError is raised.

        The method returns the SidecarQueryManager instance with the query_result and query_message. The former is the filtered DataFrame results,
        and the latter is filter_string used for the query.
        """

        if self.df is None:
            raise SidecarNoDataLoadedError(
                file_path=self.file, submethod="run_query", error_details="No data loaded. Call load_file() first."
            )

        if filter_string is not None:
            self.filter_string = filter_string

        if self.filter_string == "":
            logger.warning("Empty filter provided - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "ALL RECORDS (no filter)"
            return self

        if self.filter_string is None:
            raise SidecarInvalidFilterError("Filter cannot be None. Use an empty string ('') if you want all records.")

        key_values = self.filter_string.split(";")
        filter_conditions = []

        for key_value in key_values:
            if not key_value.strip():
                continue
            try:
                key, values = key_value.split(":", 1)
                values_list = values.split(",")

                if key in self.df.columns:
                    condition = self.df[key].isin(values_list)
                    if not condition.any():
                        logger.warning(f"None of the values {values_list} were found in column '{key}'")
                    filter_conditions.append(condition)
                else:
                    logger.warning(f"Column '{key}' not found in the TSV file. Skipping this filter condition.")
            except Exception as e:
                raise SidecarInvalidFilterError(
                    f"Invalid filter format: '{key_value}'. Expected format 'key:value1,value2'"
                ) from e

        if filter_conditions:
            combined_condition = pd.Series(True, index=self.df.index)
            for condition in filter_conditions:
                combined_condition = combined_condition & condition

            self.query_result = self.df[combined_condition].copy()
            self.query_message = self.filter_string
        else:
            logger.warning("Invalid filter conditions found - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "Invalid filter conditions - returning ALL records"

        return self

    def get_unique_values(self, column: str) -> list:
        """
        Method to fetch unique values from a specific column in the query result. Intended to be invoked on a SidecarQueryManager
        instance after a query has been run with run_query().
        """
        if self.query_result is None:
            raise SidecarColumnNotFoundError("No query result available. Run run_query() first.")

        if column in self.query_result.columns:
            return self.query_result[column].unique().tolist()
        else:
            raise SidecarColumnNotFoundError(f"Column '{column}' not found in query result")
