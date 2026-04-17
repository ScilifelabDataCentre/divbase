"""Unit tests for VCF queries modes in tasks.py"""

from unittest.mock import patch

import pytest

from divbase_api.services.queries import (
    BCFToolsCommandConfig,
    BCFToolsInput,
    BcftoolsQueryManager,
    SampleFileMapping,
    extract_region_scaffolds_from_command,
    validate_user_submitted_bcftools_command,
)
from divbase_api.worker.crud_dimensions import ProjectVCFDimensionsData, ProjectVCFDimensionsEntry
from divbase_api.worker.tasks import (
    VCFQuerySampleSelectionMode,
    _determine_sample_selection_mode,
    _resolve_inputs_for_all_samples_mode,
    _resolve_inputs_for_cli_samples_mode,
)
from divbase_lib.exceptions import BcftoolsCommandError, TaskUserError


class TestDetermineSampleSelectionMode:
    @pytest.mark.parametrize(
        "tsv_filter,samples,all_samples,expected_mode",
        [
            ("Area:Northern Portugal", None, False, VCFQuerySampleSelectionMode.SAMPLE_METADATA_QUERY),
            (None, ["S1", "S2"], False, VCFQuerySampleSelectionMode.CLI_SAMPLES),
            (None, None, True, VCFQuerySampleSelectionMode.ALL_SAMPLES),
        ],
    )
    def test_determine_sample_selection_mode(self, tsv_filter, samples, all_samples, expected_mode):
        """Test that sample selection mode is determined correctly from tsv_filter/samples inputs."""
        mode = _determine_sample_selection_mode(tsv_filter=tsv_filter, samples=samples, all_samples=all_samples)
        assert mode == expected_mode

    def test_determine_sample_selection_mode_raises_when_no_selection_mode(self):
        with pytest.raises(TaskUserError, match="No sample-selection mode was provided"):
            _determine_sample_selection_mode(tsv_filter=None, samples=None, all_samples=False)


class TestValidateUserSubmittedBcftoolsCommand:
    @pytest.mark.parametrize(
        "command",
        [
            "view -r 21:15000000-25000000",
            "view -r chr1:1-1000",
            "view --regions=chr1:1-1000",
            "view -t chr1",
            "view -s",
            "view -s -r 21:15000000-25000000",
            "view --samples -i 'QUAL>20'",
        ],
    )
    def test_validate_user_submitted_bcftools_command_accepts_valid_commands(self, command):
        """Test that valid user-submitted bcftools commands pass command validation."""
        result = validate_user_submitted_bcftools_command(command)
        assert result == command

    @pytest.mark.parametrize(
        "command,expected_msg",
        [
            ("merge -m none", "Unsupported bcftools command 'merge'"),
            ("view -S samples.txt", "-S/--samples-file"),
            ("view -Ssamples.txt", "-S/--samples-file"),
            ("view --samples-file samples.txt", "-S/--samples-file"),
            ("view --samples-file=samples.txt", "-S/--samples-file"),
            ("view -s S1,S2", "Do not provide sample names in '--command'"),
            ("view -sS1,S2", "Option '-s<LIST>' is not supported"),
            ("view --samples S1,S2", "Do not provide sample names in '--command'"),
            ("view --samples=S1,S2", "Option '--samples=<LIST>' is not supported"),
            ("view --samples=", "Option '--samples=<LIST>' is not supported"),
            ("view file.vcf.gz", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view file.vcf", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view file.bcf", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view -", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view -r file.vcf.gz", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view --regions=file.vcf.gz", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view -t file.vcf", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view --targets=file.bcf", "Do not provide VCF/BCF input filenames in '--command'"),
            ("view", "must include at least one bcftools view option flag"),
            ("view -r", "requires a non-empty region selector"),
            ("view --regions", "requires a non-empty region selector"),
            ("view --regions=", "requires a non-empty value"),
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_invalid_commands(self, command, expected_msg):
        """Test that invalid user-submitted bcftools commands raise TaskUserError with expected message."""
        with pytest.raises(TaskUserError, match=expected_msg):
            validate_user_submitted_bcftools_command(command)

    def test_validate_user_submitted_bcftools_command_rejects_duplicate_segments(self):
        """Test that duplicate command segments in the same pipe are rejected."""
        command = "view -r 1:1-100; view -r 1:1-100"
        with pytest.raises(TaskUserError, match="Duplicate bcftools command segment"):
            validate_user_submitted_bcftools_command(command)

    def test_validate_user_submitted_bcftools_command_rejects_duplicate_segments_after_normalization(self):
        """Test that duplicate command segments are detected even if they differ in superficial whitespace/quoting."""
        command = "view -r '1:1-100'; view -r 1:1-100"
        with pytest.raises(TaskUserError, match="Duplicate bcftools command segment"):
            validate_user_submitted_bcftools_command(command)

    def test_validate_user_submitted_bcftools_command_rejects_unparseable_segment(self):
        """Test that unparseable user-submitted bcftools command segments raise TaskUserError."""
        with pytest.raises(TaskUserError, match="Could not parse --command segment at position 1"):
            validate_user_submitted_bcftools_command('view -r "chr1:1-1000')

    @pytest.mark.parametrize(
        "command,expected_msg",
        [
            ("view -h", "-h/--header-only"),
            ("view --header-only", "-h/--header-only"),
            ("view -l 9", "-l/--compression-level"),
            ("view -Oz", "-O/--output-type"),
            ("view --output-type=z", "-O/--output-type"),
            ("view -o output.vcf.gz", "-o/--output"),
            ("view --output=output.vcf.gz", "-o/--output"),
            ("view -R regions.txt", "-R/--regions-file"),
            ("view --regions-file=regions.txt", "-R/--regions-file"),
            ("view -T targets.txt", "-T/--targets-file"),
            ("view --targets-file=targets.txt", "-T/--targets-file"),
            ("view --threads=4", "--threads"),
            ("view --verbosity=2", "--verbosity"),
            ("view -Wz", "-W/--write-index"),
            ("view -W=z", "-W/--write-index"),
            ("view --write-index=tbi", "-W/--write-index"),
            ("view -f PASS", "-f/--apply-filters"),
            ("view --apply-filters=PASS", "-f/--apply-filters"),
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_blacklisted_view_options(self, command, expected_msg):
        """Test that user-submitted bcftools view commands with blacklisted options raise TaskUserError with expected message."""
        with pytest.raises(TaskUserError, match=expected_msg):
            validate_user_submitted_bcftools_command(command)

    def test_validate_user_submitted_bcftools_command_aggregates_blacklisted_options(self):
        """Test that user-submitted bcftools command with multiple blacklisted options raises TaskUserError with message that aggregates all violations."""
        command = "view -h; view -Oz --threads=2"

        with pytest.raises(TaskUserError) as exc_info:
            validate_user_submitted_bcftools_command(command)

        msg = str(exc_info.value)
        assert "Unsupported bcftools view option(s) found in '--command'" in msg
        assert "segment 1" in msg
        assert "-h/--header-only" in msg
        assert "segment 2" in msg
        assert "-O/--output-type" in msg
        assert "--threads" in msg

    def test_validate_user_submitted_bcftools_command_aggregates_filename_operands(self):
        """Test that user-submitted bcftools command with filename operands in multiple segments aggregates violations."""
        command = "view -r file.vcf.gz; view sample.bcf"

        with pytest.raises(TaskUserError) as exc_info:
            validate_user_submitted_bcftools_command(command)

        msg = str(exc_info.value)
        assert "Unsupported bcftools view option(s) found in '--command'" in msg
        assert "segment 1" in msg
        assert "segment 2" in msg
        assert "Do not provide VCF/BCF input filenames in '--command'" in msg

    @pytest.mark.parametrize(
        "command",
        [
            "view",
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_view_without_flag_in_all_samples_mode(self, command):
        """Test that bare 'view' is rejected in all-samples mode as well."""
        with pytest.raises(TaskUserError, match="must include at least one bcftools view option flag"):
            validate_user_submitted_bcftools_command(command, all_samples=True)

    @pytest.mark.parametrize(
        "command",
        [
            "view -r 21:15000000-25000000",
            "view --regions=21:15000000-25000000",
            "view -t 21",
            "view --targets=21",
            "view -i 'MAF>0.05'",
            "view --include='MAF>0.05'",
            "view -e 'QUAL<20'",
            "view --exclude='QUAL<20'",
            "view -g hom",
            "view --genotype het",
            "view -q 0.01",
            "view --max-af=0.9",
            "view -v snps",
            "view --exclude-types=indels",
            "view -A",
        ],
    )
    def test_validate_user_submitted_bcftools_command_accepts_subset_filters_for_all_samples(self, command):
        """Test that user-submitted bcftools command in all-samples mode that includes at least one non-sample-selection view option passes validation."""
        result = validate_user_submitted_bcftools_command(command, all_samples=True)
        assert isinstance(result, str)

    @pytest.mark.parametrize(
        "command",
        [
            "view -s",
            "view --samples",
            "view -s S1,S2 -r 21:15000000-25000000",
            "view -r 21:15000000-25000000; view --samples=S1,S2",
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_samples_option_in_all_samples_mode(self, command):
        """Test that user-submitted bcftools command in all-samples mode that includes -s/--samples option raises TaskUserError."""
        with pytest.raises(TaskUserError, match="-s/--samples"):
            validate_user_submitted_bcftools_command(command, all_samples=True)


class TestExtractRegionScaffoldsFromCommand:
    @pytest.mark.parametrize(
        "command,expected_scaffolds",
        [
            ("view -r chr1:1-1000", ["chr1"]),
            ("view -rchr1:1-1000", ["chr1"]),
            ("view --regions chr1:1-1000", ["chr1"]),
            ("view --regions=chr1:1-1000", ["chr1"]),
            ("view -s -r 1,4,6,21,24", ["1", "4", "6", "21", "24"]),
            ("view -i 'QUAL>20'; view --regions chr2:1-10,chr3:4-5", ["chr2", "chr3"]),
            ("view -R regions.txt", []),
            ("view -t chr1", []),
            ("view -r; view -r chr2:1-1000", ["chr2"]),
        ],
    )
    def test_extract_region_scaffolds_from_command(self, command, expected_scaffolds):
        """Test that region scaffolds are extracted from supported -r/--regions command forms."""
        assert extract_region_scaffolds_from_command(command) == expected_scaffolds

    def test_extract_region_scaffolds_from_command_raises_on_malformed_segment(self):
        """Test that malformed command segments raise TaskUserError with a clear message."""
        command = 'view -r "chr1:1-1000; view -r chr2:1-1000'
        with pytest.raises(TaskUserError, match="Could not parse --command segment at position 1"):
            extract_region_scaffolds_from_command(command)


class TestSamplesPlaceholderDetectionAndInjection:
    def test_command_has_sample_placeholder_detects_placeholder_forms(self):
        """Test that _command_has_sample_placeholder correctly detects various valid forms of sample placeholders in bcftools view commands."""
        manager = BcftoolsQueryManager()

        assert manager._command_has_sample_placeholder("view -s")
        assert manager._command_has_sample_placeholder("view --samples")
        assert manager._command_has_sample_placeholder("view -s -r 1:1000-2000")
        assert manager._command_has_sample_placeholder("view --samples -i 'QUAL>20'")
        assert not manager._command_has_sample_placeholder("view --samples S1,S2")
        assert not manager._command_has_sample_placeholder("view --samples=S1,S2")
        assert not manager._command_has_sample_placeholder("view -sS1,S2")

    def test_build_commands_config_marks_pipe_and_segment_placeholder_flags(self):
        """Test that build_commands_config correctly sets flags indicating presence of sample placeholders in command segments and pipes."""
        manager = BcftoolsQueryManager()
        bcftools_inputs = BCFToolsInput(
            filenames=["file1.vcf.gz"],
            sample_and_filename_subset=[
                SampleFileMapping(sample_id="S1", filename="file1.vcf.gz"),
                SampleFileMapping(sample_id="S2", filename="file1.vcf.gz"),
            ],
            sampleIDs=["S1", "S2"],
            auto_sample_injection=True,
        )

        config = manager.build_commands_config(
            command="view -r 1:1000-2000; view -s",
            bcftools_inputs=bcftools_inputs,
            identifier="job1",
        )

        assert len(config) == 2
        assert config[0].pipe_has_sample_placeholder is True
        assert config[1].pipe_has_sample_placeholder is True
        assert config[0].command_has_sample_placeholder is False
        assert config[1].command_has_sample_placeholder is True

    def test_inject_samples_at_placeholder_preserves_position_and_other_flags(self):
        """Test that _inject_samples_at_placeholder correctly injects sample names at the placeholder position without altering other command flags or structure."""
        manager = BcftoolsQueryManager()

        injected = manager._inject_samples_at_placeholder(
            command="view -s -r 1:1000-2000",
            resolved_samples="S1,S2",
        )
        assert injected == "view -s S1,S2 -r 1:1000-2000"

        injected_long = manager._inject_samples_at_placeholder(
            command="view --samples -i 'QUAL>20'",
            resolved_samples="S1,S2",
        )
        assert injected_long == "view -s S1,S2 -i 'QUAL>20'"

    def test_run_current_command_injects_samples_at_placeholder_position(self, monkeypatch):
        """Test that run_current_command correctly injects sample names into the bcftools command at the placeholder position when auto_sample_injection is enabled."""
        manager = BcftoolsQueryManager()
        manager.ENABLE_SUBPROCESS_MONITORING = False
        executed_commands = []

        class DummyProc:
            pid = 12345

            @staticmethod
            def wait():
                return 0

            @staticmethod
            def poll():
                return 0

        def fake_run_bcftools(command: str):
            executed_commands.append(command)
            return DummyProc()

        monkeypatch.setattr(manager, "run_bcftools", fake_run_bcftools)
        monkeypatch.setattr(manager, "ensure_csi_index", lambda _file_path: None)
        monkeypatch.setattr(manager, "_log_file_size", lambda _file_path: None)

        cmd_config = BCFToolsCommandConfig(
            command="view -s -r 1:1000-2000",
            counter=1,
            input_files=["file1.vcf.gz"],
            sample_subset=[
                SampleFileMapping(sample_id="S1", filename="file1.vcf.gz"),
                SampleFileMapping(sample_id="S2", filename="file1.vcf.gz"),
            ],
            output_temp_files=["temp_subset_job1_1_0.bcf"],
            pipe_has_sample_placeholder=True,
            command_has_sample_placeholder=True,
            auto_sample_injection=True,
        )

        output_files, _metrics = manager.run_current_command(cmd_config)

        assert output_files == ["temp_subset_job1_1_0.bcf"]
        assert executed_commands == ["view -s S1,S2 -r 1:1000-2000 file1.vcf.gz -Ou -o temp_subset_job1_1_0.bcf"]


class TestResolveInputsForCliSamplesMode:
    def test_resolve_inputs_for_cli_samples_mode(self):
        """Test that CLI sample-mode inputs are resolved to files and sample mappings correctly."""
        vcf_dimensions_data = ProjectVCFDimensionsData(
            project_id=1,
            vcf_file_count=2,
            vcf_files=[
                ProjectVCFDimensionsEntry(vcf_file_s3_key="fileA.vcf.gz", s3_version_id="v1", samples=["S1", "S2"]),
                ProjectVCFDimensionsEntry(vcf_file_s3_key="fileB.vcf.gz", s3_version_id="v2", samples=["S3"]),
            ],
        )
        result = _resolve_inputs_for_cli_samples_mode(
            samples=["S3", "S1", "S3"], vcf_dimensions_data=vcf_dimensions_data
        )

        assert result.files_to_download == ["fileA.vcf.gz", "fileB.vcf.gz"]
        assert result.sample_and_filename_subset == [
            SampleFileMapping(sample_id="S1", filename="fileA.vcf.gz"),
            SampleFileMapping(sample_id="S3", filename="fileB.vcf.gz"),
        ]
        assert result.unique_sample_ids == ["S3", "S1"]
        assert result.metadata_path is None

    def test_resolve_inputs_for_cli_samples_mode_missing_samples_raises(self):
        """Test that CLI sample-mode input raises TaskUserError when requested samples are missing."""
        vcf_dimensions_data = ProjectVCFDimensionsData(
            project_id=1,
            vcf_file_count=1,
            vcf_files=[
                ProjectVCFDimensionsEntry(vcf_file_s3_key="fileA.vcf.gz", s3_version_id="v1", samples=["S1", "S2"]),
            ],
        )
        with pytest.raises(TaskUserError, match="were not found in the project's dimensions index"):
            _resolve_inputs_for_cli_samples_mode(
                samples=["S1", "DOES_NOT_EXIST"], vcf_dimensions_data=vcf_dimensions_data
            )


class TestResolveInputsForAllSamplesMode:
    def test_resolve_inputs_for_all_samples_mode(self):
        """Test that all-samples mode resolves all files and sample mappings correctly."""
        vcf_dimensions_data = ProjectVCFDimensionsData(
            project_id=1,
            vcf_file_count=2,
            vcf_files=[
                ProjectVCFDimensionsEntry(vcf_file_s3_key="fileA.vcf.gz", s3_version_id="v1", samples=["S1", "S2"]),
                ProjectVCFDimensionsEntry(vcf_file_s3_key="fileB.vcf.gz", s3_version_id="v2", samples=["S3"]),
            ],
        )
        result = _resolve_inputs_for_all_samples_mode(vcf_dimensions_data=vcf_dimensions_data)

        assert result.files_to_download == ["fileA.vcf.gz", "fileB.vcf.gz"]
        assert result.sample_and_filename_subset == [
            SampleFileMapping(sample_id="S1", filename="fileA.vcf.gz"),
            SampleFileMapping(sample_id="S2", filename="fileA.vcf.gz"),
            SampleFileMapping(sample_id="S3", filename="fileB.vcf.gz"),
        ]
        assert set(result.unique_sample_ids) == {"S1", "S2", "S3"}
        assert result.metadata_path is None


class TestBcftoolsReturnCodeHandling:
    class DummyProc:
        """Dummy process class to simulate subprocess.Popen return value for testing return code handling."""

        def __init__(self, returncode: int):
            self._returncode = returncode

        def wait(self) -> int:
            return self._returncode

    def test_wait_proc_and_check_return_code_raises_on_non_zero(self):
        """Test that _wait_proc_and_check_return_code raises BcftoolsCommandError with appropriate message when process returns non-zero exit code."""
        manager = BcftoolsQueryManager()

        with pytest.raises(BcftoolsCommandError, match="Process exited with code 1") as exc_info:
            manager._wait_proc_and_check_return_code(
                proc=self.DummyProc(returncode=1),
                command="view -h file.vcf.gz",
            )

        assert "view -h file.vcf.gz" in str(exc_info.value)

    def test_wait_proc_and_check_return_code_accepts_zero(self):
        """Test that _wait_proc_and_check_return_code does not raise an error when process returns zero exit code."""
        manager = BcftoolsQueryManager()
        manager._wait_proc_and_check_return_code(
            proc=self.DummyProc(returncode=0),
            command="view -h file.vcf.gz",
        )

    def test_ensure_csi_index_raises_on_non_zero_return_code(self):
        """Test that ensure_csi_index raises BcftoolsCommandError with appropriate message when bcftools index command returns non-zero exit code."""
        manager = BcftoolsQueryManager()
        with (
            patch(
                "divbase_api.services.queries.os.path.exists", return_value=False
            ),  # Simulate missing index file to trigger indexing
            patch.object(manager, "run_bcftools", return_value=self.DummyProc(returncode=1)),
            pytest.raises(BcftoolsCommandError, match="Process exited with code 1") as exc_info,
        ):
            manager.ensure_csi_index("input.vcf.gz")

        assert "index -f input.vcf.gz" in str(exc_info.value)

    @pytest.mark.parametrize(
        "sample_names_map,non_overlapping,identifier,failing_prefix",
        [
            (
                {"file1.bcf": ["S1"], "file2.bcf": ["S2"]},
                True,
                "job1",
                "merge --force-samples",
            ),
            (
                {"file1.bcf": ["S1"], "file2.bcf": ["S1"]},
                False,
                "job2",
                "concat -Ou -o",
            ),
            (
                {"file1.bcf": ["S1"]},
                True,
                "job3",
                "annotate -h",
            ),
        ],
    )
    def test_merge_or_concat_raises_on_step_failure(
        self,
        sample_names_map,
        non_overlapping,
        identifier,
        failing_prefix,
    ):
        """Test that merge_or_concat_bcftools_temp_files raises BcftoolsCommandError when concat/annotate steps return non-zero exit code."""
        manager = BcftoolsQueryManager()
        manager.temp_files = []

        with (
            patch.object(manager, "_get_all_sample_names_from_vcf_files", return_value=sample_names_map),
            patch.object(manager, "_check_non_overlapping_sample_names", return_value=non_overlapping),
            patch.object(manager, "_prepare_txt_with_divbase_header_for_vcf", return_value=None),
            patch.object(
                manager, "_log_file_size", return_value=None
            ),  # Patch a None return to skip actual file size logging for this test
            patch("divbase_api.services.queries.os.rename", return_value=None),
            patch.object(
                manager,
                "run_bcftools",
                side_effect=lambda command: (
                    self.DummyProc(1)
                    if command.startswith(failing_prefix)
                    and (len(command) == len(failing_prefix) or command[len(failing_prefix)] == " ")
                    else self.DummyProc(0)
                ),
            ),
            pytest.raises(BcftoolsCommandError, match="Process exited with code 1") as exc_info,
        ):
            manager.merge_or_concat_bcftools_temp_files(list(sample_names_map.keys()), identifier=identifier)

        assert failing_prefix in str(exc_info.value)
