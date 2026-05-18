"""Unit tests for VCF queries modes in tasks.py"""

from unittest.mock import MagicMock, patch

import pytest

import divbase_api.services.vcf_queries as queries_module
from divbase_api.services.vcf_queries import (
    BCFToolsCommandConfig,
    BCFToolsInput,
    BcftoolsQueryManager,
    SampleFileMapping,
    ensure_csi_index,
    extract_region_scaffolds_from_command,
    validate_user_submitted_bcftools_command,
)
from divbase_api.worker.crud_dimensions import ProjectVCFDimensionsData, ProjectVCFDimensionsEntry
from divbase_api.worker.tasks import (
    VCFQuerySampleSelectionMode,
    _calculate_pairwise_overlap_types_for_sample_sets,
    _check_if_samples_can_be_combined_with_bcftools,
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
            "view -i 'FILTER=\"A;B\"'; view -r 1:1-1000",
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
            ("view -", "Do not use stdin '-' in '--command'"),
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

        def fake_run_bcftools(command: str, capture_output: bool = False, capture_stderr: bool = False):
            executed_commands.append(command)
            return DummyProc()

        monkeypatch.setattr(queries_module, "run_bcftools", fake_run_bcftools)
        monkeypatch.setattr(queries_module, "ensure_csi_index", lambda _file_path: None)
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

    def test_ensure_csi_index_raises_on_non_zero_return_code(self, tmp_path):
        """Test that ensure_csi_index raises BcftoolsCommandError with appropriate message when bcftools index command returns non-zero exit code."""
        vcf_path = tmp_path / "input.vcf.gz"  # .csi absent — triggers indexing
        index_proc = MagicMock()
        index_proc.returncode = 1
        index_proc.communicate.return_value = ("", "")
        with (
            patch("divbase_api.services.bcftools_helpers.run_bcftools", return_value=index_proc) as mock_run_bcftools,
            pytest.raises(BcftoolsCommandError, match="Process exited with code 1") as exc_info,
        ):
            ensure_csi_index(vcf_path)

        mock_run_bcftools.assert_called_once_with(command=f"index -f {vcf_path}", capture_stderr=True)
        assert str(vcf_path) in str(exc_info.value)

    def test_regression_ensure_csi_index_raises_task_user_error_on_unsorted_positions(self, tmp_path):
        """
        Regression test (negative outcome): unsorted VCF coordinates must raise a user-facing TaskUserError.
        Why: unsorted files cannot be indexed and would break bcftools orchestration with opaque errors otherwise.
        Reference: docs/development/bcftools_task_constraints.md ("VCF files need to be sorted by position").
        """
        vcf_path = tmp_path / "input.vcf.gz"  # .csi absent — triggers indexing
        unsorted_stderr = (
            "[E::hts_idx_push] Unsorted positions on sequence #1: 22053057 followed by 17504018\n"
            'index: failed to create index for "input.vcf.gz"\n'
        )
        index_proc = MagicMock()
        index_proc.returncode = 1
        index_proc.communicate.return_value = ("", unsorted_stderr)

        with (
            patch("divbase_api.services.bcftools_helpers.run_bcftools", return_value=index_proc) as mock_run_bcftools,
            pytest.raises(TaskUserError) as exc_info,
        ):
            ensure_csi_index(vcf_path)

        mock_run_bcftools.assert_called_once_with(command=f"index -f {vcf_path}", capture_stderr=True)
        msg = str(exc_info.value)
        assert "not sorted by position" in msg
        assert "bcftools sort" in msg
        assert "input.vcf.gz" in msg

    def test_ensure_csi_index_skips_indexing_when_index_already_exists(self, tmp_path):
        """Test that ensure_csi_index does not call run_bcftools when a .csi index already exists."""
        vcf_path = tmp_path / "test.vcf.gz"
        (tmp_path / "test.vcf.gz.csi").touch()  # create index so it exists
        with patch("divbase_api.services.bcftools_helpers.run_bcftools") as mock_run_bcftools:
            ensure_csi_index(vcf_path)

        mock_run_bcftools.assert_not_called()

    @pytest.mark.parametrize("filename", ["test.vcf.gz", "test.vcf"])
    def test_ensure_csi_index_uses_full_filename_suffix_in_command(self, tmp_path, filename):
        """
        Test that the .csi index command includes the full filename suffix.
        e.g. test.vcf.gz → index -f <path>/test.vcf.gz (not index -f <path>/test.vcf).
        """
        vcf_path = tmp_path / filename  # .csi absent — triggers indexing
        index_proc = MagicMock()
        index_proc.returncode = 0
        index_proc.communicate.return_value = ("", "")

        with patch("divbase_api.services.bcftools_helpers.run_bcftools", return_value=index_proc) as mock_run_bcftools:
            ensure_csi_index(vcf_path)

        mock_run_bcftools.assert_called_once_with(command=f"index -f {vcf_path}", capture_stderr=True)

    @pytest.mark.parametrize(
        "sample_names_map,non_overlapping,identifier,failing_prefix",
        [
            (
                {"file1.bcf": ["S1"], "file2.bcf": ["S2"]},
                True,
                "job1",
                "merge -Ou",
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
            patch("divbase_api.services.vcf_queries.os.rename", return_value=None),
            patch(
                "divbase_api.services.vcf_queries.run_bcftools",
                side_effect=lambda command, capture_output=False, capture_stderr=False: (
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


class TestSampleSetOverlapHelpers:
    @pytest.mark.parametrize(
        "sample_sets_dict,expected_identical_diff_order,expected_partly,expected_non_overlap",
        [
            (
                {
                    tuple(["S1", "S2", "S3"]): [],
                    tuple(["S3", "S4"]): [],
                    tuple(["S5", "S6"]): [],
                },
                [],
                [(tuple(["S1", "S2", "S3"]), tuple(["S3", "S4"]))],
                [(tuple(["S1", "S2", "S3"]), tuple(["S5", "S6"]))],
            ),
            (
                {
                    tuple(["A"]): [],
                    tuple(["B"]): [],
                    tuple(["C"]): [],
                },
                [],
                [],
                [
                    (tuple(["A"]), tuple(["B"])),
                    (tuple(["A"]), tuple(["C"])),
                    (tuple(["B"]), tuple(["C"])),
                ],
            ),
            (
                {
                    tuple(["X", "Y"]): [],
                    tuple(["Y", "X"]): [],
                },
                [(tuple(["X", "Y"]), tuple(["Y", "X"]))],
                [],
                [],
            ),
        ],
    )
    def test_calculate_pairwise_overlap_types_for_sample_sets(
        self,
        sample_sets_dict,
        expected_identical_diff_order,
        expected_partly,
        expected_non_overlap,
    ):
        """Test that _calculate_pairwise_overlap_types_for_sample_sets correctly calculates the pairwise overlap types for sample sets."""

        result = _calculate_pairwise_overlap_types_for_sample_sets(sample_sets_dict)

        assert isinstance(result.identical_elements_different_order, list)
        assert isinstance(result.partly_overlapping, list)
        assert isinstance(result.non_overlapping, list)

        for expected in expected_identical_diff_order:
            assert any(
                (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
                for pair in result.identical_elements_different_order
            )

        for expected in expected_partly:
            assert any(
                (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
                for pair in result.partly_overlapping
            )

        for expected in expected_non_overlap:
            assert any(
                (expected[0], expected[1]) == pair or (expected[1], expected[0]) == pair
                for pair in result.non_overlapping
            )

    @pytest.mark.parametrize(
        "files_to_download,dimensions_index,should_raise_error,expected_message_part",
        [
            (
                ["file1.vcf.gz", "file2.vcf.gz"],
                {
                    "vcf_files": [
                        {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A", "B"]},
                        {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B", "A"]},
                    ]
                },
                True,
                "identical elements but different order",
            ),
            (
                ["file1.vcf.gz", "file2.vcf.gz"],
                {
                    "vcf_files": [
                        {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A", "B"]},
                        {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B", "C"]},
                    ]
                },
                True,
                "partly overlapping",
            ),
            (
                ["file1.vcf.gz", "file2.vcf.gz"],
                {
                    "vcf_files": [
                        {"vcf_file_s3_key": "file1.vcf.gz", "samples": ["A"]},
                        {"vcf_file_s3_key": "file2.vcf.gz", "samples": ["B"]},
                    ]
                },
                False,
                "No unsupported sample sets found. Proceeding with bcftools pipeline.",
            ),
        ],
    )
    def test_check_if_samples_can_be_combined_with_bcftools_param(
        self,
        files_to_download,
        dimensions_index,
        should_raise_error,
        expected_message_part,
        caplog,
    ):
        """Test that check_if_samples_can_be_combined_with_bcftools raises TaskUserError when samples have incompatible overlaps."""
        vcf_dimensions_data = ProjectVCFDimensionsData(
            project_id=1,
            vcf_file_count=len(dimensions_index["vcf_files"]),
            vcf_files=[
                ProjectVCFDimensionsEntry(
                    vcf_file_s3_key=entry["vcf_file_s3_key"],
                    s3_version_id=None,
                    samples=entry.get("samples", []),
                )
                for entry in dimensions_index["vcf_files"]
            ],
        )

        if should_raise_error:
            with pytest.raises(TaskUserError) as excinfo:
                _check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_data)
            assert expected_message_part in str(excinfo.value)
        else:
            with caplog.at_level("INFO"):
                _check_if_samples_can_be_combined_with_bcftools(files_to_download, vcf_dimensions_data)
            assert expected_message_part in caplog.text
