"""Unit tests for VCF queries modes in tasks.py"""

import pytest

from divbase_api.worker.tasks import (
    VCFQuerySampleSelectionMode,
    _determine_sample_selection_mode,
    _resolve_inputs_for_all_samples_mode,
    _resolve_inputs_for_cli_samples_mode,
    _validate_user_submitted_bcftools_command,
)
from divbase_lib.exceptions import TaskUserError


class TestDetermineSampleSelectionMode:
    @pytest.mark.parametrize(
        "tsv_filter,samples,expected_mode",
        [
            ("Area:Northern Portugal", None, VCFQuerySampleSelectionMode.SAMPLE_METADATA_QUERY),
            (None, ["S1", "S2"], VCFQuerySampleSelectionMode.CLI_SAMPLES),
            (None, None, VCFQuerySampleSelectionMode.ALL_SAMPLES),
        ],
    )
    def test_determine_sample_selection_mode(self, tsv_filter, samples, expected_mode):
        """Test that sample selection mode is determined correctly from tsv_filter/samples inputs."""
        mode = _determine_sample_selection_mode(tsv_filter=tsv_filter, samples=samples)
        assert mode == expected_mode


class TestValidateUserSubmittedBcftoolsCommand:
    @pytest.mark.parametrize(
        "command",
        [
            "view",
            "view --samples S1,S2 -r 21:15000000-25000000",
        ],
    )
    def test_validate_user_submitted_bcftools_command_accepts_valid_commands(self, command):
        """Test that valid user-submitted bcftools commands pass command validation."""
        result = _validate_user_submitted_bcftools_command(command)
        assert result is None

    @pytest.mark.parametrize(
        "command,expected_msg",
        [
            ("merge -m none", "Unsupported bcftools command 'merge'"),
            ("view -S samples.txt", "Do not use bcftools sample-file options"),
            ("view -Ssamples.txt", "Do not use bcftools sample-file options"),
            ("view --samples-file samples.txt", "Do not use bcftools sample-file options"),
            ("view --samples-file=samples.txt", "Do not use bcftools sample-file options"),
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_invalid_commands(self, command, expected_msg):
        """Test that invalid user-submitted bcftools commands raise TaskUserError with expected message."""
        with pytest.raises(TaskUserError, match=expected_msg):
            _validate_user_submitted_bcftools_command(command)

    def test_validate_user_submitted_bcftools_command_rejects_unparseable_segment(self):
        """Test that unparseable user-submitted bcftools command segments raise TaskUserError."""
        with pytest.raises(TaskUserError, match="Could not parse --command segment at position 1"):
            _validate_user_submitted_bcftools_command('view -r "chr1:1-1000')


class TestResolveInputsForCliSamplesMode:
    def test_resolve_inputs_for_cli_samples_mode(self):
        """Test that CLI sample-mode inputs are resolved to files and sample mappings correctly."""
        vcf_dimensions_data = {
            "vcf_files": [
                {"vcf_file_s3_key": "fileA.vcf.gz", "samples": ["S1", "S2"]},
                {"vcf_file_s3_key": "fileB.vcf.gz", "samples": ["S3"]},
            ]
        }
        result = _resolve_inputs_for_cli_samples_mode(
            samples=["S3", "S1", "S3"], vcf_dimensions_data=vcf_dimensions_data
        )

        assert result.files_to_download == ["fileA.vcf.gz", "fileB.vcf.gz"]
        assert result.sample_and_filename_subset == [
            {"Sample_ID": "S1", "Filename": "fileA.vcf.gz"},
            {"Sample_ID": "S3", "Filename": "fileB.vcf.gz"},
        ]
        assert result.unique_sample_ids == ["S3", "S1"]
        assert result.metadata_path is None

    def test_resolve_inputs_for_cli_samples_mode_missing_samples_raises(self):
        """Test that CLI sample-mode input raises TaskUserError when requested samples are missing."""
        vcf_dimensions_data = {
            "vcf_files": [
                {"vcf_file_s3_key": "fileA.vcf.gz", "samples": ["S1", "S2"]},
            ]
        }
        with pytest.raises(TaskUserError, match="were not found in the project's dimensions index"):
            _resolve_inputs_for_cli_samples_mode(
                samples=["S1", "DOES_NOT_EXIST"], vcf_dimensions_data=vcf_dimensions_data
            )


class TestResolveInputsForAllSamplesMode:
    def test_resolve_inputs_for_all_samples_mode(self):
        """Test that all-samples mode resolves all files and sample mappings correctly."""
        vcf_dimensions_data = {
            "vcf_files": [
                {"vcf_file_s3_key": "fileA.vcf.gz", "samples": ["S1", "S2"]},
                {"vcf_file_s3_key": "fileB.vcf.gz", "samples": ["S3"]},
            ]
        }
        result = _resolve_inputs_for_all_samples_mode(vcf_dimensions_data=vcf_dimensions_data)

        assert result.files_to_download == ["fileA.vcf.gz", "fileB.vcf.gz"]
        assert result.sample_and_filename_subset == [
            {"Sample_ID": "S1", "Filename": "fileA.vcf.gz"},
            {"Sample_ID": "S2", "Filename": "fileA.vcf.gz"},
            {"Sample_ID": "S3", "Filename": "fileB.vcf.gz"},
        ]
        assert set(result.unique_sample_ids) == {"S1", "S2", "S3"}
        assert result.metadata_path is None
