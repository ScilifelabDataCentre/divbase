"""Unit tests for VCF queries modes in tasks.py"""

import pytest

from divbase_api.worker.tasks import (
    VCFQuerySampleSelectionMode,
    _determine_sample_selection_mode,
    _resolve_inputs_for_all_samples_mode,
    _resolve_inputs_for_cli_samples_mode,
    validate_user_submitted_bcftools_command,
)
from divbase_lib.exceptions import TaskUserError


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
            "view",
            "view --samples S1,S2 -r 21:15000000-25000000",
        ],
    )
    def test_validate_user_submitted_bcftools_command_accepts_valid_commands(self, command):
        """Test that valid user-submitted bcftools commands pass command validation."""
        result = validate_user_submitted_bcftools_command(command)
        assert result is None

    @pytest.mark.parametrize(
        "command,expected_msg",
        [
            ("merge -m none", "Unsupported bcftools command 'merge'"),
            ("view -S samples.txt", "-S/--samples-file"),
            ("view -Ssamples.txt", "-S/--samples-file"),
            ("view --samples-file samples.txt", "-S/--samples-file"),
            ("view --samples-file=samples.txt", "-S/--samples-file"),
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_invalid_commands(self, command, expected_msg):
        """Test that invalid user-submitted bcftools commands raise TaskUserError with expected message."""
        with pytest.raises(TaskUserError, match=expected_msg):
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

    @pytest.mark.parametrize(
        "command",
        [
            "view",
        ],
    )
    def test_validate_user_submitted_bcftools_command_requires_non_sample_option_for_all_samples(self, command):
        """Test that user-submitted bcftools command in all-samples mode that does not include at least one non-sample-selection view option raises TaskUserError."""
        with pytest.raises(TaskUserError, match="When using all-samples mode"):
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
        assert result is None

    @pytest.mark.parametrize(
        "command",
        [
            "view -s S1,S2 -r 21:15000000-25000000",
            "view -r 21:15000000-25000000; view --samples=S1,S2",
        ],
    )
    def test_validate_user_submitted_bcftools_command_rejects_samples_option_in_all_samples_mode(self, command):
        """Test that user-submitted bcftools command in all-samples mode that includes -s/--samples option raises TaskUserError."""
        with pytest.raises(TaskUserError, match="-s/--samples"):
            validate_user_submitted_bcftools_command(command, all_samples=True)


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
