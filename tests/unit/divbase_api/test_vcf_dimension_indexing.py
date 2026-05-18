"""Unit tests for VCFDimensionCalculator in vcf_dimension_indexing.py."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy.exc import SQLAlchemyError

from divbase_api.services.bcftools_helpers import bgzip_vcf_for_indexing
from divbase_api.services.vcf_dimension_indexing import VCFDimensionCalculator, VCFDimensions
from divbase_api.worker.tasks import _remove_stale_dimensions_db_entries
from divbase_lib.exceptions import TaskUserError


@pytest.fixture
def calculator():
    return VCFDimensionCalculator()


class TestExtractSampleNamesFromVCFHeader:
    def test_normal_vcf_returns_sample_names(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator returns the sample names returned from #CHROM cols[9:]
        """
        vcf_header = (
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3\n"
        )
        mock_result = MagicMock()
        mock_result.stdout = vcf_header

        with patch("subprocess.run", return_value=mock_result):
            result = calculator._extract_sample_names_from_vcf_header(tmp_path / "test.vcf.gz")

        assert result == ["sample1", "sample2", "sample3"]

    def test_divbase_created_header_returns_none(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator returns None when the ##DivBase_created header is present.
        """
        vcf_header = (
            "##fileformat=VCFv4.2\n"
            '##DivBase_created="This is a results file created by a DivBase query; Date=Mon Oct 20 12:00:00 2025"\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        )
        mock_result = MagicMock()
        mock_result.stdout = vcf_header

        with patch("subprocess.run", return_value=mock_result):
            result = calculator._extract_sample_names_from_vcf_header(tmp_path / "divbase_result.vcf.gz")

        assert result is None

    def test_chrom_line_with_no_sample_cols_returns_empty_list(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator returns an empty list when the #CHROM line has no sample columns.
        """
        vcf_header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        mock_result = MagicMock()
        mock_result.stdout = vcf_header

        with patch("subprocess.run", return_value=mock_result):
            result = calculator._extract_sample_names_from_vcf_header(tmp_path / "sites_only.vcf")

        assert result == []

    def test_called_process_error_propagates(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator._extract_sample_names_from_vcf_header raises errors
        that are propagated to the task.py layer so that Celery properly marks the task as FAILURE.
        """
        with (
            patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "bcftools")),
            pytest.raises(subprocess.CalledProcessError),
        ):
            calculator._extract_sample_names_from_vcf_header(tmp_path / "bad.vcf.gz")

    def test_duplicate_sample_name_error_maps_to_task_user_error(self, calculator, tmp_path):
        """
        Test that duplicate sample IDs in the VCF header are mapped to a clear TaskUserError message.
        """
        with (
            patch(
                "subprocess.run",
                side_effect=subprocess.CalledProcessError(
                    1,
                    "bcftools",
                    stderr="[E::bcf_hdr_add_sample_len] Duplicated sample name 'NA00002'",
                ),
            ),
            pytest.raises(TaskUserError, match="contains duplicate sample IDs in the header"),
        ):
            calculator._extract_sample_names_from_vcf_header(tmp_path / "duplicate_header.vcf.gz")

    def test_called_process_error_with_stderr_maps_to_task_user_error(self, calculator, tmp_path):
        """
        Test that non-duplicate bcftools stderr from header parsing is surfaced as TaskUserError.
        """
        with (
            patch(
                "subprocess.run",
                side_effect=subprocess.CalledProcessError(
                    1,
                    "bcftools",
                    stderr="Failed to read from bad.vcf.gz: could not parse header",
                ),
            ),
            pytest.raises(TaskUserError, match="Details from bcftools"),
        ):
            calculator._extract_sample_names_from_vcf_header(tmp_path / "bad.vcf.gz")


class TestExtractScaffoldNamesAndVariantCount:
    def test_multi_scaffold_output_returns_correct_list_and_summed_count(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator correctly sums the variant count from the bcftools index --stats output.
        """
        # Three scaffolds with 100, 50, and 25 variants respectively
        stats_output = "chr1\t248956422\t100\nchr2\t242193529\t50\nchr3\t198295559\t25\n"
        mock_result = MagicMock()
        mock_result.stdout = stats_output

        with patch("subprocess.run", return_value=mock_result):
            scaffolds, count = calculator._extract_scaffold_names_and_variant_count_from_csi_index(
                tmp_path / "test.vcf.gz.csi"
            )

        assert scaffolds == ["chr1", "chr2", "chr3"]
        assert count == 175

    def test_single_scaffold_does_not_double_count(self, calculator, tmp_path):
        """
        Test that a single scaffold with a certain variant count returns that count without double-counting.
        """
        stats_output = "1\t22053058\t17\n"
        mock_result = MagicMock()
        mock_result.stdout = stats_output

        with patch("subprocess.run", return_value=mock_result):
            scaffolds, count = calculator._extract_scaffold_names_and_variant_count_from_csi_index(
                tmp_path / "test.vcf.gz.csi"
            )

        assert scaffolds == ["1"]
        assert count == 17

    def test_called_process_error_propagates(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator._extract_scaffold_names_and_variant_count_from_csi_index raises errors
        that are propagated to the task.py layer so that Celery properly marks the task as FAILURE.
        """
        with (
            patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "bcftools")),
            pytest.raises(subprocess.CalledProcessError),
        ):
            calculator._extract_scaffold_names_and_variant_count_from_csi_index(tmp_path / "test.vcf.gz.csi")

    def test_called_process_error_with_stderr_maps_to_task_user_error(self, calculator, tmp_path):
        """
        Test that bcftools stderr from index --stats is surfaced as TaskUserError instead of being suppressed.
        """
        with (
            patch(
                "subprocess.run",
                side_effect=subprocess.CalledProcessError(
                    1,
                    "bcftools",
                    stderr="index: failed to open test.vcf.gz",
                ),
            ),
            pytest.raises(TaskUserError, match="Details from bcftools"),
        ):
            calculator._extract_scaffold_names_and_variant_count_from_csi_index(tmp_path / "test.vcf.gz.csi")


class TestBgzipVCF:
    def test_bgzip_called_process_error_with_duplicate_sample_maps_to_task_user_error(self, calculator, tmp_path):
        """
        Test that duplicate-sample stderr during bgzip is mapped to a clear TaskUserError message.
        """
        proc = MagicMock()
        proc.returncode = 1
        proc.communicate.return_value = (
            "",
            "[E::bcf_hdr_add_sample_len] Duplicated sample name 'NA00002'",
        )
        with (
            patch("divbase_api.services.bcftools_helpers.run_bcftools", return_value=proc) as mock_run_bcftools,
            pytest.raises(TaskUserError, match="contains duplicate sample IDs in the header"),
        ):
            bgzip_vcf_for_indexing(input_vcf=tmp_path / "bad.vcf", output_vcf_gz=tmp_path / "bad.vcf.gz")
        mock_run_bcftools.assert_called_once_with(
            command=f"view -Oz -o {tmp_path / 'bad.vcf.gz'} {tmp_path / 'bad.vcf'}",
            capture_stderr=True,
        )

    def test_bgzip_called_process_error_with_stderr_maps_to_task_user_error(self, calculator, tmp_path):
        """
        Test that non-duplicate stderr during bgzip is surfaced as TaskUserError.
        """
        proc = MagicMock()
        proc.returncode = 1
        proc.communicate.return_value = (
            "",
            "Failed to read from bad.vcf: could not parse header",
        )
        with (
            patch("divbase_api.services.bcftools_helpers.run_bcftools", return_value=proc) as mock_run_bcftools,
            pytest.raises(TaskUserError, match="Details from bcftools"),
        ):
            bgzip_vcf_for_indexing(input_vcf=tmp_path / "bad.vcf", output_vcf_gz=tmp_path / "bad.vcf.gz")
        mock_run_bcftools.assert_called_once_with(
            command=f"view -Oz -o {tmp_path / 'bad.vcf.gz'} {tmp_path / 'bad.vcf'}",
            capture_stderr=True,
        )


class TestCalculateDimensions:
    def test_divbase_result_file_returns_none_without_indexing(self, calculator, tmp_path):
        """
        Test that if the VCF file is a DivBase result file (identified by the presence of the ##DivBase_created header),
        then calculate_dimensions returns None and does not attempt to index the file with CSI.
        """
        vcf_path = tmp_path / "divbase_result.vcf.gz"

        with (
            patch.object(calculator, "_extract_sample_names_from_vcf_header", return_value=None) as mock_header,
            patch("divbase_api.services.vcf_dimension_indexing.ensure_csi_index") as mock_index,
        ):
            result = calculator.calculate_dimensions(vcf_path)

        assert result is None
        mock_header.assert_called_once_with(vcf_path)
        mock_index.assert_not_called()

    def test_normal_vcf_returns_correct_vcf_dimensions_with_sorted_scaffolds(self, calculator, tmp_path):
        """
        Test that VCF files return a VCFDimensions object with the correct sample names, sample count, scaffold names, and variant count.
        """
        vcf_path = tmp_path / "test.vcf.gz"
        sample_names = ["sample_b", "sample_a"]
        unsorted_scaffold_names = ["chr3", "chr1", "chr2"]
        variant_count = 42

        with (
            patch.object(calculator, "_extract_sample_names_from_vcf_header", return_value=sample_names),
            patch("divbase_api.services.vcf_dimension_indexing.ensure_csi_index"),
            patch.object(
                calculator,
                "_extract_scaffold_names_and_variant_count_from_csi_index",
                return_value=(unsorted_scaffold_names, variant_count),
            ),
        ):
            result = calculator.calculate_dimensions(vcf_path)

        assert isinstance(result, VCFDimensions)
        assert result.sample_names == ["sample_b", "sample_a"]
        assert result.sample_count == 2
        assert result.scaffolds == ["chr1", "chr2", "chr3"]
        assert result.variants == 42


class TestStaleDimensionsCleanup:
    def test_remove_stale_dimensions_no_stale_entries_returns_empty_list(self):
        """Test that for the case of no stale entrie, helper function returns an empty list."""
        db = MagicMock()

        with (
            patch("divbase_api.worker.tasks.delete_vcf_metadata_batch") as mock_delete_indexed,
            patch("divbase_api.worker.tasks.delete_skipped_vcf_batch") as mock_delete_skipped,
        ):
            deleted = _remove_stale_dimensions_db_entries(
                indexed_vcf_keys={"a.vcf.gz"},
                skipped_vcf_keys={"b.vcf.gz"},
                current_vcf_files_in_bucket={"a.vcf.gz", "b.vcf.gz"},
                project_id=1,
                db=db,
            )

        assert deleted == []
        mock_delete_indexed.assert_not_called()
        mock_delete_skipped.assert_not_called()

    def test_remove_stale_dimensions_deletes_stale_indexed_and_skipped(self):
        """Test that stale indexed and skipped entries are both deleted by the helper function."""
        db = MagicMock()

        with (
            patch("divbase_api.worker.tasks.delete_vcf_metadata_batch") as mock_delete_indexed,
            patch("divbase_api.worker.tasks.delete_skipped_vcf_batch") as mock_delete_skipped,
        ):
            deleted = _remove_stale_dimensions_db_entries(
                indexed_vcf_keys={"keep.vcf.gz", "old_indexed.vcf.gz"},
                skipped_vcf_keys={"old_skipped.vcf.gz"},
                current_vcf_files_in_bucket={"keep.vcf.gz"},
                project_id=7,
                db=db,
            )

        assert set(deleted) == {"old_indexed.vcf.gz", "old_skipped.vcf.gz"}
        mock_delete_indexed.assert_called_once_with(db=db, vcf_file_s3_key_batch=["old_indexed.vcf.gz"], project_id=7)
        mock_delete_skipped.assert_called_once_with(db=db, vcf_file_s3_key_batch=["old_skipped.vcf.gz"], project_id=7)

    def test_remove_stale_dimensions_raises_task_user_error_when_indexed_delete_fails(self):
        """Test that indexed cleanup DB errors fail with TaskUserError."""
        db = MagicMock()

        with (
            patch(
                "divbase_api.worker.tasks.delete_vcf_metadata_batch",
                side_effect=SQLAlchemyError("db down"),
            ) as mock_delete_indexed,
            patch("divbase_api.worker.tasks.delete_skipped_vcf_batch") as mock_delete_skipped,
            pytest.raises(TaskUserError, match="Failed to clean up stale VCF dimensions entries"),
        ):
            _remove_stale_dimensions_db_entries(
                indexed_vcf_keys={"old_indexed.vcf.gz"},
                skipped_vcf_keys={"old_skipped.vcf.gz"},
                current_vcf_files_in_bucket=set(),
                project_id=99,
                db=db,
            )

        mock_delete_indexed.assert_called_once()
        mock_delete_skipped.assert_not_called()

    def test_remove_stale_dimensions_raises_task_user_error_when_skipped_delete_fails(self):
        """Test that skipped cleanup DB errors also fail with TaskUserError."""
        db = MagicMock()

        with (
            patch("divbase_api.worker.tasks.delete_vcf_metadata_batch") as mock_delete_indexed,
            patch(
                "divbase_api.worker.tasks.delete_skipped_vcf_batch",
                side_effect=SQLAlchemyError("db issue"),
            ) as mock_delete_skipped,
            pytest.raises(TaskUserError, match="Failed to clean up stale VCF dimensions entries"),
        ):
            _remove_stale_dimensions_db_entries(
                indexed_vcf_keys={"old_indexed.vcf.gz"},
                skipped_vcf_keys={"old_skipped.vcf.gz"},
                current_vcf_files_in_bucket=set(),
                project_id=42,
                db=db,
            )

        mock_delete_indexed.assert_called_once()
        mock_delete_skipped.assert_called_once()
