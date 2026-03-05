"""Unit tests for VCFDimensionCalculator in vcf_dimension_indexing.py."""

import subprocess
from unittest.mock import MagicMock, patch

import pytest

from divbase_api.worker.vcf_dimension_indexing import VCFDimensionCalculator, VCFDimensions


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


class TestIndexVCFWithCSI:
    def test_success_returns_correct_csi_path(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator returns the correct CSI path.
        Including that .csi is appended to the full suffix, not just to .vcf
        e.g. test.vcf.gz → test.vcf.gz.csi (not test.vcf.csi).
        """
        vcf_path = tmp_path / "test.vcf.gz"
        expected_csi_path = tmp_path / "test.vcf.gz.csi"

        with patch("subprocess.run", return_value=MagicMock()):
            result = calculator._index_vcf_with_csi(vcf_path)

        assert result == expected_csi_path

    def test_plain_vcf_returns_correct_csi_path(self, calculator, tmp_path):
        """Plain uncompressed .vcf → .vcf.csi (not some other extension)."""
        vcf_path = tmp_path / "test.vcf"
        expected_csi_path = tmp_path / "test.vcf.csi"

        with patch("subprocess.run", return_value=MagicMock()):
            result = calculator._index_vcf_with_csi(vcf_path)

        assert result == expected_csi_path

    def test_called_process_error_propagates(self, calculator, tmp_path):
        """
        Test that the VCFDimensionCalculator._index_vcf_with_csi raises errors
        that are propagated to the task.py layer so that Celery properly marks the task as FAILURE.
        """
        with (
            patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "bcftools")),
            pytest.raises(subprocess.CalledProcessError),
        ):
            calculator._index_vcf_with_csi(tmp_path / "bad.vcf.gz")


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


class TestCalculateDimensions:
    def test_divbase_result_file_returns_none_without_indexing(self, calculator, tmp_path):
        """
        Test that if the VCF file is a DivBase result file (identified by the presence of the ##DivBase_created header),
        then calculate_dimensions returns None and does not attempt to index the file with CSI.
        """
        vcf_path = tmp_path / "divbase_result.vcf.gz"

        with (
            patch.object(calculator, "_extract_sample_names_from_vcf_header", return_value=None) as mock_header,
            patch.object(calculator, "_index_vcf_with_csi") as mock_index,
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
            patch.object(calculator, "_index_vcf_with_csi", return_value=tmp_path / "test.vcf.gz.csi"),
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
