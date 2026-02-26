"""
Unit tests for the ClientSideClientSideMetadataTSVValidator class.

Tests the shared validation logic from SharedMetadataValidator as exercised
through the CLI's ClientSideClientSideMetadataTSVValidator wrapper.

Shared fixtures are defined in tests/unit/conftest.py.
"""

from pathlib import Path

from divbase_cli.services.sample_metadata_tsv_validator import ClientSideClientSideMetadataTSVValidator


class TestValidTSV:
    """Test that a valid TSV with list notation passes all checks."""

    def test_valid_list_tsv_passes_all_checks(self, valid_list_tsv, project_samples):
        """Test that a valid TSV with list notation passes with no errors or warnings."""
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()

        assert len(errors) == 0
        assert len(warnings) == 0
        assert stats["total_columns"] == 4
        assert stats["user_defined_columns"] == 3
        assert stats["samples_in_tsv"] == 5
        assert stats["samples_matching_project"] == 5
        assert stats["has_multi_values"] is True
        assert "Population" in stats["numeric_columns"]
        assert "Area" in stats["string_columns"]
        assert "Weight" in stats["numeric_columns"]


class TestHeaderValidation:
    """Test validation of header row."""

    def test_wrong_first_column_name(self, header_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(header_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("First column must be named '#Sample_ID'" in e for e in errors)

    def test_duplicate_column_names(self, header_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(header_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Duplicate column names" in e and "Area" in e for e in errors)

    def test_duplicate_column_names_after_stripping_hash(self, tmp_path, project_samples):
        content = "#Sample_ID\tSample_ID\tPopulation\nS1\tS1_dup\t1\nS2\tS2_dup\t2\n"
        tsv_file = tmp_path / "dup_sample_id_cols.tsv"
        tsv_file.write_text(content)
        validator = ClientSideClientSideMetadataTSVValidator(tsv_file, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Duplicate column names" in e and "Sample_ID" in e for e in errors)

    def test_empty_column_name(self, header_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(header_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Empty column name" in e for e in errors)


class TestSampleIDValidation:
    """Test validation of Sample_ID column."""

    def test_empty_sample_id(self, sample_id_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(sample_id_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Sample_ID is empty" in e for e in errors)

    def test_list_value_in_sample_id(self, sample_id_errors_tsv, project_samples):
        """List notation in Sample_ID should produce an error."""
        validator = ClientSideClientSideMetadataTSVValidator(sample_id_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("list values" in e.lower() for e in errors)

    def test_duplicate_sample_id(self, sample_id_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(sample_id_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Duplicate Sample_ID" in e and "S1" in e for e in errors)


class TestFormattingValidation:
    """Test validation of TSV formatting."""

    def test_wrong_column_count(self, format_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(format_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Expected 3 tab-separated columns" in e and "found 2" in e for e in errors)

    def test_comma_in_cell_produces_warning(self, format_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(format_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("comma" in w.lower() for w in warnings)
        assert not any("comma" in e.lower() for e in errors)

    def test_whitespace_warning(self, format_errors_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(format_errors_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("leading or trailing whitespace" in w for w in warnings)

    def test_semicolons_in_cells_produce_warning(self, semicolons_in_cells_tsv):
        validator = ClientSideClientSideMetadataTSVValidator(semicolons_in_cells_tsv, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        assert any("semicolon" in w.lower() for w in warnings)
        assert not any("semicolon" in e.lower() for e in errors)


class TestListSyntaxValidation:
    """Test validation of Python list notation in cells."""

    def test_valid_list_syntax_no_errors(self, valid_list_tsv, project_samples):
        """Test that valid list notation like [1, 2] should parse without errors."""
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert len(errors) == 0

    def test_invalid_list_syntax_produces_errors(self, invalid_list_syntax_tsv, project_samples):
        """Test that cells starting with '[' but can't be parsed produce hard errors."""
        validator = ClientSideClientSideMetadataTSVValidator(invalid_list_syntax_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("[4" in e and "invalid" in e.lower() for e in errors)
        assert any("[1 2 3]" in e and "invalid" in e.lower() for e in errors)

    def test_invalid_list_syntax_collects_all_errors(self, invalid_list_syntax_tsv, project_samples):
        """Test that all invalid list cells should be reported, not just the first one."""
        validator = ClientSideClientSideMetadataTSVValidator(invalid_list_syntax_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        list_errors = [e for e in errors if "invalid" in e.lower() and "list" in e.lower()]
        assert len(list_errors) == 2

    def test_whitespace_insensitive_list_parsing(self, whitespace_variant_lists_tsv):
        """Test that [1,2,3], [1, 2, 3], and [ 1 , 2 , 3 ] all parse identically."""
        validator = ClientSideClientSideMetadataTSVValidator(whitespace_variant_lists_tsv, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        assert len(errors) == 0
        assert "Scores" in stats["numeric_columns"]

    def test_empty_list_parses_successfully(self, empty_list_tsv):
        """Test that an empty list [] should parse without errors."""
        validator = ClientSideClientSideMetadataTSVValidator(empty_list_tsv, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        list_errors = [e for e in errors if "list" in e.lower()]
        assert len(list_errors) == 0


class TestTypeValidation:
    """Test column type classification and mixed-type detection."""

    def test_numeric_list_column_is_numeric(self, numeric_list_tsv, project_samples):
        """Test that columns with only numeric list cells (multi-value) and numeric scalars (single-values) are numeric."""
        validator = ClientSideClientSideMetadataTSVValidator(numeric_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert "Scores" in stats["numeric_columns"]
        assert "Values" in stats["numeric_columns"]
        assert "Scores" not in stats["mixed_type_columns"]

    def test_string_list_column_is_string(self, string_list_tsv):
        """Test that columns with string list cells are classified as string."""
        validator = ClientSideClientSideMetadataTSVValidator(string_list_tsv, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        assert "Areas" in stats["string_columns"]
        assert "Score" in stats["numeric_columns"]

    def test_mixed_types_within_list_cell_is_error(self, mixed_type_list_cell_tsv):
        """Test that a list cell like [1, "two", 3] with mixed element types raise an error."""
        validator = ClientSideClientSideMetadataTSVValidator(mixed_type_list_cell_tsv, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        assert any("mixed element types" in e.lower() for e in errors)

    def test_mixed_types_across_cells_is_warning(self, mixed_type_across_cells_tsv, project_samples):
        """Test that a column with both numeric and string sends a warning (not error)."""
        validator = ClientSideClientSideMetadataTSVValidator(mixed_type_across_cells_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert "Population" in stats["mixed_type_columns"]
        assert not any("Population" in e and "mixed" in e.lower() for e in errors)

    def test_mixed_types_across_cells_produces_per_cell_warnings(self, mixed_type_across_cells_tsv, project_samples):
        """Test that mixed-type columns identifies which cells are outliers."""
        validator = ClientSideClientSideMetadataTSVValidator(mixed_type_across_cells_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        cell_warnings = [w for w in warnings if "non-numeric" in w.lower() and "Population" in w]
        assert len(cell_warnings) >= 1
        assert any("abc" in w for w in cell_warnings)
        assert any("def" in w for w in cell_warnings)

    def test_mixed_types_clarification_warning(self, mixed_type_across_cells_tsv, project_samples):
        """Test that a general clarification warning is sent for mixed-type columns."""
        validator = ClientSideClientSideMetadataTSVValidator(mixed_type_across_cells_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("clarification on mixed types" in w.lower() for w in warnings)

    def test_negative_numbers_are_numeric(self, numeric_list_tsv, project_samples):
        """Test that negative numbers are classified as numeric, not string."""
        validator = ClientSideClientSideMetadataTSVValidator(numeric_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert "Temperature" in stats["numeric_columns"]
        assert "Longitude" in stats["numeric_columns"]
        assert len(stats["mixed_type_columns"]) == 0

    def test_mixed_list_and_scalar_in_same_column(self, tmp_path):
        """Test that a column with both list cells and compatible scalar (single-value) cells pass."""
        content = "#Sample_ID\tScores\n"
        content += "S1\t[1, 2]\n"
        content += "S2\t5\n"
        content += "S3\t[3, 4]\n"
        p = tmp_path / "mixed_list_scalar.tsv"
        p.write_text(content)
        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2", "S3"})
        stats, errors, warnings = validator.validate()
        assert len(errors) == 0
        assert "Scores" in stats["numeric_columns"]


class TestDimensionMatching:
    """Test validation against project dimensions."""

    def test_samples_not_in_project(self, valid_list_tsv):
        project_samples = {"S1", "S2"}
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any(
            "following samples in the TSV were not found in the DivBase project's dimensions index" in e for e in errors
        )

    def test_samples_not_in_tsv(self, valid_list_tsv):
        project_samples = {"S1", "S2", "S3", "S10", "S20"}
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert any(
            "following samples in the DivBase project's dimensions index were not found in the TSV" in w and "S10" in w
            for w in warnings
        )

    def test_large_dimension_mismatch_is_summarized(self, valid_list_tsv):
        # Create samples named S0001, S0002, ..., S0050
        project_samples = {f"S{i:04d}" for i in range(1, 51)}
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        _, _, warnings = validator.validate()
        assert any("count: 50, showing first 20" in w for w in warnings)
        assert any("--full-sample-mismatch-names" in w for w in warnings)
        assert not any("S0050" in w for w in warnings)

    def test_large_dimension_mismatch_can_show_full_list(self, valid_list_tsv):
        project_samples = {f"S{i:04d}" for i in range(1, 51)}
        validator = ClientSideClientSideMetadataTSVValidator(
            valid_list_tsv, project_samples, dimensions_sample_preview_limit=None
        )
        _, _, warnings = validator.validate()
        assert any("count: 50, samples:" in w and "S0050" in w for w in warnings)
        assert not any("showing first 20" in w for w in warnings)
        assert not any("--full-sample-mismatch-names" in w for w in warnings)


class TestStatistics:
    """Test statistics collection."""

    def test_statistics_collection(self, valid_list_tsv, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert stats["total_columns"] == 4
        assert stats["user_defined_columns"] == 3
        assert stats["samples_in_tsv"] == 5
        assert stats["samples_matching_project"] == 5
        assert stats["total_project_samples"] == 5
        assert len(stats["numeric_columns"]) == 2
        assert len(stats["string_columns"]) == 1
        assert len(stats["mixed_type_columns"]) == 0
        assert stats["has_multi_values"] is True

    def test_no_multi_values_detected(self, no_multi_values_tsv):
        """Test that has_multi_values is False when no list cells exist."""
        validator = ClientSideClientSideMetadataTSVValidator(no_multi_values_tsv, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        assert stats["has_multi_values"] is False

    def test_multi_values_detected_via_list_cells(self, valid_list_tsv, project_samples):
        """Test that has_multi_values is True when list cells exist."""
        validator = ClientSideClientSideMetadataTSVValidator(valid_list_tsv, project_samples)
        stats, errors, warnings = validator.validate()
        assert stats["has_multi_values"] is True


class TestEdgeCases:
    """Test edge cases and error conditions."""

    def test_empty_file(self, project_samples, tmp_path):
        empty_file = tmp_path / "empty.tsv"
        empty_file.write_text("")
        validator = ClientSideClientSideMetadataTSVValidator(empty_file, project_samples)
        stats, errors, warnings = validator.validate()
        assert any("File is empty" in e for e in errors)

    def test_nonexistent_file(self, project_samples):
        validator = ClientSideClientSideMetadataTSVValidator(Path("/nonexistent/file.tsv"), project_samples)
        stats, errors, warnings = validator.validate()
        assert any("Failed to read file" in e for e in errors)

    def test_non_list_bracket_strings(self, tmp_path):
        """Test that strings like 'group[1]' that don't start with '[' should not trigger list parsing."""
        content = "#Sample_ID\tCode\n"
        content += "S1\tgroup[1]\n"
        content += "S2\tnormal\n"
        p = tmp_path / "bracket_strings.tsv"
        p.write_text(content)
        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        list_errors = [e for e in errors if "list" in e.lower()]
        assert len(list_errors) == 0

    def test_cell_starting_with_bracket_but_not_list(self, tmp_path):
        """Test that a cell like '[ref]' starts with '[' -- ast.literal_eval will fail, producing an error."""
        content = "#Sample_ID\tCode\n"
        content += "S1\t[ref]\n"
        content += "S2\tnormal\n"
        p = tmp_path / "bracket_ref.tsv"
        p.write_text(content)
        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        assert any("[ref]" in e and "invalid" in e.lower() for e in errors)

    def test_tuple_notation_is_not_a_list(self, tmp_path):
        """Test that a cell like '(1, 2)' that starts with '(' not '[' does not trigger list parsing."""
        content = "#Sample_ID\tCode\n"
        content += "S1\t(1, 2)\n"
        content += "S2\tnormal\n"
        p = tmp_path / "tuple.tsv"
        p.write_text(content)
        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        list_errors = [e for e in errors if "list" in e.lower()]
        assert len(list_errors) == 0


class TestQuotedCellValues:
    """Test quoted cell values that include commas/spaces and embedded quote characters."""

    def test_quoted_cell_with_comma_space_produces_comma_warning(self, tmp_path):
        content = '#Sample_ID\tArea\nS1\t"North, South"\nS2\tEast\n'
        p = tmp_path / "quoted_comma_space.tsv"
        p.write_text(content)

        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        assert len(errors) == 0
        assert any("comma" in w.lower() and "North, South" in w for w in warnings)

    def test_quoted_cell_with_embedded_double_quotes_is_allowed(self, tmp_path):
        content = '#Sample_ID\tComment\nS1\t"He said ""Hi"""\nS2\tPlain\n'
        p = tmp_path / "embedded_quotes.tsv"
        p.write_text(content)

        validator = ClientSideClientSideMetadataTSVValidator(p, {"S1", "S2"})
        stats, errors, warnings = validator.validate()
        assert len(errors) == 0
        assert not any("list syntax" in e.lower() for e in errors)
