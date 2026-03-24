import pytest
from pydantic import ValidationError

from divbase_lib.api_schemas.queries import BcftoolsQueryKwargs, BcftoolsQueryRequest


def test_bcftools_query_request_rejects_empty_command():
    """Test that BcftoolsQueryRequest raises ValidationError when command is empty or only whitespace."""
    with pytest.raises(ValidationError, match="non-empty bcftools view string"):
        BcftoolsQueryRequest(
            command="   ",
            all_samples=True,
        )


def test_bcftools_query_kwargs_rejects_empty_command():
    """Test that BcftoolsQueryKwargs raises ValidationError when command is empty or only whitespace."""
    with pytest.raises(ValidationError, match="non-empty bcftools view string"):
        BcftoolsQueryKwargs(
            command="",
            bucket_name="bucket",
            project_id=1,
            project_name="project",
            user_id=1,
            job_id=1,
            all_samples=True,
        )


@pytest.mark.parametrize(
    "command,expected_match",
    [
        ("; view -r 21:1-1000", "empty pipeline segment at position 1"),
        ("view -r 21:1-1000;", "empty pipeline segment at position 2"),
        ("view -r 21:1-1000;;view -e 'QUAL<20'", "empty pipeline segment at position 2"),
        (" ; ", "empty pipeline segment at position 1"),
    ],
)
def test_bcftools_query_request_rejects_empty_pipeline_segments(command, expected_match):
    """Test that BcftoolsQueryRequest raises ValidationError when command includes empty pipeline segments separated by semicolons."""
    with pytest.raises(ValidationError, match=expected_match):
        BcftoolsQueryRequest(
            command=command,
            all_samples=True,
        )


@pytest.mark.parametrize(
    "command",
    [
        "view -r 21:1-1000; view -e 'QUAL<20'",
        "view -r 1,4,6,21,24; view -g hom; view -i 'MAF>0.05'",
        "view -A",
    ],
)
def test_bcftools_query_request_accepts_valid_pipeline_segments(command):
    """Test that BcftoolsQueryRequest accepts commands with valid pipeline segments separated by semicolons and does not raise ValidationError."""
    model = BcftoolsQueryRequest(
        command=command,
        all_samples=True,
    )
    assert model.command == command


def test_bcftools_query_request_rejects_semicolon_inside_quoted_segment_current_behavior():
    """Test that BcftoolsQueryRequest raises ValidationError when command includes semicolons inside quoted segments, which is currently not supported."""
    with pytest.raises(ValidationError, match="empty pipeline segment"):
        BcftoolsQueryRequest(
            command="view -i 'INFO/CSQ~\"A;B\";;QUAL>20'",
            all_samples=True,
        )
