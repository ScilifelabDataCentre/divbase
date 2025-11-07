"""
VCF dimensions (= technical metadata) DB Model.
"""

from typing import TYPE_CHECKING

from sqlalchemy import BigInteger, ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel

if TYPE_CHECKING:
    from divbase_api.models.projects import ProjectDB


class VCFMetadataDB(BaseDBModel):
    """
    DB model for the technical metadata ("dimensions") for the VCF files in the S3 buckets.
    One entry per VCF file. This table is updated by workers and read by the API.

    id (primary key), created_at and updated_at are inherited from BaseDBModel.

    To allow multiple projects to have VCF files with the same S3 key (filename),
    this model uses a UniqueConstraint (vcf_file_s3_key, project_id).
    """

    __tablename__ = "vcf_metadata"

    vcf_file_s3_key: Mapped[str] = mapped_column(String, index=True)
    project_id: Mapped[int] = mapped_column(
        ForeignKey("project.id", ondelete="CASCADE"),
        index=True,
    )
    s3_version_id: Mapped[str] = mapped_column(String, nullable=False, index=True)
    file_size_bytes: Mapped[int] = mapped_column(BigInteger)
    samples: Mapped[list[str]] = mapped_column(ARRAY(String), index=True)  # Sample names
    scaffolds: Mapped[list[str]] = mapped_column(ARRAY(String))  # Scaffold/chromosome names
    variant_count: Mapped[int] = mapped_column(BigInteger)
    sample_count: Mapped[int] = mapped_column(Integer)

    __table_args__ = (UniqueConstraint("vcf_file_s3_key", "project_id", name="unique_vcf_per_project"),)

    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="vcf_metadata")


class SkippedVCFDB(BaseDBModel):
    """
    DB model to track DivBase-generated result VCFs.
    These should ignored by queries since they contain duplicated data compared to the source VCFs.
    """

    __tablename__ = "skipped_vcf_files"

    vcf_file_s3_key: Mapped[str] = mapped_column(String, index=True)
    project_id: Mapped[int] = mapped_column(
        ForeignKey("project.id", ondelete="CASCADE"),
        index=True,
    )
    s3_version_id: Mapped[str] = mapped_column(String, nullable=False, index=True)
    skip_reason: Mapped[str | None] = mapped_column(String, nullable=True)

    __table_args__ = (UniqueConstraint("vcf_file_s3_key", "project_id", name="unique_skipped_vcf_per_project"),)

    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="skipped_vcf_files")
