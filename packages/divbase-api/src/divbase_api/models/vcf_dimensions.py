"""
VCF dimensions (= technical metadata) DB Model.
"""

from sqlalchemy import BigInteger, DateTime, ForeignKey, Integer, String, func
from sqlalchemy.dialects.postgresql import ARRAY
from sqlalchemy.orm import Mapped, mapped_column, relationship

from divbase_api.models.base import BaseDBModel
from divbase_api.models.projects import ProjectDB


class VCFMetadataDB(BaseDBModel):
    """
    DB model for the technical metadata ("dimensions") for the VCF files in the S3 buckets.
    One entry per VCF file. This table is updated by workers and read by the API.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "vcf_metadata"

    vcf_file_s3_key: Mapped[str] = mapped_column(String, primary_key=True, index=True)
    project_id: Mapped[int] = mapped_column(
        ForeignKey("project.id", ondelete="CASCADE"),  # Database-level cascade
        index=True,
    )
    s3_version_id: Mapped[str] = mapped_column(String, nullable=False, index=True)
    file_size_bytes: Mapped[int] = mapped_column(BigInteger)
    indexed_at: Mapped[DateTime] = mapped_column(DateTime, default=func.now())

    samples: Mapped[list[str]] = mapped_column(ARRAY(String), index=True)  # Sample names
    scaffolds: Mapped[list[str]] = mapped_column(ARRAY(String))  # Scaffold/chromosome names
    variant_count: Mapped[int] = mapped_column(BigInteger)
    sample_count: Mapped[int] = mapped_column(Integer)

    project: Mapped["ProjectDB"] = relationship("ProjectDB", back_populates="vcf_metadata")
