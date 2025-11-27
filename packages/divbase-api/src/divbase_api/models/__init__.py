"""
Import all models here to ensure proper initialization order.

They are imported in dependency order to avoid circular import issues.
"""

from divbase_api.models.base import Base, BaseDBModel
from divbase_api.models.projects import ProjectDB, ProjectMembershipDB, ProjectRoles
from divbase_api.models.revoked_tokens import RevokedTokenDB, TokenRevokeReason
from divbase_api.models.task_history import TaskHistoryDB
from divbase_api.models.users import UserDB
from divbase_api.models.vcf_dimensions import SkippedVCFDB, VCFMetadataDB

__all__ = [
    "Base",
    "BaseDBModel",
    "UserDB",
    "ProjectDB",
    "ProjectRoles",
    "ProjectMembershipDB",
    "VCFMetadataDB",
    "SkippedVCFDB",
    "TaskHistoryDB",
    "TokenRevokeReason",
    "RevokedTokenDB",
]
