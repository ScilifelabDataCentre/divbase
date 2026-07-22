"""Unit tests for PATPermissions in api_schemas/personal_access_tokens.py."""

import pytest
from pydantic import ValidationError

from divbase_api.models.projects import ProjectRoles
from divbase_lib.api_schemas.personal_access_tokens import ALLOWED_PROJECT_ROLES, PATPermissions


def test_allowed_project_roles_matches_project_roles_enum():
    """This test is to ensure the two different possible role definitions stay in sync"""
    lib_roles = set(ALLOWED_PROJECT_ROLES)
    api_roles = {role.value for role in ProjectRoles}
    assert lib_roles == api_roles, "Project role definitions from divbase-lib and divbase-api have drifted apart."


def test_default_pat_is_restrictive_by_default():
    """A PATPermissions with no arguments should grant access to nothing."""
    permissions = PATPermissions()
    assert permissions.all_projects is False
    assert permissions.projects == {}
    assert permissions.task_history is False


def test_each_valid_role_is_accepted():
    """Should not raise an error"""
    projects = {str(idx): role for idx, role in enumerate(ALLOWED_PROJECT_ROLES, start=1)}
    _ = PATPermissions(projects=projects)


def test_role_not_in_allow_list_is_rejected():
    # lowercase "read" would be valid, but uppercase "READ" is not
    for role in ["owner", "superadmin", "", "READ"]:
        with pytest.raises(ValidationError, match="invalid role"):
            PATPermissions(projects={"1": role})
