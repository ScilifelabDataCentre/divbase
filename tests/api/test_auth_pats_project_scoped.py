"""
Tests using Personal Access Tokens (PATs) that are scoped to specific projects in DivBase API.

We are testing permissions/authorization logic here, not the task history data returned (covered in e2e_integration tests).
"""

from divbase_api.models.projects import ProjectRoles
from tests.api.helpers import add_test_project_member, bearer_header, create_pat, create_test_project, create_test_user

QUERY_ROUTE = "/api/v1/task-history/tasks/user/projects/{name}"
MANAGE_ROUTE = "/api/v1/task-history/projects/{name}"


async def test_unscoped_pat_uses_actual_membership_role(divbase_client, db_session):
    """all_projects=True should just mean use the real user's roles."""
    user = await create_test_user(db_session)
    project = await create_test_project(db_session)
    await add_test_project_member(db_session, project_id=project.id, user_email=user.email, role=ProjectRoles.QUERY)
    raw_pat = await create_pat(db_session, user_id=user.id, permissions={"all_projects": True})

    response = await divbase_client.get(QUERY_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))
    assert response.status_code == 200


async def test_pat_with_no_project_scope_is_denied_access_outright(divbase_client, db_session):
    """A scoped PAT (all_projects=False, default) whose `projects` dict doesn't include this project is denied."""
    user = await create_test_user(db_session)
    project = await create_test_project(db_session)
    await add_test_project_member(db_session, project_id=project.id, user_email=user.email, role=ProjectRoles.MANAGE)
    raw_pat = await create_pat(db_session, user_id=user.id)

    response = await divbase_client.get(QUERY_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))

    assert response.status_code == 403
    assert "does not have access to this project" in response.json()["detail"]


async def test_pat_scoped_to_a_different_project_is_denied_access(divbase_client, db_session):
    """A PAT scoped to project A cannot be used to access project B, even with a sufficient role."""
    user = await create_test_user(db_session)
    project_a = await create_test_project(db_session, name="project-a")
    project_b = await create_test_project(db_session, name="project-b")
    await add_test_project_member(db_session, project_id=project_a.id, user_email=user.email, role=ProjectRoles.MANAGE)
    await add_test_project_member(db_session, project_id=project_b.id, user_email=user.email, role=ProjectRoles.MANAGE)
    raw_pat = await create_pat(db_session, user_id=user.id, permissions={"projects": {str(project_a.id): "manage"}})

    response = await divbase_client.get(QUERY_ROUTE.format(name=project_b.name), headers=bearer_header(raw_pat))
    assert response.status_code == 403
    assert "does not have access to this project" in response.json()["detail"]


async def test_pat_role_downgrades_user_with_higher_membership_role(divbase_client, db_session):
    """PAT restricts to QUERY even though the user's real membership role is MANAGE - effective role is QUERY."""
    user = await create_test_user(db_session)
    project = await create_test_project(db_session)
    await add_test_project_member(db_session, project_id=project.id, user_email=user.email, role=ProjectRoles.MANAGE)
    raw_pat = await create_pat(db_session, user_id=user.id, permissions={"projects": {str(project.id): "query"}})

    query_response = await divbase_client.get(QUERY_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))
    assert query_response.status_code == 200

    manage_response = await divbase_client.get(MANAGE_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))
    assert manage_response.status_code == 403


async def test_pat_role_cannot_exceed_users_actual_membership_role(divbase_client, db_session):
    """PAT claims MANAGE, but the user's real membership role is only QUERY - effective role is capped at QUERY."""
    user = await create_test_user(db_session)
    project = await create_test_project(db_session)
    await add_test_project_member(db_session, project_id=project.id, user_email=user.email, role=ProjectRoles.QUERY)
    raw_pat = await create_pat(db_session, user_id=user.id, permissions={"projects": {str(project.id): "manage"}})

    manage_response = await divbase_client.get(MANAGE_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))
    assert manage_response.status_code == 403

    query_response = await divbase_client.get(QUERY_ROUTE.format(name=project.name), headers=bearer_header(raw_pat))
    assert query_response.status_code == 200
