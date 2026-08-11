"""
API tests for routes/task_history.py,

The tests here are evaluating the permissions and authorization logic for the task history endpoints,
not the task history data returned (covered in e2e_integration tests).
"""

from divbase_api.models.projects import ProjectRoles
from tests.api.helpers import add_test_project_member, bearer_header, create_test_project, create_test_user, login


async def test_project_task_history_allowed_for_member_with_query_role(divbase_client, db_session):
    """Viewing your own task history for a project requires at least QUERY level, so READ access denied"""
    non_member_user = await create_test_user(db_session)
    read_user = await create_test_user(db_session)
    query_user = await create_test_user(db_session)
    project = await create_test_project(db_session)
    await add_test_project_member(db_session, project_id=project.id, user_email=read_user.email, role=ProjectRoles.READ)
    await add_test_project_member(
        db_session, project_id=project.id, user_email=query_user.email, role=ProjectRoles.QUERY
    )

    non_member_login_data = await login(divbase_client, non_member_user.email)
    response = await divbase_client.get(
        f"/api/v1/task-history/tasks/user/projects/{project.name}",
        headers=bearer_header(non_member_login_data.access_token),
    )
    assert response.status_code == 404
    assert response.json()["type"] == "ProjectNotFoundError"

    read_login_data = await login(divbase_client, read_user.email)
    response = await divbase_client.get(
        f"/api/v1/task-history/tasks/user/projects/{project.name}", headers=bearer_header(read_login_data.access_token)
    )
    assert response.status_code == 403

    query_login_data = await login(divbase_client, query_user.email)
    response = await divbase_client.get(
        f"/api/v1/task-history/tasks/user/projects/{project.name}", headers=bearer_header(query_login_data.access_token)
    )
    assert response.status_code == 200
    assert response.json() == []


async def test_project_tasks_endpoint_denied_for_member_without_manage_role(divbase_client, db_session):
    """This is the route that returns all tasks for a project, not just the current user's tasks. It requires MANAGE role."""
    query_user = await create_test_user(db_session)
    manage_user = await create_test_user(db_session)
    project = await create_test_project(db_session)

    await add_test_project_member(
        db_session, project_id=project.id, user_email=query_user.email, role=ProjectRoles.QUERY
    )
    await add_test_project_member(
        db_session, project_id=project.id, user_email=manage_user.email, role=ProjectRoles.MANAGE
    )
    query_login_data = await login(divbase_client, query_user.email)

    query_response = await divbase_client.get(
        f"/api/v1/task-history/projects/{project.name}", headers=bearer_header(query_login_data.access_token)
    )
    assert query_response.status_code == 403

    manager_login_data = await login(divbase_client, manage_user.email)
    manager_response = await divbase_client.get(
        f"/api/v1/task-history/projects/{project.name}", headers=bearer_header(manager_login_data.access_token)
    )
    assert manager_response.status_code == 200


async def test_project_tasks_endpoint_denied_for_admin_who_is_not_a_project_member(divbase_client, db_session):
    """An admin who is not a member of a project gets a 404, not 200."""
    admin = await create_test_user(db_session, is_admin=True)
    project = await create_test_project(db_session)
    login_data = await login(divbase_client, admin.email)

    response = await divbase_client.get(
        f"/api/v1/task-history/projects/{project.name}", headers=bearer_header(login_data.access_token)
    )
    assert response.status_code == 404
