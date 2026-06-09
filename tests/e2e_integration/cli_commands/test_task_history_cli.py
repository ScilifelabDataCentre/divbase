import datetime
import re
from contextlib import contextmanager
from csv import reader
from io import StringIO
from time import sleep
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from typer.testing import CliRunner

from divbase_cli.display_task_history import TaskHistoryDisplayManager
from divbase_cli.divbase_cli import app
from tests.e2e_integration.cli_commands.conftest import _create_logged_in_user_fixture

runner = CliRunner()

# NOTE: Getting the timing right for the below fixtures is not trivial. The test setup requires that
# multiple users log in and submit tasks before any tests are run.


@pytest.fixture(scope="module", autouse=True)
def all_users_tasks_submitted(
    edit_user_with_submitted_tasks,
    manage_user_with_submitted_tasks,
    edit_user_query_project_only_with_submitted_tasks,
    manage_user_query_project_only_with_submitted_tasks,
):
    """
    Ensure all users in the args of this fixture have logged in and submitted tasks before any test runs. It waits for tasks for the other fixtures to complete = have reached their yield statement.
    Together with the fixtures below, this results in two tasks being submitted per user. The tests in this module will then make tests based on these submitted tasks and which user they belong to.
    Since it is autouse, it will run before any tests in this module.
    """
    pass


@pytest.fixture(scope="module")
def edit_user_with_submitted_tasks(CONSTANTS):
    """
    Module-scoped fixture: login as edit user and submit tasks.

    factory is a generator object created, but not executed by _create_logged_in_user_fixture.
    The next() runs all code before the yield and then waits. This allows for tasks to be
    submitted for that logged in user. The final yield waits until all tests have been run,
    and the runs the clean up code after the yield in _create_logged_in_user_fixture.
    """
    factory = _create_logged_in_user_fixture("edit user")(CONSTANTS)
    next(factory)
    _submit_tasks_for_user(CONSTANTS)
    yield


@pytest.fixture(scope="module")
def manage_user_with_submitted_tasks(CONSTANTS):
    """
    Module-scoped fixture: login as manage user and submit tasks.

    factory is a generator object created, but not executed by _create_logged_in_user_fixture.
    The next() runs all code before the yield and then waits. This allows for tasks to be
    submitted for that logged in user. The final yield waits until all tests have been run,
    and the runs the clean up code after the yield in _create_logged_in_user_fixture.
    """
    factory = _create_logged_in_user_fixture("manage user")(CONSTANTS)
    next(factory)
    _submit_tasks_for_user(CONSTANTS)
    yield


@pytest.fixture(scope="module")
def edit_user_query_project_only_with_submitted_tasks(CONSTANTS):
    """
    Module-scoped fixture: login as edit user (query project only) and submit tasks.

    factory is a generator object created, but not executed by _create_logged_in_user_fixture.
    The next() runs all code before the yield and then waits. This allows for tasks to be
    submitted for that logged in user. The final yield waits until all tests have been run,
    and the runs the clean up code after the yield in _create_logged_in_user_fixture.
    """
    factory = _create_logged_in_user_fixture("edit user query-project only")(CONSTANTS)
    next(factory)
    _submit_tasks_for_user(CONSTANTS)
    yield


@pytest.fixture(scope="module")
def logged_in_admin_for_task_ids(CONSTANTS):
    """
    Module-scoped fixture: login as admin to query all task IDs. Not used to submit tasks, but to
    get the tasks IDs later on with the submitted_task_ids() fixture, since admins can see all tasks.
    """
    factory = _create_logged_in_user_fixture("admin")(CONSTANTS)
    next(factory)
    yield


@pytest.fixture(scope="module")
def manage_user_query_project_only_with_submitted_tasks(CONSTANTS):
    """
    Module-scoped fixture: login as manage user (query project only) and submit tasks.

    factory is a generator object created, but not executed by _create_logged_in_user_fixture.
    The next() runs all code before the yield and then waits. This allows for tasks to be
    submitted for that logged in user. The final yield waits until all tests have been run,
    and the runs the clean up code after the yield in _create_logged_in_user_fixture.
    """
    factory = _create_logged_in_user_fixture("manage user query-project only")(CONSTANTS)
    next(factory)
    _submit_tasks_for_user(CONSTANTS)
    yield


def _submit_tasks_for_user(CONSTANTS):
    """
    Helper function to submit tasks for a logged-in user.
    """
    query_string = "Area:West of Ireland,Northern Portugal;Sex:F"
    project_name = CONSTANTS["QUERY_PROJECT"]

    result_submit = runner.invoke(app, f"dimensions update --project {project_name}")
    assert result_submit.exit_code == 0
    sleep(0.5)

    result_submit = runner.invoke(app, f"query tsv '{query_string}' --project {project_name}")
    assert result_submit.exit_code == 0

    second_project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]
    result_submit = runner.invoke(
        app,
        f"query tsv '{query_string}' --project {second_project_name} --metadata-tsv-name sample_metadata_HOM_chr_split_version.tsv",
    )
    assert (
        result_submit.exit_code == 1
    )  # should fail due to missing dimensions index, this allows to test task history for fail messages


# TODO could make a login and submit task for admin, but admin user does not belong to any projects in the testing stack and can thus not submit tasks


@contextmanager
def capture_task_history_manager():
    """
    Context manager that captures the state of TaskHistoryDisplayManager instance so that
    assertions can be made on the actually content of the class and not the rich table output,
    which is sensitive to formatting errors (newlines, whitespaces, etc).
    """

    # TODO now that there is a --tsv option for divbase-cli task-history cmds, this fixture and the tests that call it could be refeactored to assert on the tsv output instead.
    captured_manager = None

    def capture_display_manager(self):
        nonlocal captured_manager
        captured_manager = self
        return original_print_task_history(self)

    original_print_task_history = TaskHistoryDisplayManager.print_task_history

    with patch.object(TaskHistoryDisplayManager, "print_task_history", capture_display_manager):
        yield lambda: captured_manager

    assert captured_manager is not None, "No TaskHistoryDisplayManager was captured"


@pytest.fixture(scope="module")
def submitted_task_ids(all_users_tasks_submitted, CONSTANTS, logged_in_admin_for_task_ids):
    """
    Fixture that provides task IDs for the tasks submitted by the test users in the setup fixtures.
    Intended for tests that use the `divbase task-history id` CLI command.
    Uses the admin test user to get all submitted tasks.
    """

    with capture_task_history_manager() as get_manager:
        runner.invoke(app, "task-history user")
        captured_manager = get_manager()

    task_ids_by_user = {}
    for task in captured_manager.task_items:
        user_email = task.submitter_email
        if user_email not in task_ids_by_user:
            task_ids_by_user[user_email] = []
        task_ids_by_user[user_email].append(task.id)
    return task_ids_by_user


def test_edit_user_can_only_see_their_own_task_history(CONSTANTS, logged_in_edit_user_with_existing_config):
    """Integration test where edit user can only see their own tasks."""

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, "task-history user")
        assert result_history.exit_code == 0

        captured_manager = get_manager()
    if not captured_manager:
        raise AssertionError("No TaskHistoryDisplayManager was captured")

    assert captured_manager.user_email == CONSTANTS["TEST_USERS"]["edit user"]["email"]

    user_emails = {task.submitter_email for task in captured_manager.task_items}

    assert CONSTANTS["TEST_USERS"]["edit user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user"]["email"] not in user_emails


def test_task_history_user_tsv_output_is_parseable_and_plain(
    logged_in_edit_user_with_existing_config,
):
    """Test that task-history user --tsv emits parseable TSV (not rich table output)."""
    result_history = runner.invoke(app, "task-history user --tsv --limit 50")
    assert result_history.exit_code == 0

    # Rich table border chars should not be present in TSV mode.
    for border_char in ("│", "┃", "┏", "┓", "┗", "┛", "┡", "┩", "┠", "┨", "┬", "┴", "┼", "┯", "┷", "━", "─"):
        assert border_char not in result_history.stdout

    rows = list(
        reader(StringIO(result_history.stdout), delimiter="\t")
    )  # use stringIO to handle newlines in stdout output
    assert len(rows) > 1  # Header + at least one data row

    expected_headers = [
        "Submitting user",
        "Task ID",
        "State",
        "Created at",
        "Started at",
        "Runtime (s)",
        "Result",
    ]
    assert rows[0] == expected_headers

    for row in rows[1:]:
        assert len(row) == len(expected_headers)

        task_id = row[1]
        state = row[2]
        result_text = row[6]

        assert task_id.isdigit()
        assert (
            "[" not in state and "]" not in state
        )  # Rich markup tags such as [green]...[/green] should not be present in TSV state cell
        assert state in {"QUEUING", "STARTED", "SUCCESS", "FAILURE", "RETRY", "REVOKED"}
        assert result_text != ""


def test_admin_user_can_see_all_task_history(CONSTANTS, logged_in_admin_with_existing_config):
    """Integration test where admin user can see all users' tasks."""

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, "task-history user")
        assert result_history.exit_code == 0

        captured_manager = get_manager()
    if not captured_manager:
        raise AssertionError("No TaskHistoryDisplayManager was captured")

    assert captured_manager.user_email == CONSTANTS["ADMIN_CREDENTIALS"]["email"]

    user_emails = {task.submitter_email for task in captured_manager.task_items}

    assert CONSTANTS["TEST_USERS"]["edit user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user"]["email"] in user_emails


def test_manage_user_can_see_all_task_history_for_a_project(CONSTANTS, logged_in_manage_user_with_existing_config):
    """Integration test where a manage user can see all tasks of a project they manage."""
    project_name = CONSTANTS["QUERY_PROJECT"]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, f"task-history project {project_name}")
        assert result_history.exit_code == 0

        captured_manager = get_manager()
    if not captured_manager:
        raise AssertionError("No TaskHistoryDisplayManager was captured")

    assert captured_manager.project_name == project_name

    user_emails = {task.submitter_email for task in captured_manager.task_items}

    assert CONSTANTS["TEST_USERS"]["edit user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user"]["email"] in user_emails


def test_task_history_project_uses_default_project_when_none_specified(
    CONSTANTS, logged_in_manage_user_with_existing_config
):
    """Test that task-history project falls back to the default project from config when --project is omitted."""
    default_project = CONSTANTS["DEFAULT_PROJECT"]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, "task-history project")
        assert result_history.exit_code == 0

        captured_manager = get_manager()

    assert captured_manager.project_name == default_project


def test_edit_user_can_filter_task_history_by_projects_they_belong_to(
    CONSTANTS, logged_in_edit_user_with_existing_config
):
    """Integration test where edit user can filter their own tasks by project."""
    project_name = CONSTANTS["QUERY_PROJECT"]
    omitted_project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, f"task-history user --project {project_name}")
        assert result_history.exit_code == 0

        captured_manager = get_manager()
    if not captured_manager:
        raise AssertionError("No TaskHistoryDisplayManager was captured")

    assert captured_manager.user_email == CONSTANTS["TEST_USERS"]["edit user"]["email"]
    assert captured_manager.project_name == project_name

    user_emails = {task.submitter_email for task in captured_manager.task_items}
    task_projects = {task.kwargs.project_name for task in captured_manager.task_items}

    assert CONSTANTS["TEST_USERS"]["edit user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user"]["email"] not in user_emails
    assert project_name in task_projects
    assert omitted_project_name not in task_projects


def test_read_user_cannot_see_task_history(CONSTANTS, logged_in_read_user_with_existing_config):
    """
    Integration test that read user cannot access the task history, since they cannot submit tasks.
    """
    project_name = CONSTANTS["QUERY_PROJECT"]
    result_history = runner.invoke(app, f"task-history user --project {project_name}")
    assert result_history.exit_code == 1
    assert "authorization_error" in str(result_history.exception)
    assert "Project not found or you don't have permission to view task history from this project" in str(
        result_history.exception
    )


def test_edit_user_cannot_see_task_history_for_project_not_member_of(
    CONSTANTS, logged_in_edit_user_query_project_only_with_existing_config
):
    """Integration test where edit user cannot see task history for a project they are not a member of."""
    non_member_project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    result_history = runner.invoke(app, f"task-history project {non_member_project_name}")
    assert result_history.exit_code == 1
    assert "project_not_found_error" in str(result_history.exception)


def test_task_history_id_tsv_can_lookup_submitted_query_task(
    CONSTANTS,
    logged_in_edit_user_with_existing_config,
    run_update_dimensions,
    project_map,
):
    """
    Smoke test that a task ID returned by `query vcf` can be looked up via
    `task-history id --tsv`.
    """
    project_name = CONSTANTS["QUERY_PROJECT"]
    project_id = project_map[project_name]
    bucket_name = CONSTANTS["PROJECT_TO_BUCKET_MAP"][project_name]
    user_id = 1
    run_update_dimensions(bucket_name=bucket_name, project_id=project_id, project_name=project_name, user_id=user_id)

    submit_command = (
        f"query vcf --tsv-filter 'Area:West of Ireland,Northern Portugal;' "
        f"--command 'view -r 21:15000000-25000000' --project {project_name}"
    )
    submit_result = runner.invoke(app, submit_command)
    assert submit_result.exit_code == 0
    task_id = submit_result.stdout.strip().split()[-1]

    max_retries = 10
    retry_delay = 0.5
    lookup_result = None
    rows = None

    for _ in range(max_retries):
        lookup_result = runner.invoke(app, f"task-history id {task_id} --tsv")
        if lookup_result.exit_code == 0:
            rows = list(reader(StringIO(lookup_result.stdout), delimiter="\t"))
            if len(rows) > 1:
                break
        sleep(retry_delay)

    assert lookup_result is not None
    assert lookup_result.exit_code == 0
    assert rows is not None and len(rows) > 1

    headers = rows[0]
    data_row = rows[1]
    assert "Task ID" in headers
    assert "State" in headers

    task_id_index = headers.index("Task ID")
    state_index = headers.index("State")
    assert data_row[task_id_index] == task_id
    assert data_row[state_index] in {"QUEUING", "STARTED", "SUCCESS", "FAILURE", "RETRY", "REVOKED"}


def test_manage_user_query_project_only_can_see_all_task_history_for_their_project(
    CONSTANTS, logged_in_manage_user_query_project_only_with_existing_config
):
    """Integration test where a manage user that only belongs to query-project can see all tasks of that project."""

    # Test that user can see all tasks for query-project
    project_name = CONSTANTS["QUERY_PROJECT"]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, f"task-history project {project_name}")
        assert result_history.exit_code == 0

        captured_manager = get_manager()

    assert captured_manager.project_name == project_name

    user_emails = {task.submitter_email for task in captured_manager.task_items}

    assert CONSTANTS["TEST_USERS"]["edit user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user"]["email"] in user_emails
    assert CONSTANTS["TEST_USERS"]["manage user query-project only"]["email"] in user_emails

    non_member_project_name = CONSTANTS["SPLIT_SCAFFOLD_PROJECT"]

    # Test that user cannot see tasks for a project they do not belong to
    result_history = runner.invoke(app, f"task-history project {non_member_project_name}")
    assert result_history.exit_code == 1
    assert "project_not_found_error" in str(result_history.exception)


def test_manage_user_can_get_task_id_from_project_even_when_they_did_not_submit_task(
    CONSTANTS,
    logged_in_manage_user_query_project_only_with_existing_config,
    submitted_task_ids,
):
    """
    Integration test where a manage user that only belongs to query-project can see task IDs of tasks
    submitted by other users in that project. Since task ID are UUID assigned at task submission, first
    get a list of all tasks submitted to the project, then extract a task submitted by an edit user. The
    Manage user should be able to access that task ID.
    """

    submitting_user = CONSTANTS["TEST_USERS"]["edit user query-project only"]["email"]
    edit_user_task_id = submitted_task_ids[submitting_user][0]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, f"task-history id {edit_user_task_id}")
        assert result_history.exit_code == 0

        captured_manager = get_manager()

    assert len(captured_manager.task_items) == 1
    assert captured_manager.task_items[0].id == edit_user_task_id


def test_edit_user_can_only_get_task_ids_they_submitted(
    CONSTANTS,
    logged_in_edit_user_query_project_only_with_existing_config,
    submitted_task_ids,
):
    """
    Integration test where an edit user that only belongs to query-project can only access tasks by ID for
    tasks that they submitted themselves.
    """
    submitting_user = CONSTANTS["TEST_USERS"]["edit user query-project only"]["email"]
    manage_user = CONSTANTS["TEST_USERS"]["manage user query-project only"]["email"]
    edit_user_task_id = submitted_task_ids[submitting_user][0]
    manage_user_task_id = submitted_task_ids[manage_user][0]

    with capture_task_history_manager() as get_manager:
        result_history = runner.invoke(app, f"task-history id {edit_user_task_id}")
        assert result_history.exit_code == 0

        captured_manager = get_manager()

    assert len(captured_manager.task_items) == 1
    assert captured_manager.task_items[0].id == edit_user_task_id
    task = next(t for t in captured_manager.task_items if t.id == edit_user_task_id)
    assert task.submitter_email == submitting_user

    result_history = runner.invoke(app, f"task-history id {manage_user_task_id}")
    assert result_history.exit_code == 1
    assert "authorization_error\nDetails" in str(result_history.exception)
    assert "Task ID not found or you don't have permission to view the history for this task ID." in str(
        result_history.exception
    )


def test_display_queuing_state_when_queue_full():
    """
    Unit test that TaskHistoryDisplayManager displays 'QUEUING' state for tasks
    that do not have a known Celery worker state (i.e., are still queued).
    """
    # Simulate a started task (should not be QUEUING)
    started_task = SimpleNamespace()
    started_task.id = 1
    started_task.status = "STARTED"
    started_task.submitter_email = "user2@example.com"
    started_task.created_at = datetime.datetime.now()
    started_task.started_at = datetime.datetime.now()
    started_task.runtime = 5.0
    started_task.result = None

    # Simulate a queued task (no status)
    queuing_task = SimpleNamespace()
    queuing_task.id = 2
    queuing_task.status = None
    queuing_task.submitter_email = "user@example.com"
    queuing_task.created_at = datetime.datetime.now()
    queuing_task.started_at = None
    queuing_task.runtime = None
    queuing_task.result = None

    # Simulate a queued task (PENDING status not in CELERY_STATES_EXCLUDING_PENDING)
    queuing_task_2 = SimpleNamespace()
    queuing_task_2.id = 3
    queuing_task_2.status = "PENDING"
    queuing_task_2.submitter_email = "user@example.com"
    queuing_task_2.created_at = datetime.datetime.now()
    queuing_task_2.started_at = None
    queuing_task_2.runtime = None
    queuing_task_2.result = None

    task_items = [started_task, queuing_task, queuing_task_2]

    manager = TaskHistoryDisplayManager(
        task_items=task_items,
        user_email="user@example.com",
        project_name=None,
        mode="user",
        display_limit=10,
    )

    with patch("rich.table.Table.add_row") as mock_add_row, patch("rich.console.Console.print"):
        manager.print_task_history()

    state_columns = [call_args[0][2] for call_args in mock_add_row.call_args_list]
    state_columns_clean = [re.sub(r"\[.*?\]", "", s).strip() for s in state_columns]

    assert "QUEUING" in state_columns_clean
    assert "STARTED" in state_columns_clean
