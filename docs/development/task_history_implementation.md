# Task History Implementation

DivBase users can view the history of their submitted tasks using the `divbase-cli task-history` CLI command. This document describes how the task history implementation works. As long as a task has been implemented as decribed in [Celery Task Implementation](celery_task_implementation.md), there is no need to make any alterations to the Task History code.

There are three subcommands for `divbase-cli task-history`: `user`, `id`, and `project`. These are described in the subheadings below. All commands take the `--limit` option which controls the number of tasks to display descending from the most recent task (default: 10).

## task-history user

The `user` subcommand fetches all tasks submitted by the current logged in user from the postgreSQL database. It is possible to filter the tasks from a specific DivBase project in the user's task history with the `--project` option.

![Task History User Sequence Diagram](../assets/diagrams/task_history_user_sequence_diagram.svg)

Figure 1: Sequence diagram of signal flow for the `task-history user` command.

## task-history id

The `id` subcommand fetches a single task history entry based on its task ID (the integer returned to the user upon submitting the task; not the Celery task UUID that is only used internally) from the postgreSQL database. The user needs to have permission to view the task ID. A user with the MANAGE use roles can view all task IDs from the project(s) they manage, otherwise users need to be the submitting user of the task ID in order to view it.

![Task History ID Sequence Diagram](../assets/diagrams/task_history_id_sequence_diagram.svg)

Figure 2: Sequence diagram of signal flow for the `task-history id` command.

## task-history project

The `project` subcommand requires that the user has a MANAGE role in the project they want to view tasks from. If so, this command will fetch all tasks history entries for the project from the postgreSQL database.

![Task History Project Sequence Diagram](../assets/diagrams/task_history_project_sequence_diagram.svg)

Figure 3: Sequence diagram of signal flow for the `task-history project` command.
