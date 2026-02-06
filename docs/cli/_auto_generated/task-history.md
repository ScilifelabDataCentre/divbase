# `divbase-cli task-history`

Get the task history of query jobs submitted by the user to the DivBase API.

**Usage**:

```console
$ divbase-cli task-history [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `user`: Check status of all tasks submitted by the...
* `id`: Check status of a specific task submitted...
* `project`: Check status of all tasks submitted for a...

## `divbase-cli task-history user`

Check status of all tasks submitted by the user. Displays the latest 10 tasks by default, unless --limit is specified. Can be filtered by project name.

**Usage**:

```console
$ divbase-cli task-history user [OPTIONS]
```

**Options**:

* `--limit INTEGER`: Maximum number of tasks to display in the terminal. Sorted by recency.  [default: 10]
* `--project TEXT`: Optional project name to filter the user&#x27;s task history by project.
* `--help`: Show this message and exit.

## `divbase-cli task-history id`

Check status of a specific task submitted by the user by its task ID.

**Usage**:

```console
$ divbase-cli task-history id [OPTIONS] TASK_ID
```

**Arguments**:

* `TASK_ID`: Task ID to check the status of a specific query job.  [required]

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli task-history project`

Check status of all tasks submitted for a project. Requires a manager role in the project. Displays the latest 10 tasks by default, unless --limit is specified.

**Usage**:

```console
$ divbase-cli task-history project [OPTIONS] PROJECT
```

**Arguments**:

* `PROJECT`: Project name to check the task history for.  [required]

**Options**:

* `--limit INTEGER`: Maximum number of tasks to display in the terminal. Sorted by recency.  [default: 10]
* `--help`: Show this message and exit.
