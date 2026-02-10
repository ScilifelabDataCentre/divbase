# `divbase-cli version`

Add, view and remove versions representing the state of all files in the entire project at the current timestamp.

**Usage**:

```console
$ divbase-cli version [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `add`: Add a new project version entry which...
* `list`: List all entries in the project versioning...
* `info`: Provide detailed information about a user...
* `delete`: Delete a version entry in the project...

## `divbase-cli version add`

Add a new project version entry which specifies the current state of all files in the project at the current timestamp.

**Usage**:

```console
$ divbase-cli version add [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the version (e.g., semantic version).  [required]

**Options**:

* `--description TEXT`: Optional description of the version.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.

## `divbase-cli version list`

List all entries in the project versioning file.

Displays version name, creation timestamp, and description for each project version.
If you specify --include-deleted, soft-deleted versions will also be shown.
Soft-deleted versions can be restored by a DivBase admin within 30 days of deletion.

**Usage**:

```console
$ divbase-cli version list [OPTIONS]
```

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--include-deleted / --no-include-deleted`: Include soft-deleted versions in the listing.  [default: no-include-deleted]
* `--tsv`: If set, will print the output in .TSV format for easier programmatic parsing.
* `--help`: Show this message and exit.

## `divbase-cli version info`

Provide detailed information about a user specified project version, including all files present and their unique hashes.

**Usage**:

```console
$ divbase-cli version info [OPTIONS] VERSION
```

**Arguments**:

* `VERSION`: Specific version to retrieve information for  [required]

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.

## `divbase-cli version delete`

Delete a version entry in the project versioning table. This does not delete the files themselves.
Deleted version entries older than 30 days will be permanently deleted.
You can ask a DivBase admin to restore a deleted version within that time period.

**Usage**:

```console
$ divbase-cli version delete [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the version (e.g., semantic version).  [required]

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.
