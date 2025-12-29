# `divbase-cli`

This tool lets you interact with your DivBase project(s) in order to:

    - Query the metadata for the VCF files stored in the project.

    - Upload/download files to/from the project. 

    - Version the state of all files in the entire project at a given timestamp.

**Usage**:

```console
$ divbase-cli [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `version`: Add, view and remove versions representing...
* `files`: Download/upload/list files to/from the...
* `config`: Manage your user configuration file for...
* `query`: Run queries on the VCF files stored in the...
* `dimensions`: Create and inspect dimensions (number of...
* `auth`: Login/logout of DivBase server.
* `task-history`: Get the task history of query jobs...