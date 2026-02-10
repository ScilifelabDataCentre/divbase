# `divbase-cli config`

Manage your user configuration file for the DivBase CLI.

**Usage**:

```console
$ divbase-cli config [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `add`: Add a new project to your user...
* `rm`: Remove a project from your user...
* `set-default`: Set your default project to use in all...
* `show-default`: Print the currently set default project to...
* `set-dload-dir`: Set the default download dir
* `show`: Pretty print the contents of your current...

## `divbase-cli config add`

Add a new project to your user configuration file.

**Usage**:

```console
$ divbase-cli config add [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to add to your config file.  [required]

**Options**:

* `-u, --divbase-url TEXT`: DivBase API URL associated with this project.  [default: http://localhost:8000/api]
* `-d, --default`: Set this project as the default project in your config file.
* `--help`: Show this message and exit.

## `divbase-cli config rm`

Remove a project from your user configuration file.

**Usage**:

```console
$ divbase-cli config rm [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to remove from your user configuration file.  [required]

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli config set-default`

Set your default project to use in all divbase-cli commands.

**Usage**:

```console
$ divbase-cli config set-default [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to add to the user configuration file.  [required]

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli config show-default`

Print the currently set default project to the console.

**Usage**:

```console
$ divbase-cli config show-default [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli config set-dload-dir`

Set the default download dir

**Usage**:

```console
$ divbase-cli config set-dload-dir [OPTIONS] DOWNLOAD_DIR
```

**Arguments**:

* `DOWNLOAD_DIR`: Set the default directory to download files to. 
        By default files are downloaded to the current working directory.
        You can specify an absolute path. 
        You can use &#x27;.&#x27; to refer to the directory you run the command from.  [required]

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli config show`

Pretty print the contents of your current config file.

**Usage**:

```console
$ divbase-cli config show [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.
