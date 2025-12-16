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

* `create`: Create a user configuration file for the...
* `add-project`: Add a new project to the user...
* `remove-project`: Remove a project from the user...
* `set-default`: Set the default project to use.
* `show-default`: Show the currently set default project.
* `set-dload-dir`: Set the default download dir
* `show`: Pretty print the contents of your current...

## `divbase-cli config create`

Create a user configuration file for the divbase-cli tool.

**Usage**:

```console
$ divbase-cli config create [OPTIONS]
```

**Options**:

* `--config-file PATH`: Where to store your config file locally on your pc.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli config add-project`

Add a new project to the user configuration file.

**Usage**:

```console
$ divbase-cli config add-project [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to add to the user configuration file.  [required]

**Options**:

* `--divbase-url TEXT`: DivBase API URL associated with this project.  [default: http://localhost:8000/api]
* `-d, --default`: Set this project as the default project in the user configuration file.
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli config remove-project`

Remove a project from the user configuration file.

**Usage**:

```console
$ divbase-cli config remove-project [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to remove from the user configuration file.  [required]

**Options**:

* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli config set-default`

Set the default project to use.

**Usage**:

```console
$ divbase-cli config set-default [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the project to add to the user configuration file.  [required]

**Options**:

* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli config show-default`

Show the currently set default project.

**Usage**:

```console
$ divbase-cli config show-default [OPTIONS]
```

**Options**:

* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
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

* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli config show`

Pretty print the contents of your current config file.

**Usage**:

```console
$ divbase-cli config show [OPTIONS]
```

**Options**:

* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
