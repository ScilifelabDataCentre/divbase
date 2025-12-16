# `divbase-cli version`

Version the state of all files in the entire projects storage bucket at a given timestamp.

**Usage**:

```console
$ divbase-cli version [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `create`: Create a bucket versioning file that is...
* `add`: Add an entry to the bucket versioning file...
* `list`: List all entries in the bucket versioning...
* `delete`: Delete an entry in the bucket versioning...
* `info`: Provide detailed information about a user...

## `divbase-cli version create`

Create a bucket versioning file that is stored inside the project&#x27;s storage bucket.

**Usage**:

```console
$ divbase-cli version create [OPTIONS]
```

**Options**:

* `--name TEXT`: Name of the version (e.g., semantic version).  [default: v0.0.0]
* `--description TEXT`: Optional description of the version.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli version add`

Add an entry to the bucket versioning file specfying the current state of all files in the project&#x27;s storage bucket.

**Usage**:

```console
$ divbase-cli version add [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the version (e.g., semantic version).  [required]

**Options**:

* `--description TEXT`: Optional description of the version.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli version list`

List all entries in the bucket versioning file.

**Usage**:

```console
$ divbase-cli version list [OPTIONS]
```

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli version delete`

Delete an entry in the bucket versioning file specfying a specific state of all files in the project&#x27;s storage bucket.
Does not delete the files themselves.

**Usage**:

```console
$ divbase-cli version delete [OPTIONS] NAME
```

**Arguments**:

* `NAME`: Name of the version (e.g., semantic version).  [required]

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
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
* `-c, --config PATH`: Path to your user configuration file. By default it is stored at ~/.config/.divbase_tools.yaml.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
