# `divbase-cli files`

Download/upload/list files to/from the project&#x27;s store on DivBase.

**Usage**:

```console
$ divbase-cli files [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `ls`: # TODO - paginate results list all files...
* `download`: Download files from the project&#x27;s store on...
* `upload`: Upload files to your project&#x27;s store on...
* `rm`: Remove files from the project&#x27;s store on...

## `divbase-cli files ls`

# TODO - paginate results
list all files in the project&#x27;s DivBase store.

To see files at a user specified project version (controlled by the &#x27;divbase-cli version&#x27; subcommand),
you can instead use the &#x27;divbase-cli version info [VERSION_NAME]&#x27; command.

**Usage**:

```console
$ divbase-cli files ls [OPTIONS]
```

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files download`

Download files from the project&#x27;s store on DivBase. This can be done by either:
    1. providing a list of files paths directly in the command line
    2. providing a directory to download the files to.

**Usage**:

```console
$ divbase-cli files download [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space separated list of files/objects to download from the project&#x27;s store on DivBase.

**Options**:

* `--file-list PATH`: Text file with list of files to upload.
* `--download-dir TEXT`: Directory to download the files to. 
If not provided, defaults to what you specified in your user config. 
If also not specified in your user config, downloads to the current directory.
You can also specify &quot;.&quot; to download to the current directory.
* `--disable-verify-checksums`: Turn off checksum verification which is on by default. Checksum verification means all downloaded files are verified against their MD5 checksums.It is recommended to leave checksum verification enabled unless you have a specific reason to disable it.
* `--project-version TEXT`: User defined version of the project&#x27;s at which to download the files. If not provided, downloads the latest version of all selected files.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files upload`

Upload files to your project&#x27;s store on DivBase by either:
    1. providing a list of files paths directly in the command line
    2. providing a directory to upload
    3. providing a text file with or a file list.

**Usage**:

```console
$ divbase-cli files upload [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space seperated list of files to upload.

**Options**:

* `--upload-dir PATH`: Directory to upload all files from.
* `--file-list PATH`: Text file with list of files to upload.
* `--disable-safe-mode`: Turn off safe mode which is on by default. Safe mode adds 2 extra bits of security by first calculating the MD5 checksum of each file that you&#x27;re about to upload:(1) Checks if any of the files you&#x27;re about to upload already exist (by comparing name and checksum) and if so stops the upload process.(2) Sends the file&#x27;s checksum when the file is uploaded so the server can verify the upload was successful (by calculating and comparing the checksums).It is recommended to leave safe mode enabled unless you have a specific reason to disable it.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files rm`

Remove files from the project&#x27;s store on DivBase by either:
    1. providing a list of files paths directly in the command line
    2. providing a text file with or a file list.

&#x27;dry_run&#x27; mode will not actually delete the files, just print what would be deleted.

**Usage**:

```console
$ divbase-cli files rm [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space seperated list of files/objects in the project&#x27;s store on DivBase to delete.

**Options**:

* `--file-list PATH`: Text file with list of files to upload.
* `--dry-run`: If set, will not actually delete the files, just print what would be deleted.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
