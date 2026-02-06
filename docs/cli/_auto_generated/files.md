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

* `ls`: list all currently available files in the...
* `info`: Get detailed information about a specific...
* `download`: Download files from the project&#x27;s store on...
* `stream`: Stream a file&#x27;s content to standard output.
* `upload`: Upload files to your project&#x27;s store on...
* `rm`: Soft delete files from the project&#x27;s store...
* `restore`: Restore soft deleted files from the...

## `divbase-cli files ls`

list all currently available files in the project&#x27;s DivBase store.

You can optionally filter the listed files by providing a prefix.
By default, DivBase query results files are hidden from the listing. Use the --include-results-files option to include them.
To see information about the versions of each file, use the &#x27;divbase-cli files info [FILE_NAME]&#x27; command instead

**Usage**:

```console
$ divbase-cli files ls [OPTIONS]
```

**Options**:

* `--tsv`: If set, will print the output in .TSV format for easier programmatic parsing.
* `-p, --prefix TEXT`: Optional prefix to filter the listed files by name (only list files starting with this prefix).
* `-r, --include-results-files`: If set, will also show DivBase query results files which are hidden by default.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files info`

Get detailed information about a specific file in the project&#x27;s DivBase store.

This includes all versions of the file and whether the file is currently marked as soft deleted.

**Usage**:

```console
$ divbase-cli files info [OPTIONS] FILE_NAME
```

**Arguments**:

* `FILE_NAME`: Name of the file to get information about.  [required]

**Options**:

* `--tsv`: If set, will print the output in .TSV format for easier programmatic parsing.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files download`

Download files from the project&#x27;s store on DivBase.

This can be done by either:
    1. providing a list of files paths directly in the command line
    2. providing a text file with a list of files to download (new file on each line).

To download the latest version of a file, just provide its name. &quot;file1&quot; &quot;file2&quot; etc.
To download a specific/older version of a file, use the format: &quot;file_name:version_id&quot;
You can get a file&#x27;s version id using the &#x27;divbase-cli file info [FILE_NAME]&#x27; command.
You can mix and match latest and specific versions in the same command.
E.g. to download the latest version of file1 and version &quot;3xcdsdsdiw829x&quot;
of file2: &#x27;divbase-cli files download file1 file2:3xcdsdsdiw829x&#x27;

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

## `divbase-cli files stream`

Stream a file&#x27;s content to standard output.

This allows your to pipe the output to other tools like &#x27;less&#x27;, &#x27;head&#x27;, &#x27;zcat&#x27; and &#x27;bcftools&#x27;.

Examples:
- View a file: divbase-cli files stream my_file.tsv | less
- View a gzipped file: divbase-cli files stream my_file.vcf.gz | zcat | less
- Run a bcftools command: divbase-cli files stream my_file.vcf.gz | bcftools view -h -  # The &quot;-&quot; tells bcftools to read from standard input

**Usage**:

```console
$ divbase-cli files stream [OPTIONS] FILE_NAME
```

**Arguments**:

* `FILE_NAME`: Name of the file you want to stream.  [required]

**Options**:

* `--version-id TEXT`: Specify this if you want to look at an older/specific version of the file. If not provided, the latest version of the file is used. To get a file&#x27;s version ids, use the &#x27;divbase-cli file info [FILE_NAME]&#x27; command.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files upload`

Upload files to your project&#x27;s store on DivBase:

To provide files to upload you can either:
    1. provide a list of files paths directly in the command line
    2. provide a directory to upload
    3. provide a text file with or a file list.

**Usage**:

```console
$ divbase-cli files upload [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space separated list of files to upload.

**Options**:

* `--upload-dir PATH`: Directory to upload all files from.
* `--file-list PATH`: Text file with list of files to upload.
* `--disable-safe-mode`: Turn off safe mode which is on by default. Safe mode adds 2 extra bits of security by first calculating the MD5 checksum of each file that you&#x27;re about to upload:(1) Checks if any of the files you&#x27;re about to upload already exist (by comparing name and checksum) and if so stops the upload process.(2) Sends the file&#x27;s checksum when the file is uploaded so the server can verify the upload was successful (by calculating and comparing the checksums).It is recommended to leave safe mode enabled unless you have a specific reason to disable it.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files rm`

Soft delete files from the project&#x27;s store on DivBase

To provide files to delete you can either:
    1. provide a list of file names directly in the command line
    2. provide a text file with a list of files to delete.

Note that deleting a non existent file will be treated as a successful deletion.

**Usage**:

```console
$ divbase-cli files rm [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space seperated list of files/objects in the project&#x27;s store on DivBase to delete.

**Options**:

* `--file-list PATH`: Text file with list of files to delete.
* `--dry-run`: If set, will not actually delete the files, just print what would be deleted.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli files restore`

Restore soft deleted files from the project&#x27;s store on DivBase

To provide files to restore you can either:
    1. provide a list of files directly in the command line.
    2. provide a text file with a list of files to restore (new file on each line).

NOTE: Attempts to restore a file that is not soft deleted will be considered successful and the file will remain live. This means you can repeatedly run this command on the same file and get the same response.

**Usage**:

```console
$ divbase-cli files restore [OPTIONS] [FILES]...
```

**Arguments**:

* `[FILES]...`: Space seperated list of files/objects in the project&#x27;s store on DivBase to restore.

**Options**:

* `--file-list PATH`: Text file with list of files to restore.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
