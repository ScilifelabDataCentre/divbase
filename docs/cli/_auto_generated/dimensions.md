# `divbase-cli dimensions`

Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project

**Usage**:

```console
$ divbase-cli dimensions [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `update`: Calculate and add the dimensions of a VCF...
* `show`: Show the dimensions index file for a project.
* `create-metadata-template`: Use the samples index in a projects...
* `validate-metadata-file`: Validate a sidecar metadata TSV file...

## `divbase-cli dimensions update`

Calculate and add the dimensions of a VCF file to the dimensions index file in the project.

**Usage**:

```console
$ divbase-cli dimensions update [OPTIONS]
```

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.

## `divbase-cli dimensions show`

Show the dimensions index file for a project.
When running --unique-scaffolds, the sorting separates between numeric and non-numeric scaffold names.

**Usage**:

```console
$ divbase-cli dimensions show [OPTIONS]
```

**Options**:

* `--filename TEXT`: If set, will show only the entry for this VCF filename.
* `--unique-scaffolds`: If set, will show all unique scaffold names found across all the VCF files in the project.
* `--unique-samples`: If set, will show all unique sample names found across all the VCF files in the project.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.

## `divbase-cli dimensions create-metadata-template`

Use the samples index in a projects dimensions cache to create a TSV metadata template file
that has the sample names as pre-filled as the first column.

**Usage**:

```console
$ divbase-cli dimensions create-metadata-template [OPTIONS]
```

**Options**:

* `-o, --output TEXT`: Name of the output TSV file to create. Defaults to sample_metadata_&lt;project_name&gt;.tsv. If a file with the same name already exists in the current directory, you will be prompted to confirm if you want to overwrite it.
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.

## `divbase-cli dimensions validate-metadata-file`

Validate a sidecar metadata TSV file against DivBase formatting requirements and project dimensions.

Validation is run client-side to keep sensitive metadata local during validation.

Validation checks:
- File is properly tab-delimited
- First column is named &#x27;#Sample_ID&#x27;
- No commas in cells
- Sample_ID has only one value per row (no semicolons)
- No duplicate sample IDs
- Invalid characters
- Basic type consistency in user-defined columns. But not Pandas type inference,
  as we want to avoid having the user install Pandas just for validation. So just check that numeric columns have only numeric values (excluding header).
- All samples in the TSV exist in the project&#x27;s dimensions index

Returns errors for critical issues and warnings for non-critical issues.

**Usage**:

```console
$ divbase-cli dimensions validate-metadata-file [OPTIONS] INPUT_FILENAME
```

**Arguments**:

* `INPUT_FILENAME`: Name of the input TSV file to validate.  [required]

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `--help`: Show this message and exit.
