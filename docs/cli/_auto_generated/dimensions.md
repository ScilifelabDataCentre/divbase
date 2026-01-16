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

## `divbase-cli dimensions update`

Calculate and add the dimensions of a VCF file to the dimensions index file in the project.

**Usage**:

```console
$ divbase-cli dimensions update [OPTIONS]
```

**Options**:

* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
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
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
