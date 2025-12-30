# `divbase-cli query`

Run queries on the VCF files stored in the project&#x27;s data store on DivBase. Queries are run on the DivBase API

**Usage**:

```console
$ divbase-cli query [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `tsv`: Query the tsv sidecar metadata file for...
* `bcftools-pipe`: Submit a query to run on the DivBase API.

## `divbase-cli query tsv`

Query the tsv sidecar metadata file for the VCF files in the project&#x27;s data store on DivBase.
Returns the sample IDs and filenames that match the query.

TODO: it perhaps be useful to set the default download_dir in the config so that we can
look for files there? For now this code just uses file.parent as the download directory.
TODO: handle when the name of the sample column is something other than Sample_ID

**Usage**:

```console
$ divbase-cli query tsv [OPTIONS] FILTER
```

**Arguments**:

* `FILTER`: String consisting of keys:values in the tsv file to filter on.
    The syntax is &#x27;Key1:Value1,Value2;Key2:Value3,Value4&#x27;, where the key
    are the column header names in the tsv, and values are the column values.
    Multiple values for a key are separated by commas, and multiple keys are
    separated by semicolons. When multple keys are provided, an intersect query
    will be performed. E.g. &#x27;Area:West of Ireland,Northern Portugal;Sex:F&#x27;.
      [required]

**Options**:

* `--show-sample-results / --no-show-sample-results`: Print sample_ID and Filename results from the query.  [default: no-show-sample-results]
* `--metadata-tsv-name TEXT`: Name of the sample metadata TSV file in the project&#x27;s data store on DivBase.  [default: sample_metadata.tsv]
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.

## `divbase-cli query bcftools-pipe`

Submit a query to run on the DivBase API. A single, merged VCF file will be added to the project on success.

TODO Error handling for subprocess calls.
TODO: handle case empty results are returned from tsv_query()
TODO what if the user just want to run bcftools on existing files in the bucket, without a tsv file query first?
TODO what if a job fails and the user wants to re-run it? do we store temp files?
TODO be consistent about input argument and options. when are they optional, how is that indicated in docstring? etc.
TODO consider handling the bcftools command whitelist checks also on the CLI level since the error messages are nicer looking?
TODO consider moving downloading of missing files elsewhere, since this is now done before the celery task

**Usage**:

```console
$ divbase-cli query bcftools-pipe [OPTIONS]
```

**Options**:

* `--tsv-filter TEXT`: String consisting of keys:values in the tsv file to filter on.
The syntax is &#x27;Key1:Value1,Value2;Key2:Value3,Value4&#x27;, where the key
are the column header names in the tsv, and values are the column values. 
Multiple values for a key are separated by commas, and multiple keys are 
separated by semicolons. When multple keys are provided, an intersect query
will be performed. E.g. &#x27;Area:West of Ireland,Northern Portugal;Sex:F&#x27;.
* `--command TEXT`: String consisting of the bcftools command to run on the files returned by the tsv query.  [required]
* `--metadata-tsv-name TEXT`: Name of the sample metadata TSV file in the project&#x27;s data store on DivBase.  [default: sample_metadata.tsv]
* `--project TEXT`: Name of the DivBase project, if not provided uses the default in your DivBase config file
* `-c, --config PATH`: Path to your user configuration file. If you didn&#x27;t specify a custom path when you created it, you don&#x27;t need to set this.  [default: /home/roryc/.config/divbase/config.yaml]
* `--help`: Show this message and exit.
