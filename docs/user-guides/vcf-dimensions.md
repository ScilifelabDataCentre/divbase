# VCF Dimensions caching

The DivBase query system is built around the DivBase server having cached key technical metadata from each VCF file in a project. In DivBase, this information as the "VCF dimensions", and for instance includes the names of samples and scaffolds in each VCF file and the version ID of the file in the project data store. This allows the DivBase server to make quick checks against against the project VCF Dimensions when users submit a query or validate a sidecar metadata TSV file instead of having to read each VCF file every time a query is submitted. The VCF Dimensions is a snapshot of the VCF files in the object store at the time the command to update the VCF Dimensions cache for the project was last run.

**An updated VCF Dimensions cache for DivBase project is a prerequisite step before submitting any queries**. Updating it is done with the command `divbase-cli dimensions update`. The command needs to be run every time the VCF files in a DivBase project has changed, that is:

- When any new VCF file is uploaded to the project.

- When an existing VCF file is replaced with a new version.

!!! Warning
    DivBase system does not automatically run `divbase-cli dimensions update` when files are uploaded. This needs to be manually done by one project member that has at least an EDIT role in the project.

There are four CLI commands sorted under `divbase-cli dimensions`:

- [update](#dimensions-update)

- [show]((#dimensions-show))

- [create-metadata-template](sidecar-metadata.md#creating-a-sidecar-sample-metadata-tsv-for-a-divbase-project)

- [validate-metadata-file](sidecar-metadata.md#validating-a-sidecar-metadata-tsv-with-divbase-cli)

The former two relate generating/updating and viewing the VCF dimensions of a project and will be described on this page. The latter two relate to the user-defined sidecar metadata TSV and are described in the guide on [Sidecar sample metadata TSV queries](sidecar-metadata.md). The commands are sorted under `divbase-cli dimensions` since they rely on calling the VCF dimensions cache of the project. A list of all DivBase CLI commands that require an up-to-date VCF Dimensions cache is found in the section [CLI Commands that rely on that the project's VCF dimensions cache is up to date](#cli-commands-that-rely-on-that-the-projects-vcf-dimensions-cache-is-up-to-date)

## Dimensions update

The `divbase-cli dimensions update` pre-reads the VCF files in the project's data store so that DivBase does not need to fetch and read each VCF every time the user submits commands that require knowledge about the VCF dimensions, such as queries.

The `divbase-cli dimensions update` will look for `.vcf.gz` files in the project's data store. If there is no VCF dimensions cache for the project, it will create it. If not, it will compare the existing record in the cache with the current status of the object store. If any new VCF files have been added or if any VCF file version have been updated, it will update the VCF dimensions cache with that information.

The update of the VCF Dimensions is scheduled as a job in the DivBase job queue system since it can a potentially take little time to update a VCF dimensions for projects that contain many or large VCF files that has not been previously cached. This is an up-front time investment: the time it takes to run the update command is saved on every subsequent command that needs to check the VCF dimensions.

VCF results files that have been produced by DivBase are not indexed in the VCF Dimensions cache. The main reason for this is that they result files contain a subset of the VCF data in the project, and will thus contain duplicate data. The results files are not used for any queries: only the source VCF files uploaded by the users are. DivBase recognizes its results files on two levels: the files names have a `result_of_job_` prefix, and their VCF headers contain a row with `##DivBase_created`.

The VCF Dimensions cache stores this information for each VCF file in the project's data store:

- filename (`*.vcf.gz`)
- version ID in the data store (each file name is versioned and has a version history. The VCF Dimensions cache uses the latest version of all files in the data store.)
- time stamp of when the VCF dimensions for this file was updated
- sample names contained in this VCF file
- sample count
- scaffold names contained in this VCF file
- scaffold count
- variant count
variant_count, sample_count, file_size_bytes
updated_at (timestamp of last indexing — this answers the TODO in the current draft: it is when dimensions was last run for that file, not when the file was last modified in S3)
Related rows in vcf_metadata_samples and vcf_metadata_scaffolds (normalized in the latest migration)

## Dimensions show

The skipped DivBase result files are tracked by DivBase and will be displayed

--filename
--unique-scaffolds
--unique-samples
--sample-names-limit
--sample-names-output
--sample-names-stdout

## CLI Commands that rely on that the project's VCF dimensions cache is up to date

Several DivBase CLI commands required that the VCF dimensions cache of the project is up-to-date with the current versions of the VCF files in the project's data store:

`divbase-cli dimensions show`
`divbase-cli dimensions create-metadata-template`
`divbase-cli dimensions validate-metadata-file`
`divbase-cli query tsv` Raises VCFDimensionsEntryMissingError if index is empty; raises stale-data error if version IDs don't match bucket
`divbase-cli query bcftools-pipe` if run as a combined metadata VCF query, the above. Also uses dimensions to route the query to the correct subset of VCF files per sample, and to filter out irrelevant files based on sample-filename mapping and on scaffolds for -r region queries
