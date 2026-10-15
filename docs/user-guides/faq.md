# Frequently Asked Questions

Short answers to common questions, with links to more details.

## Projects and files

### Can one project hold VCF files for different genome assemblies?

No. All VCF files in a project must be for the same version of the same reference genome. For another assembly, you need another project. See [VCF file requirements](vcf-files.md#1-vcf-file-requirements-for-divbase).

### How should I organise my VCF files?

Split large VCF files by chromosome (or by groups of scaffolds). Across files, sample sets must be either identical (same samples, same order) or not overlap at all. See [How to organise multiple VCF files](vcf-files.md#2-how-to-organize-multiple-vcf-files-in-a-divbase-project). You can divide files up into folders and subfolders, see [Uploading files](files.md#uploading-files).

### Are my files versioned?

Yes. Uploading a file with the same name adds a new version and keeps the old ones. Commands use the latest version unless you ask for an older one. See [How do I work with file versions?](files.md#how-do-i-work-with-file-versions)

### Can I save a snapshot of my whole project?

Yes. A project version records every file in the project (excluding results files) at that moment. You can download files as they were at that version later on. Files in a project version are never permanently deleted, unless the project version is deleted. See [Managing project versions](project-versioning.md).

### What happens when I delete a file?

It is soft deleted (aka hidden), and can be restored within 30 days. After 30 days it is permanently deleted, unless it is part of a project version. See [What if I want to delete a file permanently?](files.md#what-if-i-want-to-delete-a-file-permanently)

### Who can do what in a project?

Each member has one of four roles: Read, Query, Edit or Manage. Managers add and remove members on the website. See [Project roles and permissions](account-management.md#project-roles-and-permissions).

## Queries

### Why does my VCF query fail after I upload or change a VCF file?

DivBase keeps an index of what is in each VCF file (the dimensions cache). Run `divbase-cli dimensions update` after adding, replacing or restoring a VCF file, and wait for it to finish before querying. See [Dimensions update](vcf-dimensions.md#dimensions-update).

### Which bcftools commands can I use?

Only `bcftools view`, see [Writing the bcftools command](vcf-query-syntax.md#4-writing-the-bcftools-command-argument). Not all bcftools options are supported (for example, output and thread options are handled by DivBase), see [unsupported options](vcf-query-syntax.md#42-bcftools-view-commands-not-supported-by-divbase).

### Can I select samples using their metadata?

Yes. Upload a sample metadata TSV file, then filter on its columns with `divbase-cli query tsv` or `divbase-cli query vcf --tsv-filter`. `divbase-cli dimensions create-metadata-template` and `divbase-cli dimensions validate-metadata-file` help you make and check the file. See [Sidecar metadata TSV files](sidecar-metadata.md).

### Where do my query results go?

Into the project, as `result_of_job_<TASK_ID>.vcf.gz` plus a log file with format `result_of_job_<TASK_ID>.log`. Every project member can see them. Results are never used as input for other queries. See [What happens after query submission?](running-queries-overview.md#what-happens-after-query-submission)

### How long can a job run?

VCF queries and dimensions updates can run for up to 10 hours. Time spent waiting in the queue does not count. See [Maximum run times](running-queries-overview.md#maximum-run-times).

## Using DivBase

### Can I use DivBase in scripts, pipelines or HPC jobs?

Yes. Use a personal access token instead of logging in, and `divbase-cli query get-vcf-results` to wait for a query and download its result. See [Using DivBase programmatically](using-divbase-programmatically.md).

### How can I best get help from AI when using DivBase?

The entire documentation site is available as plain text for LLMs and agents. Pointing the agent at our docs site's llms.txt <https://scilifelabdatacentre.github.io/divbase/llms.txt> will give it the context it needs.

### Where can I get help?

Run `divbase-cli <command> --help`, check [Troubleshooting](troubleshooting.md), or email [datacentre@scilifelab.se](mailto:datacentre@scilifelab.se). You can also open a [GitHub issue](https://github.com/ScilifelabDataCentre/divbase/issues).
