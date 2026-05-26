# Key design decisions related to VCF files

This document outlines the key underlying design choices of how DivBase handles VCF files. Some, but not all of these decisions might also be covered in the  architecture decision records (ADRs) folder (`adr/`) of the GitHub repository.

## 1. All data lives in VCF files

Genetic variant data is traditionally stored in the VCF file format. There is an ecosystem of workflows for generating and analysing the data in VCF format. To align with the established upstream and downstream workflows, DivBase was designed on the core assumption that the input AND the output of the system will be VCF files, and that the overhead for the VCF I/O should be as small as possible.

However, there are several other systems, commercial and open-source, that aim to do similar things as DivBase but are based on ingesting data from VCF files into the system's internal data structure and perform all analysis on that structure. Examples include [TileDB](https://github.com/TileDB-Inc/TileDB), [XetaBase](https://zettagenomics.com/xetabase/)/[OpenGCA](https://github.com/opencb/opencga), and [VCF Zarr](https://github.com/sgkit-dev/vcf-zarr-spec). During the pre-study to what became DivBase, it was found that the performance/ wall time overhead for VCF I/O in these systems was too high for the original DivBase use cases.

This led to the choice to make all genetic variant data in DivBase live in VCF files, stored in an S3 object store with Role Based Access Control. This means that there should not be any format conversion when files pass in and out from DivBase. This differentiates DivBase from the other systems listed above. If you believe that ingestion to an internal data format is preferable to the DivBase strategy of working directly on the VCF files, we suggest you look at the open source project [VCF Zarr](https://github.com/sgkit-dev/vcf-zarr-spec).

## 2. The VCF files in DivBase are parsed by `bcftools`

There is a general bioinformatics adage that one should not attempt to write custom VCF parsers, since the VCF format has many quirks and potential pitfalls. The origin of this is likely the following documentation of one of the main tools for generating and analysing genetic variants, [GATK](https://gatk.broadinstitute.org/):

> No, really, don't write your own parser if you can avoid it. This is not a comment on how smart or how competent we think you are -- it's a comment on how annoyingly obtuse and convoluted the VCF format is.
>
> <cite>— [GATK FAQ, What is a VCF and how should I interpret it? Section 6](https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-faqs/What_is_a_VCF_and_how_should_I_interpret_it%3F.md)</cite>

The reason for this is that the VCF file format is a very flexible format. It's a tabular format that allows users to specify their own fields, making it effectively an n-dimensional format: the `FORMAT` column defines a per-row schema for the sample data, meaning each sample cell is an n-element tuple (for example, `GT:DP:GQ` gives three sub-values per sample per variant). The `INFO` column similarly allows user-defined annotations at the variant level. There is a _de facto_ [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf), but it cannot be assumed that all VCF files will adhere to that.

For this reason, all VCF processing in DivBase is done with the [`bcftools` toolkit](https://samtools.github.io/bcftools/bcftools.html). This is a compiled C binary that is considered a gold standard tool for interacting with VCF files (and its binary compressed version, BCF). Benchmarks with the DivBase deployment show that the way the `bcftools` runs are orchestrated makes it an almost completely CPU bound process, and thus the (wall) time it takes to process a VCF file with `bcftools` should be considered the performance baseline for a VCF. (Other formats, like VCF Zarr might be performance optimized in other ways, but for raw VCF file processing, assume that `bcftools` is one of the fastest tools).

The main use case for `bcftools` in DivBase is the so-called VCF queries. It is a checkout of data from all VCF files in a given DivBase project that fulfills a `bcftools view` subsetting operation. Unlike regular `bcftools` that only acts on a single file at a time, DivBase dynamically builds a workflow based on the existing VCF files in a bucket and the user's query, runs the required `bcftools` operations on all VCF files sequentially, and then creates a final result file via `bcftools merge` and/or `bcftools concat` followed by `bcftools sort`. (More about this in Section 4 below.)

## 3. Concurrency: long-running VCF jobs are executed asynchronously via a Celery job queue

VCF files can be large, not uncommonly containing 100 million variants over hundreds or thousands of samples. The processing time scales with the size of the file and the exact `bcftools view` subsetting operation. Therefore, running a DivBase VCF query can take considerable time.

To handle this, DivBase uses Celery as an asynchronous job management system, with RabbitMQ as the message broker and PostgreSQL as the results backend. When a user submits a VCF query from their command line, the DivBase API enqueues a Celery task and immediately assigns and returns a task ID to the user. The actual computation runs on a separate Celery worker pod. `worker-quick` is for tasks that are assumed to finish in less than a minute or so (e.g. sample metadata queries), and `worker-long` is for tasks that act on VCF files (VCF dimensions update, VCF queries). There is a custom task history implementation so that users can see current and historical jobs status using `divbase-cli task-history`.

To keep worker resource management simple and predictable in the k8s deployment, each Celery worker pod is configured with `concurrency=1`, which means that it will only pick up and handle one job at a time. This means that job failures are isolated: a failing pod or a failing job only affects itself. Note that `task_acks_late=False` (set in `tasks.py`) means tasks are acknowledged when picked up, not when completed; if a worker pod crashes mid-execution the task will not be automatically requeued and must be manually resubmitted. The Celery `concurrency=1` also opens up for using k8s horizontal pod autoscaling to spin up additional pod replicas if the queue is long (not implemented yet at the time of writing).

Since the S3 object store cannot be mounted to the Celery workers, each VCF related job will need to download the files it needs to the worker volume. Internal benchmarks show that the transfer time overhead is fairly small thanks to multi-part downloads from S3 and decent bandwidth on the deployment hardware. Your mileage may vary if you deploy DivBase outside of the original SciLifeLab hardware.

Side-note: DivBase also uses the Celery equivalent of cronjobs (Celery Beat) for scheduled system maintenance tasks, such as cleaning up stale task history entries. These are separate from user-submitted tasks and are not visible in the user-facing task history.

## 4. Multi-file VCF data queries are made possible by the dimensions cache

To orchestrate the `bcftools` workflow for a given VCF query, the DivBase servers need to know which sample and scaffold names that are contained in which VCF file. Since the files are stored on S3, they will need to be downloaded to the worker before they can be processed. In the most basic implementation, the DivBase server would have to download each VCF file to the worker each time a user submits a query just to find out its samples and scaffolds. This will not scale well: projects with large and/or numerous VCF files will waste resources on downloads on files that the queries will potentially never need.

To avoid this, a strategy named the "VCF dimensions cache" was implemented. This is a set of normalized PostgreSQL tables that store key technical metadata (the "dimensions") for each VCF file in a DivBase project's data store: filename, S3 version ID, sample names, sample count, scaffold names, scaffold count, variant count, file size, and a timestamp of when the entry was last updated.

The cache is populated (and updated) by submitting a "dimensions update" Celery task with `divbase-cli dimensions update` which downloads each previously unencountered VCF file, extracts the required metadata with `bcftools`, and writes the results to the Postgres table. The dimensions update task will only act on file versions not already cached in the Postgres table, meaning that if a single new file (version) has been uploaded since the last time the task was run, only that file will be processed.

Parsing the VCF files like this will take some time since `bcftools` needs to parse each file, but it is an up-front time investment that pays off on every subsequent command that needs to know the structure of the project's VCF files. At the time of writing, the current implementation requires that a user manually runs `divbase-cli dimensions update` every time that a new VCF file or a new version of a VCF file is uploaded to the DivBase project.

An up-to-date dimensions cache is a prerequisite for all query commands (sample metadata queries, and VCF queries). The server will return an error if the cached version IDs do not match the current files in the bucket and ask the user to run `divbase-cli dimensions update` to resync the dimensions cache of the DivBase project. The dimensions cache is used in several ways to ensure the data integrity and correctness of the VCF queries:

1. It guards against drift between the cache and the actual files in the bucket by exiting the query early and asking the user to run `divbase-cli dimensions update`.

2. It reduces processing overhead by only downloading the exact VCF files the query needs from S3 to the worker.

3. It determines, before any files are downloaded from S3, whether the VCF files used for the query are compatible with the `bcftools` requirements for `bcftools merge` and/or `bcftools concat` used to combine the VCF files into the single result file. This is a very important, but quite complex, decision chain that is described in more detail in [Bcftools Task Constraints](bcftools_task_constraints.md).
