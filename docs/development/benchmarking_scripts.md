# Benchmarking scripts

In addition to the [`pytest` testing suite](testing.md), load testing experiments have also been performed using scripts. This document describes the historical use-cases of those scripts. Please bear in mind that the code in these scripts reflect the state of the codebase at the time these experiments were performed, and that the scripts may most likely not work with the latest version of DivBase. They are presented here as a piece of history that might inform and inspire future implementation of load testing methods.

The scripts in [Section 1](#1-factorial-design-scripts-august-2025) and [Section 2](#2-cpu-and-ram-assesment-on-local-docker-compose-stack-and-on-kubernetes-deployment-january-2026) were specifically designed to investigate how VCF file sizes affect the DivBase query subsystem.

At the time of writing, scripts related to (historical) benchmarking experiments are stored in `./scripts/benchmarking`.

> **Shared helper script:** The [scripts/generate_mock_sample_metadata.py](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/generate_mock_sample_metadata.py)
    script can be used to generate a mock DivBase sidecar sample metadata TSV file from an input VCF file. It will use the Sample IDs from the VCF and generate mock metadata that can be used to make VCF queries based on sample metadata filters. Versions of this script was used to generate the metadata TSV in [Section 1](#1-factorial-design-scripts-august-2025) and [Section 2](#2-cpu-and-ram-assesment-on-local-docker-compose-stack-and-on-kubernetes-deployment-january-2026), as well as for the user guide at [Tutorial: Running a DivBase query on a public dataset](../user-guides/tutorial-query-on-public-data.md)

## 1. Factorial design scripts (August 2025)

(This experiment was also described in [Pull Request 19](https://github.com/ScilifelabDataCentre/divbase/pull/19))

At the time in the DivBase development that this experiment was performed, we had been able to show proof-of-concept of the `bcftools` orchestration logic with toy VCF. We had no results that could tell us how it would scale. In the pre-study for DivBase, there had been discussion about the expected sample pool size contained in the VCF files in  a DivBase project: the estimate was 100-1000 samples at that point in time, but that with improvements in technology it could scale to 10,000 or 100,000 sample scale. (Larger than that was consider human biobank scale and outside of the foreseeable applications of DivBase). The pre-study had also estimated that Whole-Genome Sequencing variant calling could result in many millions of variants (depending on the genome size and the population size and diversity used in the comparison). 1 million variants was considered a baseline estimate of the average DivBase use-case at the time, and with upper ranges estimated at 50-80 million variants based on the the largest public VCF datasets we were able to find at the time.

This experiment set out to investigate the following (quoted from [PR19](https://github.com/ScilifelabDataCentre/divbase/pull/19)):

- "to learn more about how the dimensions (samples x variants) of a VCF file affects performance in our current implementation add scripts to create mock VCFs and mock sample metadata so that files of varying dimensions can be tested"
- "stress-test the job system by queueing hundreds of jobs at once to identify bugs or flaws in the current implementation of the bcftools-query tasks, including interactions with buckets and the jobsystem"

To do this, a factorial design/design-of-experiments methodology was used. In short, [factorial design](https://en.wikipedia.org/wiki/Factorial_experiment) is a statistical method to see how multiple factors influence a response, and also how the factors influence and interact with each other. A Full Factorial Design will require experiments with all combinations of factors to be performed, but there are fractured factorial design methods that reduce the experiment number while still attempting to cover the same response space. A full factorial design was used here, albeit for a rather small design space.

The experiment was run on the local Docker Compose DivBase stack on a MacBook Pro M3 laptop with 32 GB RAM.

**List of files related to this experiment:**

- [`scripts/benchmarking/generate_mock_vcf.sh`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/generate_mock_vcf.sh)

    Generate mock VCF files for a full factorial design experiment. Uses the [`fake-vcf`](https://github.com/endast/fake-vcf) repository to generate mock VCF files, containerized in `docker/benchmarking.dockerfile`. For this experiment generate files for 10 different sample ranges (10-1000 samples; linspace 10 levels) and 20 different variants ranges (10-1 000 000 variants; logspace 20 levels). In total 200 mock VCF files.

- [`docker/benchmarking.dockerfile`](https://github.com/ScilifelabDataCentre/divbase/blob/main/docker/benchmarking.dockerfile)

    Docker image to run the [`fake-vcf`](https://github.com/endast/fake-vcf) tool in a container.

- [`scripts/benchmarking/factorial_design_submit_jobs.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/factorial_design_submit_jobs.py)

    For each generated mock VCF in the design space, submit the same query and measure the query wall time (time waiting in queue excluded). Run 3 replicates per query for an in total of 200 VCF files x 3 replicates = 600 jobs

- [`scripts/benchmarking/factorial_design_analyze_results.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/factorial_design_analyze_results.py)

    Fetch the results of the Celery tasks from the DivBase API, calculate the elapsed runtime, and generate plots.

- [`scripts/benchmarking/_benchmarking_shared_utils.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/_benchmarking_shared_utils.py)

    Shared logic used by the scripts of this experiment.

**Results:**

See plots and extended discussion of the results at [PR19](https://github.com/ScilifelabDataCentre/divbase/pull/19).

For this DivBase VCF query and VCF file design space, it was clear that (quoting from PR19):

> "The factorial design experiment showed that the number of variants affected the runtime more than the number of samples when it comes to the type of sample subsetting we have currently been using as an example query. Queueing 600 jobs (most of which took <1s took complete) on `celery-long` worker went well, but the queuing itself took some time."

It is perhaps not that surprising that millions of variants dominated over up-to-a-thousand samples, but it does reveal something about the nature of how `bcftools` parses VCF files.

> Note! Be careful with extrapolating these results to different queries. Depending on the subset operations and how the `bcftools` commands of the queries are writted by the user, these results may or may not hold.

## 2. CPU and RAM assesment on local Docker Compose stack and on Kubernetes deployment (January 2026)
