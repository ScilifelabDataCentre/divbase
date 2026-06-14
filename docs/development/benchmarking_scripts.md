# Benchmarking scripts

In addition to the [`pytest` testing suite](testing.md), load testing experiments have also been performed using scripts. This document describes the historical use-cases of those scripts. Please bear in mind that the code in these scripts reflect the state of the codebase at the time these experiments were performed, and that the scripts will most likely not work with the latest version of DivBase. They are presented here as a piece of history that might inform and inspire future implementation of load testing methods.

The scripts in [Section 1](#1-factorial-design-scripts-august-2025) and [Section 2](#2-cpu-and-ram-assessment-on-local-docker-compose-stack-and-on-kubernetes-deployment-january-2026) were specifically designed to investigate how VCF file sizes affect the DivBase query subsystem.

At the time of writing, scripts related to (historical) benchmarking experiments are stored in `./scripts/benchmarking`.

!!! Note
    The [scripts/generate_mock_sample_metadata.py](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/generate_mock_sample_metadata.py) script can be used to generate a mock DivBase sidecar sample metadata TSV file from an input VCF file. It will use the Sample IDs from the VCF and generate mock metadata that can be used to make VCF queries based on sample metadata filters. Versions of this script were used to generate the metadata TSV in [Section 1](#1-factorial-design-scripts-august-2025) and [Section 2](#2-cpu-and-ram-assessment-on-local-docker-compose-stack-and-on-kubernetes-deployment-january-2026), as well as for the user guide at [Tutorial: Running a DivBase query on a public dataset](../user-guides/tutorial-query-on-public-data.md)

## 1. Factorial design scripts (August 2025)

(This experiment was also described in [Pull Request 19](https://github.com/ScilifelabDataCentre/divbase/pull/19))

At the time in the DivBase development that this experiment was performed, we had been able to show proof-of-concept of the `bcftools` orchestration logic with toy VCF. We had no results that could tell us how it would scale. In the pre-study for DivBase, there had been discussion about the expected sample pool size contained in the VCF files in a DivBase project: the estimate was 100-1000 samples at that point in time, but that with improvements in technology it could scale to 10,000 or 100,000 sample scale. (Larger than that was considered human biobank scale and outside of the foreseeable applications of DivBase). The pre-study had also estimated that Whole-Genome Sequencing variant calling could result in many millions of variants (depending on the genome size and the population size and diversity used in the comparison). 1 million variants was considered a baseline estimate of the average DivBase use-case at the time, and with upper ranges estimated at 50-80 million variants based on the largest public VCF datasets we were able to find at the time.

This experiment set out to investigate the following (quoted from [PR19](https://github.com/ScilifelabDataCentre/divbase/pull/19)):

- "to learn more about how the dimensions (samples x variants) of a VCF file affects performance in our current implementation add scripts to create mock VCFs and mock sample metadata so that files of varying dimensions can be tested"
- "stress-test the job system by queueing hundreds of jobs at once to identify bugs or flaws in the current implementation of the bcftools-query tasks, including interactions with buckets and the jobsystem"

To do this, a factorial design/design-of-experiments methodology was used. In short, [factorial design](https://en.wikipedia.org/wiki/Factorial_experiment) is a statistical method to see how multiple factors influence a response, and also how the factors influence and interact with each other. A Full Factorial Design will require experiments with all combinations of factors to be performed, but there are fractional factorial design methods that reduce the experiment number while still attempting to cover the same response space. A full factorial design was used here, albeit for a rather small design space.

The experiment was run on the local Docker Compose DivBase stack on a MacBook Pro M3 laptop with 32 GB RAM.

**List of files related to this experiment:**

- [`scripts/benchmarking/generate_mock_vcf.sh`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/generate_mock_vcf.sh)

    Generate mock VCF files for a full factorial design experiment. Uses the [`fake-vcf`](https://github.com/endast/fake-vcf) repository to generate mock VCF files, containerized in `docker/benchmarking.dockerfile`. For this experiment generate files for 10 different sample ranges (10-1000 samples; linspace 10 levels) and 20 different variants ranges (10-1 000 000 variants; logspace 20 levels). In total 200 mock VCF files.

- [`docker/benchmarking.dockerfile`](https://github.com/ScilifelabDataCentre/divbase/blob/main/docker/benchmarking.dockerfile)

    Docker image to run the [`fake-vcf`](https://github.com/endast/fake-vcf) tool in a container.

- [`scripts/benchmarking/factorial_design_submit_jobs.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/factorial_design_submit_jobs.py)

    For each generated mock VCF in the design space, submit the same query and measure the query wall time (time waiting in queue excluded). Run 3 replicates per query for a total of 200 VCF files x 3 replicates = 600 jobs. In this experiment, the script was used to run a full factorial design, but a Latin Hypercube method to reduce sample space was also toyed with, but not used for generating results described in [Pull Request 19](https://github.com/ScilifelabDataCentre/divbase/pull/19).

- [`scripts/benchmarking/factorial_design_analyze_results.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/factorial_design_analyze_results.py)

    Fetch the results of the Celery tasks from the DivBase API, calculate the elapsed runtime, and generate plots.

- [`scripts/benchmarking/_benchmarking_shared_utils.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/_benchmarking_shared_utils.py)

    Shared logic used by the scripts of this experiment.

**Results:**

See plots and extended discussion of the results at [PR19](https://github.com/ScilifelabDataCentre/divbase/pull/19).

For this DivBase VCF query and VCF file design space, it was clear that (quoting from PR19):

> "The factorial design experiment showed that the number of variants affected the runtime more than the number of samples when it comes to the type of sample subsetting we have currently been using as an example query. Queueing 600 jobs (most of which took <1s took complete) on `celery-long` worker went well, but the queuing itself took some time."

It is perhaps not that surprising that millions of variants dominated over up-to-a-thousand samples, but it does reveal something about the nature of how `bcftools` parses VCF files.

> Note! Be careful with extrapolating these results to different queries. Depending on the subset operations and how the `bcftools` commands of the queries are written by the user, these results may or may not hold.

## 2. CPU and RAM assessment on local Docker Compose stack and on Kubernetes deployment (January 2026)

This experiment investigated the resource demand of running a DivBase VCF query on a medium-sized public VCF file. The aim was to investigate CPU and RAM resource usage during VCF queries in detail and use that to identify appropriate Kubernetes resource `requests` and `limits` settings.

A hard requirement was to identify Kubernetes memory settings that would not lead to the pod running the query task being terminated with an OOMKill (out-of-memory-killed) signal. Such termination would mean that the jobs could never be finished within the given memory constraints even if they were resubmitted.

A secondary requirement was to try to optimize the wall time of the running job (i.e. the time it takes to start and complete the query job; time spent in queue excluded). This was primarily done by investigating how the `bcftools` compression format during file I/O affects performance. Secondarily, the effect on query wall time of splitting a large monolithic VCF into smaller VCF files by chromosome was also investigated. This idea is based on the VCF Dimensions Cache mapping sample and scaffolds (here: chromosomes) to filename and how that is used to only download and process the required files during query runtime.

The experiment was run on two environments: locally using the Docker Compose DivBase development stack on a MacBook Pro M3 Pro laptop with 36 GB RAM, and remotely on the `scilifelab-2-dev` Kubernetes cluster.

### 2.1. Inputs (data and query)

The publicly deposited _Mus musculus_ SNP dataset `mgp.v3.snps.rsIDdbSNPv137.vcf.gz` was
downloaded from the European Nucleotide Archive (ENA) at accession [ERZ022025](https://www.ebi.ac.uk/ena/browser/view/ERZ022025). It is 5.5 GB in size and contains 18 samples and 66,007,044 variants. This is henceforth referred to as the _monolithic_ version of this VCF file.

A chromosome-split version of the same file was created for comparison using [`scripts/benchmarking/split_mouse_vcf_per_scaffold.sh`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/split_mouse_vcf_per_scaffold.sh),
producing 20 per-chromosome files (`mgp.v3.snps.rsIDdbSNPv137.<chromosome>.vcf.gz`, chromosomes
1–20). These are referred to as the _split_ versions of the `mgp.v3.snps.rsIDdbSNPv137.vcf.gz` file.

A mock sample metadata TSV (`mock_metadata_mgpv3snps.tsv`) was generated for `mgp.v3.snps.rsIDdbSNPv137.vcf.gz` with [`scripts/generate_mock_sample_metadata.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/generate_mock_sample_metadata.py), assigning each of the 18 samples a mock `Population`, `Area`, and `Sex` field. The same file was used for the monolithic and split files, since they contained the same samples. (Note that this script has been refactored since the experiment was conducted, and that newer versions of the script will produce slightly different fields than the ones described here).

The files were uploaded to two separate buckets: one with the monolithic file, and one with all the split files. The same sample metadata file was uploaded to both buckets. The reason for the separation into two buckets is that the VCF dimensions caching will not accept data duplication in the same bucket, and in this case, the data is exactly the same in the two buckets but distributed across VCFs differently.

The same query was used throughout the experiment for all the factors varied. It subsets the VCF data down to 10 samples with the `Area:North,East` sample metadata filter, and further subsets the variants to only include those with coordinates within: chromosome 1, range 15–25 Mbp. This produces a results file containing approximately 297,000 variants. The CLI input to enqueue the query was:

```bash
# Note! The CLI syntax has been slightly refactored since this experiment
divbase-cli query bcftools-pipe \
  --tsv-filter 'Area:North,East' \
  --command 'view -s SAMPLES; view -r 1:15000000-25000000' \
  --metadata-tsv-name mock_metadata_mgpv3snps.tsv \
  --project <project_name>
```

### 2.2. Monitoring setup

For full metric definitions (wall time; CPU time; peak and average RSS memory; broken down per Python process, `bcftools` subprocesses, and S3 download step), see [Monitoring: Custom metrics collection for Celery worker](worker_metrics.md).

Prior to this experiment there had been some initial experiments with measuring CPU and memory using Prometheus globally for each service in the DivBase stack. The insight from that was that per-task granularity would be more useful for understanding query performance, so the VCF query task in [`packages/divbase-api/src/divbase_api/worker/tasks.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/packages/divbase-api/src/divbase_api/worker/tasks.py) was expanded to contain per-task metrics collection logic using the [`prometheus-client`](https://github.com/prometheus/client_python) and `psutil` Python libraries. Implementation details for the per-task metrics monitoring are described in [Monitoring: Custom metrics collection for Celery worker](worker_metrics.md). Below are the details on how this was run for the experiment.

The per-task metrics are emitted by the Celery worker at a custom `prometheus-client` endpoint on port 8101. Two environment variables control this:

```yaml
ENABLE_WORKER_METRICS=1          # enable the metrics server and basic worker metrics
ENABLE_WORKER_METRICS_PER_TASK=1 # also collect per-task CPU, memory, and wall time
```

**For the local Docker Compose environment**, the monitoring stack (Prometheus, Grafana, cAdvisor) is started separately from the main stack using the dedicated monitoring Compose file. The Docker network `divbase-observability` allows the monitoring stack to listen to the main stack.

First, enable the per-task metrics in the main stack. Both `ENABLE_WORKER_METRICS` and `ENABLE_WORKER_METRICS_PER_TASK` are set to `0` in `docker/divbase_compose.yaml` by default and need to be changed to `1` before starting the stack.

Then spin up the main stack, and when that is running, spin up the monitoring stack:

```bash
# Start the main stack
docker compose -f docker/divbase_compose.yaml up -d

# Start the monitoring stack
docker compose -f docker/monitoring_compose.yaml up -d
```

This will start Prometheus directly on `localhost:9090`.

**For the Kubernetes dev cluster environment**, Prometheus is deployed inside the cluster as a `ClusterIP` service (not publicly exposed). To run the metrics collection script against it, port-forward the cluster Prometheus to localhost as described in the [private deployment repository](https://github.com/ScilifelabDataCentre/argocd-divbase).

### 2.3. Factors varied

- Kubernetes CPU and memory resource settings (`requests`/`limits`). Multiple configurations were tested, from an initial low estimate upwards.
- Kubernetes temporary file storage: `tmpfs` (in-memory) vs the default node filesystem for intermediate query files.
- `bcftools` intermediate compression format: gzip-compressed BCF (`-Oz`), BGZF-compressed BCF (`-Ob`), and uncompressed BCF (`-Ou`).
- VCF data structure: single monolithic VCF file vs chromosome-split VCFs (20 files; one per chromosome).

### 2.4. List of files related to this experiment

- [`scripts/benchmarking/run_mouse_vcf_job_in_docker_compose.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/run_mouse_vcf_job_in_docker_compose.py)

End-to-end workflow script for the **local Docker Compose runs**. Creates the MinIO bucket and DivBase project if needed, downloads the mouse VCF from ENA, generates mock sample metadata using [`scripts/generate_mock_sample_metadata.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/generate_mock_sample_metadata.py), uploads both files via `divbase-cli files upload`, runs a dimensions update, and submits the query job. The script does not poll for the task to complete: the user has to check for task completion manually. The query took approximately 10 minutes on the laptop used for the local tests. Usage:

```bash
docker compose -f docker/divbase_compose.yaml up -d
python scripts/benchmarking/run_mouse_vcf_job_in_docker_compose.py
```

!!! Note
    For the Kubernetes side of this experiment, all commands were run manually rather than to write a reproducible script. This was considered more convenient as the Kubernetes deployment was undergoing a lot of refactoring at the time.

- [`scripts/benchmarking/split_mouse_vcf_per_scaffold.sh`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/split_mouse_vcf_per_scaffold.sh)

Script that downloads the mouse VCF from ENA (if not already present), creates a CSI index (required for the splitting), and splits the file by chromosome using `bcftools view -r`. For convenience, also outputs a file list (`split_scaffold_files.txt`) for use with `divbase-cli files upload`. Current state of the script requires that `bcftools` is installed locally. Usage:

```bash
# Note that the mouse VCF URL is hardcoded.
# This is not a general-purpose script to split VCFs by chromosomes (but it could be refactored to become one)
bash scripts/benchmarking/split_mouse_vcf_per_scaffold.sh
```

- [`scripts/benchmarking/fetch_per_task_metrics_from_prometheus.py`](https://github.com/ScilifelabDataCentre/divbase/blob/main/scripts/benchmarking/fetch_per_task_metrics_from_prometheus.py)

Script that fetches the per-task metrics from the Prometheus data store for a list of completed job IDs and dates specified in a YAML file. Requires Prometheus to be reachable on localhost:9090: either via `docker/monitoring_compose.yaml` (Docker Compose) or via port-forwarding to the Kubernetes cluster.

To make a task metrics lookup in the Prometheus database, the DivBase Task ID and the date the task was run is needed. To conveniently pull several entries from Prometheus at the same time, the script was designed to take a YAML file as input. Example of a input YAML file for this script:

```yaml
jobs:
  - job_id: 23
    date: 2026-01-28
    comment: "Optional, but will be saved to the output JSON if included"
  - job_id: 28
    date: 2026-01-29
```

Assuming the tasks and dates are correct, the script will then fetch the metrics results from all jobs and write them to `scripts/benchmarking/results/fetched_metrics_<yaml_name>.json`.

Example command for an input YAML named `task_metrics_docker_local.yaml`:

```bash
python scripts/benchmarking/fetch_per_task_metrics_from_prometheus.py \
    --yaml scripts/benchmarking/task_metrics_docker_local.yaml
```

### 2.5. Results

The absolute results from this experiment are considered private for staff only. But in summary, the key findings from the experiment were:

- The `bcftools view` process seem to be mainly CPU-bound. Setting too low `requests`/`limits` CPU values in the Kubernetes deployment led to throttled jobs with longer wall time. Sufficient CPU values for a VCF file of this size was identified. Since Kubernetes will throttle - but not terminate - jobs that exceed the CPU limit, the current values were considered acceptable for the time being. Also worthy of noting is that in the current DivBase implementation, `bcftools` subprocesses run sequentially and single-threaded, so allocating more CPU cores will not reduce wall time.

- Some Kubernetes memory resource configurations were too low and caused a OOMKill signal and resulted in pod termination. Sufficient memory `requests`/`limits` values were implemented to avoid this. This mainly happened when trying to optimize the run by using the `tmpfs` in-memory file system to decrease the wall time of the tasks. The `tmpfs` did decrease the wall time, but uses a fixed value set at deployment and will thus not be dynamically scalable. The conclusion is that `tmpfs` should not be used for this purpose in DivBase deployments.

- Changing the `bcftools` intermediate file output from gzip-compressed (`-Oz`) to uncompressed BCF (`-Ou`) cut wall time by roughly half on both the local Docker Compose stack and on the Kubernetes cluster. The cost for this was higher peak memory usage and more disk storage needed for the uncompressed temp files (compared to their compressed dito), but the quantities were not excessive and thus were considered an acceptable price. The reason for the improvement seem to be that `bcftools` will always operate on uncompressed BCF, so by excluding the compression step, substantial computational resources can be saved. Note that the input and output VCFs to the query task still are gzip-compressed (`-Oz`) VCF.GZ since this is the standard way to store VCF data; only the temp file compression has been changed to uncompressed BCF.

- The DivBase overhead processes during the task (i.e. all subprocesses after excluding `bcftools` subprocesses) accounted for approximately 4% of total wall time in the best-case condition. Almost all of this overhead consisted of the time needed for file transfer from S3, which is a known bottleneck in the system design. That most wall time was spent on `bcftools` operations means that the implementation can be considered efficient: a DivBase VCF query can never be faster than the time it takes to run its required `bcftools` steps.
