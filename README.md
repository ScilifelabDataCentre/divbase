# divBase

**DivBase's goal is to be a service which is comprised of:**

- An S3 object store to store researchers VCF files and associated metadata.
- A CLI installed by users to interact with their DivBase projects and the webAPI (next bullet).
- A webAPI, installed on SciLifeLab cluster to run query jobs and create "results files" from these queries.
  - The results files are added to the object store by and can be downloaded by the researchers.
- As you may have guessed, the CLI can therefore requests to the webAPI including making queries.submitting query jobs and can download/upload files from the object store.

## Table of Contents

1. [Folders overview](#folders-overview)
2. [Minio Deployment on Dev cluster](#minio-deployment-on-dev-cluster)
3. [divbase CLI](#divbase-cli)
4. [Developer Setup](#developer-setup)
5. [Queries](#queries)

## Folders overview

- *docker:* bcftools docker image.
- *docs:* Documentation
- *kustomize:* Deploy the MinIO instance for the DivBase prototype using Kustomize + Helm.
- *src:* Source code for divbase's CLI and webAPI.
- *tests*: Tests (with Pytest).

## Minio Deployment on Dev cluster

Deployment for testing/development purposes. The deployment and how it was set up is covered in more detail at [docs/bitnami-minio-setup.md](docs/bitnami-minio-setup.md).

## divBase CLI

The CLI tool `divbase-cli` is designed to help you interact with your DivBase project. The CLI is built using [Typer](https://typer.tiangolo.com).

- See the basic usage/documentation at [docs/using_divbase_cli.md](docs/using_divbase_cli.md)
- See the [Developer Setup](#developer-setup) section to install the CLI tool in develop mode.

## Developer Setup

### 1. Setup your python virtual environment

You can use either [uv](https://github.com/astral-sh/uv) or something like `venv` or `pip`.

#### Using `uv` (Recommended)

- Install `uv`, follow the docs: <https://docs.astral.sh/uv/>

- In the root of the repository, run:

```bash
uv sync
```

This will create a virtual environment at the root of the repo and install all dependencies and the package (with all code in the src folder installed in "editable" mode).

#### Using pip and venv

```bash
python3 -m venv .venv
source .venv/bin/activate
# The `-e .` flag installs the package in "editable" mode.
pip install -e .
```

#### Available commands

You can now run the cli tool and webAPI/server with the following commands.

```bash
divbase-cli
divbase-api # NOTE: We typically do not run the API like this, instead we use docker compose, keep reading below to see how.
```

### 2. Install pre-commit hooks

We also use [pre-commit hooks](https://pre-commit.com/). pre-commit runs on every commit.

```bash
pre-commit install
```

### 3. CLI Enviroment variables

The DivBase project has 2 groups of environment variables:

1. Those needed by the user for running CLI commands and
2. Those needed by the server to interact with itself.

For the CLI commands we have the following environments:

1. local - for local development
2. test - used by pytest
3. scilifelab2dev - to interact with instance of scilifelab2dev cluster.
4. scilifelab2prod (not available yet)

To specify the environment to use you can prepend each command (or set in your shell) the envrioment you want to use:

```bash
DIVBASE_ENV=local divbase-cli files list --project a-local-project
```

**Note:** You do not need to create `.env` files for the `local` (e.g. `.env.local`) or `test`(e.g. `.env.test` ) environments as these settings are not secret. They will be automatically provided when running.

To use the `scilifelab2dev` environment you'll need to create the following file `.env.scilifelab2dev` which should never be committed to source control.

```bash
# .env.scilifelab2dev
DIVBASE_S3_ACCESS_KEY=your_access_key_here
DIVBASE_S3_SECRET_KEY=your_secret_key_here
```

Once created you can run cli commands by passing the environment as shown below.

```bash
DIVBASE_ENV=scilifelab2dev divbase-cli files list --project a-project-in-the-cloud
```

**Note:** Default behaviour if `DIVBASE_ENV` is not set is to get enviroment varialbes from `.env`. This would be used by actual users of the service who will not have to deal with having multiple environments like us.

### 4. Run DivBase backend locally using Docker compose watch

```bash
docker compose -f docker/divbase_compose.yaml down && docker compose -f docker/divbase_compose.yaml watch
```

This will deploy the API, job system and MinIO (S3) instance locally. Using `compose watch` means changes in the `/src` folder will automatically be sycned to the container and if needed trigger a restart or rebuild of them.

Once the backend is up and running you can then run `divbase-cli` commands against it or go to for example <http://localhost:8000/docs> to test out some of API commands directly.

### 5. Running tests

We use docker-compose to setup a testing environment. The testing stack contains a MinIO instance which is populated with some default buckets and data and provided to each test.

```bash
pytest # you may want to append the -s flag to print standard output.
```

**The tests will be slower the first time you run them as the docker images will need to be downloaded and built.**

## Queries

The data and metadata stored in a DivBase project can be queried two main ways:

### 1. Query sidecar metadata file for sample IDs and VCF filenames

Attributes of the samples that are typically not part of the data in the VCF file can be stored in a tabular sidecar file where each row describe a sample. The file need to start with a header row that contains two mandatory columns: `Sample_ID` and `Filename`. Other columns can be added based on the needs of the project. To accomodate queries of user-defined columns, the following query syntax is used:

```bash
divbase-cli query tsv [...] --filter “key1: value1, value2; key2: value3, value4 […]”
```

where `key` is the column name and “value” is the column value to filter on. Values can be stacked by commas, and (intersect) queries across multiple keys can be done by adding a semicolon.

Example:

```bash
query tsv --file sample_metadata.tsv --filter "Area:West of Ireland,Northern Portugal;Sex:F"
```

### 2. Query the data in the VCF files

DivBase uses bcftools to query data contained in the VCF files. One of the idea of DivBase is that data in a project can split over multiple VCF files. If the VCF file query is combined with a sidecar metadata query (using `--tsv-filter` which accepts the same the query syntax the sidecar queries described above), DivBase will perform the query on all VCF files identified by the sidecar query, apply any chain bcftools subset/filter operation given by the user, and return a single merged VCF file with the results.

Bcftools subset/filters are set by the `--command` option. The syntax is experimental: it uses existing bcftools commands as defined in the [bcftools manual](https://samtools.github.io/bcftools/bcftools.html) but automatically handles the filenames based on the sidecar metadata query results. It is possible to create pipes of bcftools operations by semicolon separation:

```bash
--command "view -s SAMPLES; view -r 21:15000000-25000000"
```

Example of a full command that queries the sidecar metadata file and uses the results to perform bcftools subsets:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F" \
--command "view -s SAMPLES; view -r 21:15000000-25000000"
```

## Running query jobs asynchronously using Celery

Queries that require bcftools operations can potentially take long time to complete. Therefore, these queries are submitted to job manager based on Celery.

Unless specified, the bcftools query jobs are run synchronously in the job manager, meaning that they will be executed directly.

It is also possible to submit asyncronous jobs to Celery from the CLI by adding the `--async` flag to `query bcftools-pipe` commands.
For instance, to query of the toy data found in `./tests/fixtures/`:

```bash
divbase-cli query bcftools-pipe --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F" \
--command "view -s SAMPLES; view -r 21:15000000-25000000" --async
```

This will return a Celery task-id to `stdout`. To check the status of the task:

```bash
divbase-cli query task-status --task-id <ID>
```
