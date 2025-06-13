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

First setup your python virtual environment.

You can use either [uv](https://github.com/astral-sh/uv) or something like `venv` or `pip`.

### Using `uv` (Recommended)

- Install `uv`, follow the docs: https://docs.astral.sh/uv/

- In the root of the repository, run:

```bash
uv sync
```

This will create a virtual environment at the root of the repo and install all dependacies and the package (with all code in the src folder installed in "editable" mode).

### Using pip and venv

```bash
python3 -m venv .venv
source .venv/bin/activate
# The `-e .` flag installs the package in "editable" mode.
pip install -e .
```

### Available commands

You can now run the cli tool and webAPI/server with the following commands.

```bash
divbase-cli
divbase-api
```

### Pre-commit

We also use [pre-commit hooks](https://pre-commit.com/). pre-commit runs on every commit. 

```bash
pre-commit install
```