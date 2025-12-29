# Developer Overview

## Folders overview

The DivBase codebase is structured as a monorepo with the following main folders:

- **`adr/`** - Architecture Decision Records
- **`docker/`** - dockerfiles and docker compose files for local development and testing (k8s used in production).
- **`docs/`** - Documentation
- **`packages/`** - Source code for DivBase packages (including CLI, API and worker and a shared library)
- **`scripts/`** - Helper scripts for development tasks (e.g. generating mock data for local development, running benchmarks, building docs).
- **`tests/`** - Tests for all DivBase packages

## Architecture overview

DivBase is structured as a monorepo and contains 2 packages for end use, and 1 library module:

1. `divbase-cli`: A command-line interface for users to interact with a DivBase instance built using [Typer](https://typer.tiangolo.com).
2. `divbase-api`: An API (built with [FastAPI](https://fastapi.tiangolo.com/)) deployed on the server which handles all communication between the user, job system and S3.
    - `divbase-worker` (a submodule called "worker") A [celery worker](https://docs.celeryq.dev/en/stable/) deployed on the server that performs long running tasks like queries on the VCF files and associated sample metadata.
3. `divbase-lib`: A library that allows us to share code between the CLI and API packages. One example use case is the API schemas used for the request and response models between the API and CLI.

The API package has 2 entrypoints (fastapi and celery worker) and the CLI package has a single entrypoint (the command line interface).

We use docker compose in local development and testing to create a matching local environment with all the services we need. Taking a look at the `docker/divbase_compose.yaml` file is a good way to get an overview of the different services that make up DivBase and how they interact.

You can also take a look at our architecture decision records (ADRs) folder - `adr/` (and their associated PRs) for some more context surronding choices we have made over time.

## Some good to know things

1. Whilst we expose an API, all user interactions with DivBase are intended to be done via `divbase-cli`. The CLI handles authentication, configuration management and making the appropriate API calls to divbase.

2. In production/deployed enviroments, we have an S3 tenant provided by KTH NetApp. We mimic this in local development/testing with a [MinIO](https://min.io/) instance running in the docker compose stack.

3. User data is stored in S3 buckets, one bucket for one project. The metadata about the projects, users and jobs is stored in a Postgres database. A user should not need to understand what a bucket is in order to use DivBase.
