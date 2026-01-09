# Developer Setup Guide

## Prerequisites

- Python3.12+
- docker + docker compose
- pre-commit

## 1. Setup your python virtual environment

You can use either [uv](https://github.com/astral-sh/uv) or something like `venv` or `pip`.

### Using `uv` (Recommended)

- Install `uv`, follow the docs: <https://docs.astral.sh/uv/>
- In the root of the repository, run:

```bash
uv sync
```

This will create a virtual environment at the root of the repo and install all dependencies and the package (with all code in the src folder installed in "editable" mode).

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
divbase-api # NOTE: We typically do not run the API like this, instead we use docker compose, keep reading below to see how.
```

## 2. Install pre-commit hooks

We also use [pre-commit hooks](https://pre-commit.com/). pre-commit runs on every commit.

```bash
pre-commit install
```

If you ever need to skip a hook for a specific commit you can use:

```bash
git commit --no-verify
```

Or consider updating the hook config inside `.pre-commit-config.yaml`.

We have a GH action that checks that all pre-commit hooks pass on every PR and commit to main.

## 3. CLI Setup

divbase-cli relies on 2 local files to preserve state between runs:

1. A config file stored in your home directory at `~/.config/divbase/config.yaml`.
    This config contains information about your DivBase projects and which divbase server (if any) you're logged into.
2. A secrets file stored in your home directory at `~/.config/divbase/.secrets`.
    This file contains your access/refresh tokens for the DivBase server you're logged into.

You can read [the user guides](../user-guides/index.md) for more information on how to setup the CLI tool for the first time use.

For convenience, we provide a local development setup script, but in order to run this you need the docker compose stack running (see next section).

## 4. Run DivBase server locally using docker compose watch

```bash
docker compose -f docker/divbase_compose.yaml down && docker compose -f docker/divbase_compose.yaml watch
```

This will deploy the API, job system and a MinIO (S3) instance locally. Using `compose watch` means changes in the `/src` folder will automatically be synced to the container and if needed trigger a restart or rebuild of them.

Once the backend is up and running you can then run `divbase-cli` commands against it or go to for example <http://localhost:8000/api/v1/docs> to test out some of API commands directly.

The DivBase frontend is running on <http://localhost:8000>. The frontend is part of the FastAPI app deployment. All API routes are appended with `/api/v1` to avoid collisions. The frontend uses Jinja2 templating.

To help with setup, we provide a local development setup script which will:

- Create users and projects on your locally running DivBase instance and add users with different roles to each project.
- create a user config file with these projects added to it.
- Upload some initial data to the different projects (sample metadata files and VCFs).

You can run the local_dev_setup script like this:

```bash
uv run scripts/local_dev_setup.py # make sure the compose stack is running
```

## Extras

If using vscode you may want to consider some/all of the following for your `.vscode/settings.json`:

```json
{
  "[python]": {
    "editor.formatOnSave": true,
    "editor.codeActionsOnSave": {
      "source.fixAll": "explicit"
    },
    "editor.defaultFormatter": "charliermarsh.ruff"
  },
  "python.analysis.autoImportCompletions": true,
  "python.testing.unittestEnabled": false,
  "python.testing.pytestEnabled": true
}
```
