# Installing divbase-cli

## Recommended options: pipx or uv tool

### 1. With pipx

```bash
pipx install divbase-cli
```

If you do not have pipx installed, you can install it by following [the official instructions from pipx](https://pipx.pypa.io/stable/installation/).

!!! info "why pipx and not pip?"
    Pipx is recommended for command-line tools that need to be run from the terminal. With pipx, the tool is installed in an isolated environment and made available globally, which avoids dependency conflicts with other Python packages on your system.

To upgrade your prior install of divbase-cli, use:

```bash
pipx upgrade divbase-cli
```

!!! info "Which version is installed?"
    To check which version of divbase-cli is currently installed, run:

    ```bash
    divbase-cli --version
    ```

### 2. With uv tool

If you use the package manager [uv](https://docs.astral.sh/uv/), you can install `divbase-cli` using UV's equivalent to pipx (`uv tool`):

```bash
uv tool install divbase-cli
```

And to upgrade an existing installation:

```bash
uv tool upgrade divbase-cli
```

## Alternative Installation Methods

### pip with a Virtual Environment

If you would prefer to manually create a virtual environment with for example `conda`/`mamba`/`venv` then you can install `divbase-cli` with pip after activating your virtual environment.

For example, with conda:

```bash
conda create -n divbase python=3.13
conda activate divbase
pip install divbase-cli
```

And to upgrade later:

```bash
conda activate divbase
pip install --upgrade divbase-cli
```
