# Installing divbase-cli

## Recommended options: pipx or uv tool

### 1. With pipx

```bash
pipx install divbase-cli
```

If you do not have pipx installed, you can install it by following [the official instructions from pipx](https://pipx.pypa.io/stable/installation/).

!!! info "why pipx and not pip?"
    Pipx is recommended for command-line tools that need to be run from the terminal. With pipx, the tool is installed in an isolated environment and made available globally, which avoids dependency conflicts with other Python packages on your system.

To upgrade later, use:

```bash
pipx upgrade divbase-cli
```

### 2. With uv tool

If you use the package manager [uv](https://docs.astral.sh/uv/), you can install `divbase-cli` using UV's equivalent to pipx (`uv tool`):

```bash
uv tool install divbase-cli
uv tool upgrade divbase-cli
```

## Alternative Installation Methods

### pip with a Virtual Environment

If you would prefer to manually create a virtual enviroment with for example `conda`/`mamba`/`venv` then you can install `divbase-cli` with pip after activating your virtual enviroment.

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

## Upgrading divbase-cli

To check your currently installed version:

```bash
divbase-cli --version
```

To upgrade to the latest version:

**With pipx (recommended):**

```bash
pipx upgrade divbase-cli
```

**With uv (also recommended):**

```bash
uv tool upgrade divbase-cli
```

**With pip:**

```bash
pip install --upgrade divbase-cli
```
