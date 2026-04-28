# Installing divbase-cli

## Recommended: uv

The recommended way to install `divbase-cli` is via [uv](https://docs.astral.sh/uv/){target="_blank"}, a fast Python package manager. uv manages isolated tool environments automatically, so there are no dependency conflicts and no need to manually create/activate a virtual environment

**Step 1 — Install uv** (if you don't have it already):

=== "Linux / macOS / Windows Subsystem for Linux"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

=== "Windows"

    ```powershell
    powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
    ```

See the [uv installation docs](https://docs.astral.sh/uv/getting-started/installation/){target="_blank"} for more details.

**Step 2 — Install divbase-cli:**

```bash
uv tool install divbase-cli
```

`divbase-cli` is now available in your path and you can start using it immediately.

To upgrade later:

```bash
uv tool upgrade divbase-cli
```

!!! info "Which version is installed?"
    To check which version of divbase-cli is currently installed, run:

    ```bash
    divbase-cli --version
    ```

## Alternative: pipx

If you already use [pipx](https://pipx.pypa.io/stable/), it works just as well:

```bash
pipx install divbase-cli
```

To upgrade:

```bash
pipx upgrade divbase-cli
```

---

!!! question "Why not create a virtual environment?"
    `divbase-cli` is a command-line tool, not a library. Installing it into a virtual environment (venv, conda, mamba etc.) is unnecessary and means you have to activate that environment every time you want to use it. Use `uv tool` or `pipx` instead. The benefits of isolation from creating a virtual environment are already provided by these tools.

    If you really want to install divbase-cli into a virtual environment you can do so with pip
    ```bash
    # make sure you are using python 3.12 or higher
    pip install divbase-cli
    ```
