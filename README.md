# DivBase: Share and query genetic variant data at scale

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://scilifelabdatacentre.github.io/divbase/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/divbase-cli.svg)](https://pypi.org/project/divbase-cli/)

DivBase is a service that enables life science researchers at Swedish research institutions and their collaborators to easily manage, explore, and query genomic variants in VCF format, along with associated sample metadata.

> [!NOTE]
> This service is built and maintained by [SciLifeLab Data Centre](https://www.scilifelab.se/) and is currently in pre-release for a limited number of users for testing and feedback collection. We are actively seeking feedback from the community to help shape DivBase. If you would like to be involved in testing, have suggestions, or want to contribute, please reach out to us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues).

---

## Table of Contents

1. [What is DivBase?](#what-is-divbase)
2. [Documentation](#documentation)
3. [Install `divbase-cli`](#install-divbase-cli)
4. [(Very) Quick Start Guide](#very-quick-start-guide)
5. [Get support](#get-support)
6. [Developers/Contributing](#developerscontributing)

## What is DivBase?

TODO - Explain what divbase looks like, aka cli, web interface, and how it works at a high level.

![Overview of DivBase Features](docs/assets/img/divbase_key_feats_light.webp#gh-light-mode-only)
![Overview of DivBase Features](docs/assets/img/divbase_key_feats_dark.webp#gh-dark-mode-only)

## Documentation

For full guides, tutorials, and command references, [visit our documentation website](https://scilifelabdatacentre.github.io/divbase/).

## Install `divbase-cli`

To manage files, submit queries, and interact with DivBase, you can use our command line tool `divbase-cli`. The recommended way to install the `divbase-cli` is using [uv](https://docs.astral.sh/uv/) or [pipx](https://pipx.pypa.io/stable/). To install with uv or pipx, run:

```bash
uv tool install divbase-cli
# or with pipx
pipx install divbase-cli
```

For detailed instructions and alternative methods, see the [Installation Guide](https://scilifelabdatacentre.github.io/divbase/user-guides/installation/).

## (Very) Quick Start Guide

> [!TIP]
> Go to our documentation website for a proper quick start guide [Quick Start Guide](https://scilifelabdatacentre.github.io/divbase/user-guides/quick-start/).

1. **Create an Account:** Sign up at the [DivBase Web Interface](https://divbase.scilifelab-2-prod.sys.kth.se).
2. **Configure:** Add your project to the CLI:

   ```bash
   divbase-cli config add MY_PROJECT_NAME --default
   ```

3. **Login:**

   ```bash
   divbase-cli auth login your.email@example.com
   ```

4. **Upload Data:**

   ```bash
   divbase-cli files upload data/my_samples.vcf.gz
   ```

5. **Run a Query:**

   ```bash
   # Subset a chromosomal region
   divbase-cli query vcf --command "view -r 21:15000000-25000000"
   ```

## Get support

- **Need help?** Contact us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues).

- **Have feedback or spotted a bug?** Please report any issues or suggestions on our [GitHub Issues page](https://github.com/ScilifelabDataCentre/divbase/issues).

We welcome contributions! Please see our [Contributing Guide](docs/CONTRIBUTING.md) and [Security Policy](docs/SECURITY.md) to get started.

---

## Developers/Contributing

We welcome contributions! Feel free to checkout our contributing guide and developer setup guide to get started. You are welcome to reach out to us if you have any questions or want to contribute but are not sure where to start.

See the [Developer Setup Guide](docs/development/setup.md) (or view the developer documentation live at <https://scilifelabdatacentre.github.io/divbase/development/setup/>) to get started.
