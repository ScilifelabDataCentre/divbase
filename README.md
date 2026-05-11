# DivBase: Share and query genetic variant data at scale

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://scilifelabdatacentre.github.io/divbase/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/divbase-cli.svg)](https://pypi.org/project/divbase-cli/)

DivBase is a service built and maintained by [SciLifeLab Data Centre](https://www.scilifelab.se/) that enables life science researchers at Swedish institutions and their collaborators to manage, explore, and query genomic variants in VCF format alongside associated sample metadata. The service provides a secure platform for managing genomic variants and metadata files for non-human and non-sensitive data.

> [!NOTE]
> DivBase is currently in pre-release for a limited number of users. We are actively seeking feedback to help shape the service. If you would like to be involved in testing or have suggestions, please reach out at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues).

---

## Want to try out DivBase?

- **Join an existing project:** Create an account at the [DivBase web interface](https://divbase.scilifelab-2-prod.sys.kth.se) and ask your project manager to add you.
- **Start your own project:** Reach out to us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se).

## Key Features

![Overview of DivBase Features](docs/assets/img/divbase_key_feats_light.webp#gh-light-mode-only)
![Overview of DivBase Features](docs/assets/img/divbase_key_feats_dark.webp#gh-dark-mode-only)

## Why use DivBase?

- A **single, centralised store** of your project's variant data and metadata.
- Easy to **collaborate and share data** with your colleagues and collaborators - and **control who has access** to what.
- Queries are run on **all (or a subset of your choosing) VCF files** stored in the project.
- Possible to **filter both on variant data and sample metadata** in the same query.
- You can use **DivBase programmatically** in for example **pipelines/HPC jobs**.
- Files are **versioned and backed up**.
- You can **version/checkpoint the state of your project's files** to refer back to at a later date - **making your research more easily reproducible**.

## Documentation

For guides, tutorials, and command references, [visit our documentation website](https://scilifelabdatacentre.github.io/divbase/).

## Install `divbase-cli`

To manage files, submit queries, and interact with DivBase, install our command line tool `divbase-cli` using [uv](https://docs.astral.sh/uv/) or [pipx](https://pipx.pypa.io/stable/):

```bash
uv tool install divbase-cli
# or with pipx
pipx install divbase-cli
```

For detailed instructions and alternative methods, see our [Installation Guide](https://scilifelabdatacentre.github.io/divbase/user-guides/installation/).

## Quick Start Guide

> [!TIP]
> Go to our documentation website for the proper [Quick Start Guide](https://scilifelabdatacentre.github.io/divbase/user-guides/quick-start/).

1. **Create an Account:** Sign up at the [DivBase Web Interface](https://divbase.scilifelab-2-prod.sys.kth.se).

2. **Configure:** Add your project to your CLI config and set it as default:

   ```bash
   divbase-cli config add MY_PROJECT_NAME --default
   ```

3. **Login:**

   ```bash
   divbase-cli auth login your.email@example.com
   ```

4. **Upload Data and sync your data**

   ```bash
   divbase-cli files upload data/my_samples.vcf.gz
   divbase-cli files upload --upload-dir data/
   divbase-cli dimensions update
   ```

5. **Run a Query:**

   ```bash
   # Filter samples based on metadata and subset a chromosomal region in one
   divbase-cli query vcf \
      --tsv-filter "Area:Northern Portugal" \
      --command "view -r 21:15000000-25000000"
   ```

This will submit a job to DivBase and once the job is complete, a new vcf.gz file containing the subset of data you requested will be available for download/streaming in your downstream analysis.

## Get Support

- **Need help?** Contact us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues).
- **Found a bug?** Please report it on our [GitHub Issues page](https://github.com/ScilifelabDataCentre/divbase/issues).

---

## Developers/Contributing

We welcome contributions! See our [Contributing Guide](docs/CONTRIBUTING.md) and [Developer Setup Guide](docs/development/setup.md) to get started, or view the [developer documentation](https://scilifelabdatacentre.github.io/divbase/development/setup/) online. Feel free to reach out if you have questions or aren't sure where to start.
