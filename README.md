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

## DivBase allows you to

- **Store all your variant data and metadata in one place** - a single, centralised store for all VCF files and sample metadata from your research project.
- **Collaborate and share data** with colleagues and collaborators, with full **control over who has access** to what.
- **Query across all your VCF files at once**, or narrow down to a subset of your choosing.
- **Filter on variant data and sample metadata in the same query** - no need to join results together manually.
- **Integrate the system into your pipelines and HPC jobs** - use DivBase programmatically wherever your workflow needs it.
- **Keep all files under version control and backed up**.
- **Checkpoint the state of your project's files** to refer back to at a later date - **making your research more easily reproducible**.

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

We welcome contributions! See our [Contributing](https://scilifelabdatacentre.github.io/divbase/CONTRIBUTING/) and [Developer Setup](https://scilifelabdatacentre.github.io/divbase/development/setup/) guides to get started. Feel free to reach out if you have questions or aren't sure where to start.
