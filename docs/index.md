# DivBase

!!! warning "DivBase is in pre-release stage and under active development. Expect changes to features, commands, and documentation before the official release. If you want to be involved in testing or have suggestions, please reach out at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues)."

---

DivBase is a service built and maintained by [SciLifeLab Data Centre](https://www.scilifelab.se/) that enables life science researchers at Swedish institutions and their collaborators to manage, explore, and query genomic variants in VCF format alongside associated sample metadata. The service provides a secure platform for managing genomic variants and metadata files for non-human and non-sensitive data.

## Want to access DivBase?

- **Join an existing project:** Create an account at the [DivBase web interface](https://divbase.scilifelab-2-prod.sys.kth.se){:target="_blank"} and ask the project manager to add you (tell them the email you used to sign up with).
- **Start your own project:** contact us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se).

## Key Features

![Overview of DivBase Features](assets/img/divbase_key_feats_light.webp#only-light)
![Overview of DivBase Features](assets/img/divbase_key_feats_dark.webp#only-dark)

## DivBase allows you to

- **Store all your variant data and metadata in one place** - a single, centralised store for all VCF files and sample metadata from your research project.
- **Collaborate and share data** with colleagues and collaborators, with full **control over who has access** to what.
- **Query across all your VCF files at once**, or narrow down to a subset of your choosing.
- **Filter on variant data and sample metadata in the same query** - no need to join results together manually.
- **Integrate the system into your pipelines and HPC jobs** - use DivBase programmatically wherever your workflow needs it.
- **Keep all files under version control and backed up**.
- **Checkpoint the state of your project's files** to refer back to at a later date - **making your research more easily reproducible**.

## Install the CLI tool

Interact with DivBase using the `divbase-cli` command line tool. Install it with [uv](https://docs.astral.sh/uv/) or [pipx](https://pipx.pypa.io/stable/):

```bash
uv tool install divbase-cli
# or with pipx
pipx install divbase-cli
```

See the [Installation Guide](user-guides/installation.md) for detailed instructions and alternative methods.

## Getting Started

See our [Quick Start Guide](user-guides/quick-start.md) to get up and running in minutes or follow our tutorial on [Running a query on a public dataset](user-guides/tutorial-query-on-public-data.md).

## Getting Support

- **Need help?** Contact us at [dsn-eb@scilifelab.se](mailto:dsn-eb@scilifelab.se) or open a [GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues).
- **Found a bug?** Please report it on our [GitHub Issues page](https://github.com/ScilifelabDataCentre/divbase/issues).
