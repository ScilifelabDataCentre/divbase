
# Setup DivBase CLI configuration guide

This section covers how to configure divbase-cli to know which projects you're a member of and log you in and out of the DivBase server.

## Prerequisites

- Installed divbase-cli, see the [quickstart guide](quick-start.md) or [installation section](installation.md) of the documentation.
- Have created a [DivBase account via the website](https://divbase.scilifelab.se) and have access to at least one DivBase project. See the [account management guide](account-management.md) for more information on how to get access to a project.

## 1. Create a user config file and add your project(s)

After installing `divbase-cli`, you can add a new project to your user config file by running:

```bash
divbase-cli config add [PROJECT_NAME] # append `--default` to make it the default project
```

!!! info
    On the DivBase website in the project's page, you'll see the command you can run to add that specific project to your user config file.

If you want to change your default project later, you can run the same command with the `--default` flag again or:

```bash
divbase-cli config set-default [PROJECT_NAME] # This must be a project that is already in your config file
```

To view the contents of your user config file, you can run:

```bash
divbase-cli config show
```

## 2. Logging into DivBase server

With the config setup you can now log into DivBase server via `divbase-cli`. This will keep you logged in for 1 week.

```bash
divbase-cli auth login [YOUR_EMAIL] # you'll be prompted for your password
# divbase-cli auth logout  # run to log out
```

You will now be able to run commands against your DivBase projects. Here are some examples:

```bash
divbase-cli files ls # runs against the default project
divbase-cli files download file1.vcf.gz file2.vcf.gz --project my-non-default-project
divbase-cli files upload file1.vcf.gz file2.vcf.gz --project my-non-default-project
```

Notice above how we sometimes specified the project name with `--project` flag since we wanted to use a non-default project. You can do this for all DivBase CLI commands that relate to working with a specific project.

If you're ever unsure which project is the default project you can run:

```bash
divbase-cli config show
# or
divbase-cli config show-default
```

---

## How divbase-cli stores state

(You don't need to read this section unless you're curious!)

`divbase-cli` relies on 2 files stored locally to preserve state between commands/sessions:

1. A config file stored in your home directory at `~/.config/divbase/config.yaml`.
    This config contains information about your DivBase projects and which divbase server (if any) you're logged into.
2. A secrets file stored in your home directory at `~/.config/divbase/.secrets`.
    This file contains your personal tokens for the DivBase server you're logged into that grant access for up to 1 week. You should never share this file with anyone.

*You do not need to manually modify these files. The commands above will create and update them for you as needed.*

!!! info "Access from multiple workstations"
    If you plan to access divbase from say both your laptop and HPC, you will need to install `divbase-cli` in both places and set up the config and login in both places separately.
