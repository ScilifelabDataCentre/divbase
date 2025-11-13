
# Using DivBase CLI tool

Please make sure you have followed the install instructions in the main README.md to install divbase-cli before continuing.

The CLI tool relies on 2 local files to preserve state between commands:

1. A config file stored in your home directory at `~/.config/divbase/config.yaml`.
    This config contains information about your DivBase projects and which divbase server (if any) you're logged into.
2. A secrets file stored in your home directory at `~/.config/divbase/.secrets`.
    This file contains your access/refresh tokens for the DivBase server you're logged into.

You do not need to manually modify these files and you should not share your .secrets file with anyone as it contains tokens that allow access to your DivBase projects for up to 1 week. If you have, please contact us.

## 1. Creating your user config file and adding your first project

After installing the divbase-cli package, you can create a user config file

```bash
divbase-cli config create
```

To view the contents of the config file you can run:

```bash
divbase-cli config show
```

To add your first project to the config file you can run:

```bash
divbase-cli config add-project [PROJECT_NAME] # append `--default`  or `-d` to make it the default project
```

**Tip:** On the DivBase website in the project's page, you'll see the command you can run to add that specific project to your user config file.

If you want to change your default project later, you can run the same command with the `--default` flag again or:

```bash
divbase-cli config set-default [PROJECT_NAME] # This must be a project that is already in your config file
```

## 2. Logging into DivBase server

With the config setup you can now log into DivBase server via the CLI tool. This will keep you logged in for 1 week.

```bash
divbase-cli auth login [YOUR_EMAIL] # you'll be prompted for your password or you can pass it via the `--password` flag
# run `divbase-cli auth logout` to log out
```

You will now be able to run commands against your DivBase projects. Here are some examples:

```bash
divbase-cli files list --project my-non-default-project
divbase-cli files download file1.vcf.gz file2.vcf.gz --project my-non-default-project
divbase-cli files upload file1.vcf.gz file2.vcf.gz --project my-non-default-project
```

Notice above how we specified the project name with `--project` flag since it's not our default project. You can do this for all DivBase CLI commands that relate to working with a specific project.

If you're ever unsure which project is the default project you can run:

```bash
divbase-cli config show
# or
divbase-cli config show-default
```

## 3. Versioning the projects state

The DivBase server uses S3 buckets to store project files. Each project has its own assigned bucket. These buckets support versioning of files natively via S3. So individual files can be restored to prior versions if needed.

On top of this we also support versioning of the entire project's state via a special file stored in the bucket. This allows users to create named versions of the entire project at a given timepoint.

One potential use case for this could be:

- Marking a time point when you did analysis for the a publication.

At a later date, you could then download/upload all files from/to the bucket as they were at that timepoint to ensure reproducibility of results.

If you're in a new project, you can create the bucket versioning file by running

```bash
divbase-cli version create
```

This will create a file in the bucket called `.bucket_versions.yaml`

An already existing bucket likely has this file, we can see the contents of this file by running:

```bash
divbase-cli version list
```

If after working with the bucket for a while you want to version the current state of the bucket, you can run:

```bash
divbase-cli version add [OPTIONS] NAME
```

To download files from bucket at a specific bucket version/state, we can use the --bucket-version option and specify the version name we want to download from:

```bash
divbase-cli files download file1.vcf.gz file2.vcf.gz --bucket-version=v0.1.0
```
