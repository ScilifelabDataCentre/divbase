
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
divbase-cli config add [PROJECT_NAME] # append `--default`  or `-d` to make it the default project
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

DivBase allows you to create named versions of your project's state at specific points in time. This is useful for tracking changes, ensuring reproducibility, and marking important milestones (e.g., when running analysis for a publication).

This means you can save the current state of your project as a named version, and later retrieve files as they were at that specific version (so you don't have to worry about having later updated files).

### Add a new version

To save the current state of your project run:

```bash
divbase-cli version add NAME [OPTIONS]
```

- Replace `NAME` with a unique name for the version (e.g., v1.0.0).
- You can optionally add a description using the --description flag.

**Note:** Versions are project wide, so you share them with all other members of the same project.

### Listing Versions

To see all existing versions of your project, run:

```bash
divbase-cli version list
```

This will display a list of all saved versions, including their names, descriptions, and creation dates.

### Specific Version Details

To get detailed information about a specific version, use:

```bash
divbase-cli version info NAME
```

This view will also include all the files associated with that version.

### Deleting Versions

To delete a specific version from your project, run:

```bash
divbase-cli version delete NAME
```

This will delete the version entry from the project. Deleted versions older than 30 days will be permanently deleted. You can ask a DivBase admin to restore a deleted version within that time period.

**NOTE:** The files associated with the version are never deleted by this operation.

### Downloading files from a specific version

To download files from your project as they were at a specific version, use the --project-version option:

```bash
divbase-cli files download file1.vcf.gz file2.vcf.gz --project-version=NAME
```

Replacing `NAME` with the name of the version you want to download files from.
