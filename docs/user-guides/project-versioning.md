# Managing Project Versions

DivBase allows you to create named versions of your project's state at specific points in time. This is useful for tracking changes, ensuring reproducibility, and marking important milestones (e.g., when running analysis for a publication).

This means you can save the current state of your project as a named version, and later retrieve files as they were at that specific version (so you don't have to worry about having later updated the files).

Files that are part of a project version are protected from hard deletion, to ensure you can restore them at any point in the future.

!!! info "What actually is a (project) version entry?"
    A version entry is stored on our server and alongside some metadata (like name, description, created at) it contains references to every file in your project at that time point like this:

    ```json
    "file_name1.vcf.gz": {
        "version_id": "unique_file_identifier_12345",
        "etag": "md5 checksum of file content",
        "size": file_size_in_bytes
    },
    ```
    The unique file identifier means uploading a new version of the file with the same name doesn't matter. DivBase will know that the entry refers to the older version of the file.
---

## Add a new version

To save the current state of your project run:

```bash
divbase-cli version add NAME --description "[OPTIONAL DESCRIPTION]"
```

- Replace `NAME` with a unique name for the version (e.g., v1.0.0).
- You can optionally add a description using the --description flag.

**Note:** Versions are project wide, so you share them with all other members of the same project.

## Listing versions

To see all existing versions for your project, run:

```bash
divbase-cli version ls # --project [PROJECT_NAME] to specify a non-default project
```

This will display a list of all saved versions, including their names, descriptions, and creation dates.

## Specific version details

To get detailed information about a specific version, use:

```bash
divbase-cli version info NAME
```

This view will also include all the files associated with that version.

## Update a version

You can change the name and/or description of an existing project version, via the update subcommand:

```bash
divbase-cli version update NAME [OPTIONS]
```

- Replace `NAME` with the current name of the version you want to modify.
- Use `--new-name` (or `-n`) to rename the version.
- Use `--new-description` (or `-d`) to update the description.

At least one of `--new-name` or `--new-description` must be provided.

Example:

```bash
divbase-cli version update "my old name" \
--new-name "spruce pine study" \
--new-description "Represents state of data prior to running pangenome analysis for paper X"
```

!!! Info "Only the name and the description can be modified"
    The files and timestamp associated with a version entry are immutable and cannot be changed.
    This is by design to ensure version entries are representations of the project's state at the timepoint they were created at.
    If you need to capture a new state, create a new version instead with `divbase-cli version add`.

## Deleting versions

To delete a specific version from your project, run:

```bash
divbase-cli version rm NAME
```

This will delete the version entry from the project. Deleted versions older than 30 days will be permanently deleted. You can ask a DivBase admin to restore a deleted version within that time period.

**NOTE:** The files associated with the version are never deleted by this operation.

## Downloading files from a specific version

To download files from your project as they were at a specific version, use the --project-version option:

```bash
divbase-cli files download file1.vcf.gz file2.vcf.gz --project-version=NAME
```

Replacing `NAME` with the name of the version you want to download files from.

## Download all files from a specific version

To download all files from a specific project version, use the `files download-all` command and specify the project version:

```bash
divbase-cli files download-all --project-version=NAME
```

!!! Info "`files download-all` tips"
    Three useful flags that you may want to use:

    1. `--download-dir` - Where to download to
    1. `--resume` - Skip already download files (same name & MD5 checksum) that are already downloaded in the download-dir.
    2. `--dry-run` - see what would be downloaded, but don't actually perform the download.

    More information about the `files download-all` command can be found in the [files download-all documentation](./files.md#downloading-all-files).
