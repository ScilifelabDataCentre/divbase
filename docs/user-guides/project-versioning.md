# Managing Project Versions

DivBase allows you to create named versions of your project's state at specific points in time. This is useful for tracking changes, ensuring reproducibility, and marking important milestones (e.g., when running analysis for a publication).

This means you can save the current state of your project as a named version, and later retrieve files as they were at that specific version (so you don't have to worry about having later updated files).

!!! info "What actually is a version entry?"
    A version entry is stored on our server and alongside some metadata (like name, description, created at) it contains references to every file in your project at that time point like this:

    ```json
    files = {
        "file1.vcf.gz": "unique_file_hash_12345",
        "file2.vcf.gz": "unique_file_hash_67890",
        ...
    }
    ```
    The unique files hashes mean uploading a new version of the file with the same name doesn't matter. DivBase will now that entry refers to the older version of the file.

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
divbase-cli version list # --project [PROJECT_NAME] to specify a non-default project
```

This will display a list of all saved versions, including their names, descriptions, and creation dates.

## Specific version details

To get detailed information about a specific version, use:

```bash
divbase-cli version info NAME
```

This view will also include all the files associated with that version.

## Deleting versions

To delete a specific version from your project, run:

```bash
divbase-cli version delete NAME
```

This will delete the version entry from the project. Deleted versions older than 30 days will be permanently deleted. You can ask a DivBase admin to restore a deleted version within that time period.

**NOTE:** The files associated with the version are never deleted by this operation.

## Downloading files from a specific version

To download files from your project as they were at a specific version, use the --project-version option:

```bash
divbase-cli files download file1.vcf.gz file2.vcf.gz --project-version=NAME
```

Replacing `NAME` with the name of the version you want to download files from.
