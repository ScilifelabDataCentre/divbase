# S3 File Transfers

Each DivBase project is assigned it's own S3 bucket and divbase interacts with S3 object storage in two primary ways:

1. Through `divbase-cli`. A user receives pre-signed URLs from the DivBase API to upload and download files directly to/from S3.

2. Through the `divbase-api` celery worker. The worker has a service account and make use of `boto3` library to interact with S3.

Both can handle single part and multipart transfers.

## Single vs Multipart Transfers

To handle large files efficiently (and allow file sizes > 5 GB), DivBase implements the S3 multipart transfer protocol. A file is considered "large" if it exceeds the `S3_MULTIPART_UPLOAD_THRESHOLD` (currently set to 96 MiB), defined in `divbase-lib`.

Unlike singlepart uploads, where the file is uploaded in one go, multipart uploads break the file into smaller parts/chunks and uploads them independently. The boto3 library does this for you, when using pre-signed URLs, the process is more involved.

The steps in a multipart upload for divbase-cli are approximately as follows:

1. Initiate the multipart upload with S3 to get an `upload_id`.
2. Split the file into parts of size `S3_MULTIPART_CHUNK_SIZE`.
3. For each part:
    - Generate a pre-signed URL for that part using the `upload_id` and part
    - Upload the part to S3 using the pre-signed URL.
4. After all parts are uploaded, complete the multipart upload by informing S3 to assemble the parts into the final object.

## Checksum Validation (ETag)

DivBase server performs integrity checks on files uploaded/download to/from the worker via MD5 checksums. The CLI has this functionality too, but a user can turn off checksum validation with a flag.

Essentially, you compare the locally computed checksum of the file against the checksum provided by S3 (the `ETag` in the files metadata).

!!! Info "The format of the `ETag` depends on how the file was uploaded."
    -   **Single-Part Uploads:** For smaller files uploaded in one go, the `ETag` is simply the **hex-encoded MD5 hash** of the file.

    -   **Multipart Uploads:** For large files, the `ETag` is a **composite checksum**. It is calculated by taking the MD5 hash of each individual part, concatenating them, and then taking the MD5 hash of that combined string. The final `ETag` is this composite hash followed by a hyphen and the total number of parts (e.g., `a1b2c3d4-56`). So 56 parts were used in the upload.

The `divbase-lib` contains the necessary logic (`s3_checksums.py`) to correctly calculate and verify both types of checksums. If a downloaded file's checksum does not match the expected `ETag`, the operation fails, and the corrupted local file is deleted.

!!! Warning "Multipart Uploads Validation Gotcha"
    - The above means that for multipart uploads, you need to know the part/chunk size to correctly compute the checksum. If the part size used during upload is different to what you expect it to be when you download it, then your computed checksum will not match the `ETag`, leading to validation failures. Therefore, it's crucial to use the same `S3_MULTIPART_CHUNK_SIZE` for all file uploads.

    - This concern does not apply to downloads (can use whatever chunk you want), only when you actually calculate the etag for the downloaded file do you need to know the chunk size.

    - To make sure we do this consistently, we set some shared constants, see below...

### Shared Constants

To avoid this gotcha DivBase enforces the same:

- `S3_MULTIPART_UPLOAD_THRESHOLD`: The file size (in bytes) above which multipart transfers are triggered.
- `S3_MULTIPART_CHUNK_SIZE`: The size (in bytes) of each part for an upload.

These are defined in the shared `divbase-lib` package to ensure consistency between the CLI and the worker (`packages/divbase-lib/src/divbase_lib/api_schemas/divbase_constants.py`).

## Concurrency and Retries

- **Concurrency:** Both the CLI and the worker use thread pools to perform multipart transfers in parallel.
- **Retries:** For the CLI, S3 operations have retry logic (using the `stamina` library) to automatically retry downloading /uploading file parts that fail due to transient issues (e.g., network connection problems). The worker uses built-in `boto3` retry mechanisms.
