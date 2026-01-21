"""
Constants that both divbase-api and divbase-cli need to agree on.
"""

ONE_MiB = 1024 * 1024

# When you download a file that has been uploaded in parts, you have
# to know the part/chunk size used in order to correctly calculate the composite checksum
S3_MULTIPART_CHUNK_SIZE = 32 * ONE_MiB

# At what point you swap from single part to multipart upload to S3.
# If server and client used the same threshold then makes life easier
# when validating the checksums of files in s3 as single part and multipart uploads use different ETag formats.
# (No benefit in constraining the download threshold, so not done here)
S3_MULTIPART_UPLOAD_THRESHOLD = 96 * ONE_MiB

# Max number of items that can be processed in a single API call to divbase-api's S3 routes
# covers e.g pre-signed urls for upload/download, soft delete and checksum comparisons
# client has to batch requests if exceeding this limit
MAX_S3_API_BATCH_SIZE = 100

# How long the pre-signed URLs divbase-api creates are valid for
SINGLE_PART_UPLOAD_URL_EXPIRATION_SECONDS = 3600  # 1 hour
MULTI_PART_UPLOAD_URL_EXPIRATION_SECONDS = 36000  # 10 hours
DOWNLOAD_URL_EXPIRATION_SECONDS = 36000  # 10 hours

# (Not used anywhere, just making it explicit)
# This is limited by our fixing of the chunk size and S3's limit to the number of chunks allowed (10,000)
# 320 GiB if using 32 MiB chunks
LARGEST_FILE_UPLOADABLE_TO_DIVBASE_BYTES = 10_000 * S3_MULTIPART_CHUNK_SIZE
