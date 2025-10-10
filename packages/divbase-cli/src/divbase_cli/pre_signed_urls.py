"""
Helper functions to handle using pre-signed URLs for file download and upload.

TODO: Consider adding retries, error handling, progress bars, etc.
"""

from pathlib import Path

import httpx


def download_multiple_pre_signed_urls(pre_signed_urls: list[dict], download_dir: Path) -> list[Path]:
    """
    Download files using pre-signed URLs.
    Returns a list of Paths to the downloaded files.
    """
    downloaded_files = []
    with httpx.Client(timeout=30.0) as client:
        for obj in pre_signed_urls:
            object_name = obj["object_name"]
            pre_signed_url = obj["pre_signed_url"]
            out_file_path = download_dir / object_name

            with client.stream("GET", pre_signed_url) as response:
                response.raise_for_status()

                with open(out_file_path, "wb") as file:
                    for chunk in response.iter_bytes(chunk_size=8192):
                        file.write(chunk)

            downloaded_files.append(out_file_path)

    return downloaded_files


def upload_multiple_pre_signed_urls(pre_signed_urls: list[dict], all_files: list[Path]) -> dict[str, Path]:
    """
    Upload files using pre-signed POST URLs.
    Returns a dict mapping the uploaded file names to their Paths.
    """
    file_map = {file.name: file for file in all_files}
    uploaded_files = {}

    with httpx.Client(timeout=30.0) as client:
        for obj in pre_signed_urls:
            object_name = obj["object_name"]
            post_url = obj["post_url"]
            fields = obj["fields"]

            file_path = file_map[object_name]
            with open(file_path, "rb") as file:
                files = {"file": (object_name, file, "application/octet-stream")}
                response = client.post(post_url, data=fields, files=files)
                response.raise_for_status()

            uploaded_files[object_name] = file_path

    return uploaded_files
