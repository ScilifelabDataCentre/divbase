# Shared bcftools helpers used by both VCF dimension indexing and query orchestration in DivBase.

import logging
import os
import shlex
import subprocess
from pathlib import Path

from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    TaskUserError,
)

logger = logging.getLogger(__name__)

BCFTOOLS_CONTAINER_NAME = (
    "divbase-worker-quick-1"  # for synchronous tasks in local dev, use this container name to find the container ID
)


def _is_in_kubernetes() -> bool:
    return "KUBERNETES_SERVICE_HOST" in os.environ


def get_container_id(container_name: str) -> str:
    """
    Return the container ID of a running Docker container by name.
    Raises BcftoolsEnvironmentError if the container is not found or the command fails.
    """
    try:
        result = subprocess.run(
            ["docker", "ps", "--filter", f"name={container_name}", "--format", "{{.ID}}"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except subprocess.SubprocessError as e:
        logger.error(f"Docker command failed: {e}")
        raise BcftoolsEnvironmentError(container_name) from e


def run_bcftools(command: str, capture_output: bool = False, capture_stderr: bool = False) -> subprocess.Popen:
    """
    Run a bcftools command, routing through docker exec when called outside a Docker/k8s container.
    Returns a subprocess.Popen so the metrics server can monitor the process or call communicate().
    """
    logger.info(f"Running: bcftools {command}")
    try:
        command_args = shlex.split(command)
    except ValueError as e:
        raise BcftoolsCommandError(
            command=command, error_details=f"Could not parse bcftools command arguments: {e}"
        ) from None

    if not command_args:
        raise BcftoolsCommandError(command=command, error_details="Empty bcftools command after parsing") from None

    in_docker = os.path.exists("/.dockerenv")
    in_k8s = _is_in_kubernetes()
    if capture_output:
        popen_kwargs = {"stdout": subprocess.PIPE, "stderr": subprocess.PIPE, "text": True}
    elif capture_stderr:
        popen_kwargs = {"stdout": subprocess.DEVNULL, "stderr": subprocess.PIPE, "text": True}
    else:
        popen_kwargs = {}

    if in_docker or in_k8s:
        logger.debug("Running inside Celery worker Docker container, executing bcftools directly")
        try:
            return subprocess.Popen(["bcftools"] + command_args, **popen_kwargs)
        except Exception as e:
            logger.error(f"Failed to run bcftools directly: {e}")
            raise BcftoolsCommandError(command=command, error_details=str(e)) from e
    else:
        try:
            container_id = get_container_id(BCFTOOLS_CONTAINER_NAME)
            logger.debug(f"Executing command in container with ID: {container_id}")
            docker_cmd = ["docker", "exec", container_id, "bcftools"] + command_args
            return subprocess.Popen(docker_cmd, **popen_kwargs)
        except BcftoolsEnvironmentError:
            raise
        except Exception as e:
            logger.error(f"Failed to run bcftools in container: {e}")
            raise BcftoolsCommandError(command=command, error_details=str(e)) from e


def _raise_task_user_error_from_bcftools_stderr(stderr: str, operation: str, target: str) -> None:
    """
    Raise user-facing TaskUserError for known/important bcftools stderr patterns.

    This centralizes bcftools stderr classification so both dimensions indexing and query orchestration
    can produce consistent user-facing errors.
    """
    stderr_clean = (stderr or "").strip()
    if not stderr_clean:
        return

    stderr_lower = stderr_clean.lower()

    if "unsorted positions" in stderr_lower:
        raise TaskUserError(
            f"{target} is not sorted by position and cannot be indexed by bcftools.\n"
            "DivBase requires VCF files to be sorted by position per scaffold for bcftools orchestration.\n"
            "Please sort the file (for example with 'bcftools sort'), upload the file to the DivBase project, and submit the query again."
        ) from None

    if "duplicated sample name" in stderr_lower:
        raise TaskUserError(
            f"{target} contains duplicate sample IDs in the header and is not valid for DivBase.\n"
            "Please ensure all sample names in the file header are unique, re-upload the corrected file, and run dimensions update again."
        ) from None

    if "could not parse header" in stderr_lower:
        raise TaskUserError(
            f"bcftools failed while {operation} for {target} due to an invalid VCF header.\n"
            f"Details from bcftools:\n{stderr_clean}"
        ) from None

    if "unknown file type" in stderr_lower:
        raise TaskUserError(
            f"{target} has invalid VCF header/content and could not be parsed by bcftools.\n"
            "Please validate the file format/header, upload a corrected VCF, and run dimensions update again."
        ) from None

    if (
        "wrong number of fields" in stderr_lower
        or "could not parse the line" in stderr_lower
        or ("number of columns" in stderr_lower and "does not match the number of samples" in stderr_lower)
    ):
        raise TaskUserError(
            f"{target} contains malformed VCF record lines and cannot be processed by bcftools.\n"
            f"Details from bcftools:\n{stderr_clean}"
        ) from None

    if "not compressed with bgzip" in stderr_lower or "cannot be usefully indexed" in stderr_lower:
        raise TaskUserError(
            f"{target} is gzip-compressed but not bgzip-compressed, so bcftools cannot create an index.\n"
            "Please re-compress the file with bgzip, upload the corrected file to the DivBase project, and run dimensions update again."
        ) from None

    if (
        "truncated" in stderr_lower
        or ("bgzf" in stderr_lower and "failed to read" in stderr_lower)
        or "eof marker is absent" in stderr_lower
    ):
        raise TaskUserError(
            f"{target} appears to be corrupted or truncated and cannot be processed by bcftools.\n"
            "Please re-create/re-upload the file and run dimensions update again."
        ) from None

    if "failed to create index" in stderr_lower:
        raise TaskUserError(
            f"bcftools failed while creating an index for {target}.\nDetails from bcftools:\n{stderr_clean}"
        ) from None

    raise TaskUserError(
        f"bcftools failed while {operation} for {target}.\nDetails from bcftools:\n{stderr_clean}"
    ) from None


def ensure_csi_index(file: Path) -> None:
    """
    Ensure a .csi index exists for the given VCF file, creating it if absent.
    Used by both BcftoolsQueryManager (query execution) and VCFDimensionCalculator (dimension indexing).

    This function also is the main guard against unsorted VCF in DivBase, since 'bcftools index' requires
    sorted filed and will fail with a specific error if not.
    """
    index_file = file.with_suffix(file.suffix + ".csi")
    if not index_file.exists():
        index_command = f"index -f {file}"
        proc = run_bcftools(command=index_command, capture_stderr=True)
        _, stderr = proc.communicate()

        if proc.returncode != 0:
            stderr = (stderr or "").strip()
            _raise_task_user_error_from_bcftools_stderr(
                stderr=stderr,
                operation="creating a CSI index",
                target=f"VCF file '{file}'",
            )

            error_details = stderr or f"Process exited with code {proc.returncode}"
            raise BcftoolsCommandError(command=index_command, error_details=error_details) from None


def bgzip_vcf_for_indexing(input_vcf: Path, output_vcf_gz: Path) -> None:
    """
    Bgzip a plain .vcf file to .vcf.gz using bcftools view.
    Shared helper for workflows that need CSI indexing, since bcftools index --csi requires bgzipped input.

    The whitelist for file extensions in packages/divbase-lib/src/divbase_lib/divbase_constants.py
    does currently not allow .vcf but required (bgzipped) .vcf.gz
    As such, this function is a guard against future changes to the whilelist,
    since bcftools index requires bgzip-compressed VCFs (https://samtools.github.io/bcftools/bcftools.html#index)
    """
    bgzip_command = f"view -Oz -o {shlex.quote(str(output_vcf_gz))} {shlex.quote(str(input_vcf))}"
    proc = run_bcftools(command=bgzip_command, capture_stderr=True)
    _, stderr = proc.communicate()

    if proc.returncode != 0:
        stderr = (stderr or "").strip()
        _raise_task_user_error_from_bcftools_stderr(
            stderr=stderr,
            operation="bgzipping the VCF for indexing",
            target=f"VCF file '{input_vcf}'",
        )
        error_details = stderr or f"Process exited with code {proc.returncode}"
        raise BcftoolsCommandError(command=bgzip_command, error_details=error_details) from None

    logger.info(f"Bgzipped {input_vcf} to {output_vcf_gz} for CSI indexing.")
