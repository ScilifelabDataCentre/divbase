import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def tsv_query_command(file: Path, filter: str) -> tuple[pd.DataFrame, str]:
    """
    Query a TSV file for specific key-value pairs. Assumes that the TSV file has a header row
    and a column named Filename. The user can use any column header as a key, and the values
    can be any string. Returns the unique filenames that match the query.
    """
    # TODO if key not in df.columns, there is currently no warning or error
    # TODO if value not in df[key], there is no warning or error. only if both values are missing
    # is there is message saying that an empty df is returned
    # TODO what if the user wants to make queries on the Filename column? It should work, but might result in wierd edge-cases?

    df = pd.read_csv(file, sep="\t")
    df.columns = df.columns.str.lstrip("#")

    key_values = filter.split(";")
    filter_conditions = []
    for key_value in key_values:
        key, values = key_value.split(":", 1)
        values_list = values.split(",")

        if key in df.columns:
            condition = f"{key} in {values_list}"
            filter_conditions.append(condition)

    query_string = " and ".join(filter_conditions)
    query_result = df.query(query_string)

    return query_result, query_string


def pipe_query_command(command: str, bcftools_inputs: dict) -> None:
    """
    Ensure that the bcftools Docker image is available, then pass "query" commands to bcftools.

    Current implementation is a "merge-last" strategy. In short, all subsetting and filtering commands
    are executed on each input files, temp files are created, and the results are merged at the end.
    The benefit of this approach is that merge operations can be memory intensive, and by subsetting first,
    the row-columns of the input files for the eventual merge operation are smaller.
    Downside is that the logic is more complex than in a "merge-first" strategy.
    """
    IMAGE_NAME = "bcftools-image"

    image_exists, container_id = check_bcftools_docker_image(image_name=IMAGE_NAME)

    if not image_exists:
        logger.info(f"Docker image '{IMAGE_NAME}' not found. Building it now...")
        build_bcftools_docker_image(image_name=IMAGE_NAME)

    filenames = bcftools_inputs.get("filenames")
    sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")
    all_temp_files = []
    command_list = command.split(";")
    commands_config_structure = []
    current_inputs = filenames
    serializable_samples = sample_and_filename_subset.to_dict(orient="records")

    for c_counter, cmd in enumerate(command_list):
        command_details = {
            "command": cmd,
            "counter": c_counter,
            "input_files": current_inputs,
            "sample_subset": serializable_samples,
        }

        temp_files = [f"temp_subset_{c_counter}_{f_counter}.vcf.gz" for f_counter, _ in enumerate(current_inputs)]

        command_details["temp_files"] = temp_files
        commands_config_structure.append(command_details)

        all_temp_files.extend(temp_files)
        current_inputs = temp_files

    execute_bcftools_job_in_container(commands_config_structure=commands_config_structure, container_id=container_id)


def execute_bcftools_job_in_container(commands_config_structure: list[dict], container_id: str) -> str:
    with tempfile.NamedTemporaryFile(prefix="bcftools_config_", suffix=".json", mode="w", delete=False) as f:
        json.dump(commands_config_structure, f, indent=2)
        temp_config_file = f.name
    try:
        container_config_file = "/app/bcftools_divbase_job_config.json"
        subprocess.run(["docker", "cp", temp_config_file, f"{container_id}:{container_config_file}"], check=True)

        logger.info("Executing bcftools job with commands structure in container...")
        subprocess.run(
            [
                "docker",
                "exec",
                "-w",
                "/app",
                container_id,
                "python",
                "/app/src/divbase_tools/bcftools_runner_for_container.py",
                "--config",
                container_config_file,
            ],
            check=True,
        )

        logger.info("Job completed successfully.")
    finally:
        os.remove(temp_config_file)


def check_bcftools_docker_image(image_name: str) -> bool:
    """
    Check if the bcftools Docker image is available locally.
    The docker comand returns the image ID if the image exists,
    or an empty string if it does not.
    """

    def get_container_id():
        """Helper function to get the container ID of the running bcftools Docker container."""
        result = subprocess.run(
            ["docker", "ps", "--filter", "name=docker-tools", "--format", "{{.ID}}"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    result = subprocess.run(
        ["docker", "images", "-q", image_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True
    )
    image_exists = bool(result.stdout.strip())

    container_id = get_container_id()

    if not container_id:
        logger.warning("Starting bcftools-image...")
        subprocess.run(["docker", "compose", "-f", "docker/docker-compose.yaml", "up", "-d"], check=True)
        container_id = get_container_id()

    logger.info(f"Using existing docker-compose container: {container_id}")

    return image_exists, container_id


def build_bcftools_docker_image(image_name: str) -> None:
    dockerfile_path = "./docker/tools.dockerfile"

    subprocess.run(["docker", "build", "-f", dockerfile_path, "-t", image_name, "."], check=True)


def dummy_pipe_query_command() -> None:
    """
    Dummy function that copies over a precompiled results file instead of generating it with bcftools.

    it mocks this command
    python -m divbase_tools query bcftools-pipe --tsv-filter "Area:West of Ireland,Northern Portugal;Sex:F"
    --comand "[MERGE, SUBSET ON SAMPLES, FILTER on -r 21:15000000-25000000]"
    (the last flag is not implemented yet)
    """
    cmd = ["cp", "tests/fixtures/subset.vcf.gz", "./result.vcf.gz"]
    subprocess.run(cmd, check=True)
    logger.info("Dry-run results are found in 'result.vcf.gz'.")
