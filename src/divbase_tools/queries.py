import logging
import os
import subprocess
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

    for c_counter, cmd in enumerate(command_list):
        command_details = {
            "command": cmd,
            "counter": c_counter,
            "input_files": current_inputs,
            "sample_subset": sample_and_filename_subset,
        }

        temp_files = [f"temp_subset_{c_counter}_{f_counter}.vcf.gz" for f_counter, _ in enumerate(current_inputs)]

        command_details["temp_files"] = temp_files
        commands_config_structure.append(command_details)

        all_temp_files.extend(temp_files)
        current_inputs = temp_files

    for cmd_details in commands_config_structure:
        logger.info(f"Executing command: {cmd_details['command']}")
        temp_files = process_bcftools_command(
            cmd=cmd_details["command"],
            current_inputs=cmd_details["input_files"],
            c_counter=cmd_details["counter"],
            sample_and_filename_subset=cmd_details["sample_subset"],
            temp_files=cmd_details["temp_files"],
            container_id=container_id,
        )

    merge_bcftools_temp_files(temp_files=temp_files, container_id=container_id)
    delete_temp_files(all_temp_files)


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


def run_bcftools_docker(command: str, container_id: str) -> None:
    """
    Run a bcftools command in a Docker container.
    """

    command_args = command.split()
    try:
        cmd = ["docker", "exec", "-w", "/app", container_id, "bcftools"] + command_args
        logger.info(f"Using existing container {container_id} to run: bcftools {command}")
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"bcftools command failed with return code {e.returncode}")
        raise


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


def ensure_csi_index(file: str, container_id: str) -> None:
    """
    Ensure that the given VCF file has a .csi index. If not, create it using bcftools.

    bcftools can often work on VCF files that lack an index file, but for consistency
    it is better to create an index file for all VCF files.
    """
    index_file = f"{file}.csi"
    if not os.path.exists(index_file):
        index_command = f"index -f {file}"
        run_bcftools_docker(command=index_command, container_id=container_id)


def process_bcftools_command(
    cmd: str,
    current_inputs: list,
    c_counter: int,
    sample_and_filename_subset: pd.DataFrame,
    temp_files: list,
    container_id: str,
) -> list:
    """
    Helper function for pipe_query_command. For each file in current_inputs,
    Process a single command for all input files and return the list of temporary files generated.
    """
    for f_counter, file in enumerate(current_inputs):
        temp_file = temp_files[f_counter]
        samples_in_file = sample_and_filename_subset[sample_and_filename_subset["Filename"] == file][
            "Sample_ID"
        ].tolist()
        samples_in_file_bcftools_formatted = ",".join(samples_in_file)

        cmd_with_samples = cmd.strip().replace("SAMPLES", samples_in_file_bcftools_formatted)
        formatted_cmd = f"{cmd_with_samples} {file} -Oz -o {temp_file}"
        run_bcftools_docker(command=formatted_cmd, container_id=container_id)
        ensure_csi_index(temp_file, container_id)
    return temp_files


def merge_bcftools_temp_files(temp_files: list, container_id: str) -> None:
    """
    Merge all temporary files produced by pipe_query_command into a single output file.
    """
    if len(temp_files) > 1:
        merge_command = f"merge --force-samples -Oz -o merged.vcf.gz {' '.join(temp_files)}"
        run_bcftools_docker(command=merge_command, container_id=container_id)
        logger.info("Merged all temporary files into 'merged.vcf.gz'.")


def delete_temp_files(temp_files: list) -> None:
    """
    Delete all temporary files and their associated .csi index files
    generated during the pipe_query_command execution.
    """
    for temp_file in temp_files:
        try:
            if os.path.exists(temp_file):
                os.remove(temp_file)
            index_file = f"{temp_file}.csi"
            if os.path.exists(index_file):
                os.remove(index_file)
        except Exception as e:
            logger.error(f"Failed to delete temporary file or index {temp_file}: {e}")
