from pathlib import Path
import pandas as pd
import subprocess

import logging
logger = logging.getLogger(__name__)

def tsv_query_command(file: Path, filter: str) -> tuple[pd.DataFrame, str]:
    """
    Query a TSV file for specific key-value pairs. Assumes that the TSV file has a header row 
    and a column named Filename. The user can use any column header as a key, and the values 
    can be any string. Returns the unique filenames that match the query.
    """

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

    #TODO if key not in df.columns, there is currently no warning or error
    #TODO if value not in df[key], there is no warning or error. only if both values are missing 
    # is there is message saying that an empty df is returned
    # TODO what if the user wants to make queries on the Filename column? It should work, but might result in wierd edge-cases?

def pipe_query_command(command: str) -> None:
    """
    Ensure that the bcftools Docker image is available, then pass "query" commands to bcftools.
    """
    IMAGE_NAME = "bcftools-image"

    if not check_bcftools_docker_image(image_name=IMAGE_NAME):
        logger.info(f"Docker image '{IMAGE_NAME}' not found. Building it now...")
        build_bcftools_docker_image(image_name=IMAGE_NAME)

    run_bcftools_docker(command=command)


def check_bcftools_docker_image(image_name: str) -> bool:
    """
    Check if the bcftools Docker image is available locally.
    The docker comand returns the image ID if the image exists, 
    or an empty string if it does not.
    """
    result = subprocess.run(
        ["docker", "images", "-q", image_name],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    return bool(result.stdout.strip())

def build_bcftools_docker_image(image_name: str) -> None:
    dockerfile_path = "./docker/tools.dockerfile"

    subprocess.run(
        ["docker", "build", "-f", dockerfile_path, "-t", image_name, "."],
        check=True
    )

def run_bcftools_docker(command: str, input_files: str = None, output_file: str = None) -> None:
    """
    Run a bcftools command in a Docker container.
    """
    
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{Path.cwd()}:/data",
        "bcftools-image",
        "bcftools", command #, input_files, "-o", output_file
    ]
    logger.info(f"Using Docker image to run the command: bcftools {command}")
    subprocess.run(cmd, check=True)
    logger.info(f"bcftools {command} completed successfully.")