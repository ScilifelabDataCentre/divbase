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

def pipe_query_command(command: str, bcftools_inputs: dict) -> None:
    """
    Ensure that the bcftools Docker image is available, then pass "query" commands to bcftools.
    """
    IMAGE_NAME = "bcftools-image"

    sampleIDs = bcftools_inputs.get("sampleIDs")
    filename = bcftools_inputs.get("filenames")

    if sampleIDs:
        sampleIDs_bcftools_formatted = ",".join(sampleIDs)
        #TODO: Hard-coded input file for now, should be changed to a user input
        command = f"view -s {sampleIDs_bcftools_formatted} HOM_20ind_17SNPs.vcf.gz -Oz -o subset_samples.vcf.gz"

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

def run_bcftools_docker(command: str) -> None:
    """
    Run a bcftools command in a Docker container.
    """

    #TODO /app is used for now since this is where bcftools is installed, might want to change?

    command_args = command.split()
    try:
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{Path.cwd()}:/app",
            "bcftools-image",
            "bcftools"
        ] + command_args
        logger.info(f"Using Docker image to run the command: bcftools {command}")
        subprocess.run(cmd, check=True)
        logger.info(f"the bcftools operation completed successfully.")
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