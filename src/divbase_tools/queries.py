from pathlib import Path
import pandas as pd
import subprocess
import os
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

    Current implementation is a "merge-last" strategy. In short, all subsetting and filtering commands
    are executed on each input files, temp files are created, and the results are merged at the end.
    The benefit of this approach is that merge operations can be memory intensive, and by subsetting first, 
    the row-columns of the input files for the eventual merge operation are smaller. 
    Downside is that the logic is more complex than in a "merge-first" strategy.
    """
    IMAGE_NAME = "bcftools-image"

    if not check_bcftools_docker_image(image_name=IMAGE_NAME):
        logger.info(f"Docker image '{IMAGE_NAME}' not found. Building it now...")
        build_bcftools_docker_image(image_name=IMAGE_NAME)

    filenames = bcftools_inputs.get("filenames")
    sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")
    all_temp_files = []
    command_list = command.split(";")

    current_inputs = filenames
    for c_counter, cmd in enumerate(command_list):
        temp_files = process_bcftools_command(cmd, current_inputs, c_counter, sample_and_filename_subset)
        all_temp_files.extend(temp_files)
        current_inputs = temp_files

    merge_bcftools_temp_files(temp_files)
    delete_temp_files(all_temp_files)


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
        logger.info(f"the bcftools operation completed successfully.\n")
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


def ensure_csi_index(file: str) -> None:
    """
    Ensure that the given VCF file has a .csi index. If not, create it using bcftools.

    bcftools can often work on VCF files that lack an index file, but for consistency
    it is better to create an index file for all VCF files.
    """
    index_file = f"{file}.csi"
    if not os.path.exists(index_file):
        logger.info(f"Index file {index_file} not found. Creating index...")
        index_command = f"index -f {file}"
        run_bcftools_docker(command=index_command)


def process_bcftools_command(cmd: str, current_inputs: list, c_counter: int, sample_and_filename_subset: pd.DataFrame) -> list:
    """
    Helper function for pipe_query_command. For each file in current_inputs,
    Process a single command for all input files and return the list of temporary files generated.
    """
    temp_files = []
    for f_counter, file in enumerate(current_inputs):
        samples_in_file = sample_and_filename_subset[sample_and_filename_subset["Filename"] == file]["Sample_ID"].tolist()
        samples_in_file_bcftools_formatted = ",".join(samples_in_file)
        temp_file = f"temp_subset_{c_counter}_{f_counter}.vcf.gz"
        temp_files.append(temp_file)

        cmd_with_samples = cmd.strip().replace("SAMPLES", samples_in_file_bcftools_formatted)
        formatted_cmd = f"{cmd_with_samples} {file} -Oz -o {temp_file}"
        run_bcftools_docker(command=formatted_cmd)
        ensure_csi_index(temp_file)
    return temp_files


def merge_bcftools_temp_files(temp_files: list) -> None:
    """
    Merge all temporary files produced by pipe_query_command into a single output file.
    """
    if len(temp_files) > 1:
        merge_command = f"merge --force-samples -Oz -o merged.vcf.gz {' '.join(temp_files)}"
        run_bcftools_docker(command=merge_command)
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