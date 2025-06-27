import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from divbase_tools.services import download_files_command
from divbase_tools.utils import resolve_bucket_name

logger = logging.getLogger(__name__)


class BcftoolsQueryManager:
    """
    A class that manages the execution of querys that require bcftools.

    The main entry point is the `execute_pipe` method, which takes a semicolon-separated string of
    bcftools commands and a dictionary of inputs. It builds a configuration structure for the commands,
    ensure that the relevant Docker image is available (for local runs), and processes the commands in
    a structured manner.

    Current implementation is a "merge-last" strategy. In short, all subsetting and filtering commands
    are executed on each input files, temp files are created, and the results are merged at the end.
    The benefit of this approach is that merge operations can be memory intensive, and by subsetting first,
    the row-columns of the input files for the eventual merge operation are smaller.
    Downside is that the logic is more complex than in a "merge-first" strategy.

    """

    VALID_BCFTOOLS_COMMANDS = ["view"]
    IMAGE_NAME = "docker-worker-1"
    CONTAINER_NAME = "docker-worker-1"

    def execute_pipe(self, command: str, bcftools_inputs: dict, run_local_docker: bool = False) -> str:
        """
        Main entrypoint for executing executing divbase queries that require bcftools.
        Calls on the sub-methods to build the command structure, ensure the Docker image is available,
        and process the commands.
        """
        commands_config_structure = self.build_commands_config(command, bcftools_inputs)

        if run_local_docker:
            self.ensure_docker_image()

        output_file = self.process_bcftools_commands(commands_config_structure)
        return output_file

    def build_commands_config(self, command: str, bcftools_inputs: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Build a configuration structure for the bcftools commands based on the provided command string and inputs.
        The command string is expected to be a semicolon-separated list of bcftools commands.
        Each command is processed to create a list of dictionaries containing the command details,
        input files, sample subsets, and output temporary files.
        """
        filenames = bcftools_inputs.get("filenames")
        sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")

        if not command or command.strip() == ";" or command.strip() == "":
            raise ValueError("Empty command provided. Please specify at least one valid bcftools command.")
        command_list = command.split(";")
        commands_config_structure = []
        current_inputs = filenames

        for c_counter, cmd in enumerate(command_list):
            cmd = cmd.strip()

            if not cmd:
                logger.warning(f"Skipping empty command at position {c_counter + 1} in command pipeline")
                continue

            cmd_name = cmd.split()[0] if cmd and " " in cmd else cmd
            if cmd_name not in self.VALID_BCFTOOLS_COMMANDS:
                error_msg = (
                    f"Unsupported bcftools command '{cmd_name}' at position {c_counter + 1}. "
                    f"Only the following commands are supported for DivBase queries: {', '.join(self.VALID_BCFTOOLS_COMMANDS)}"
                )
                logger.error(error_msg)
                raise ValueError(error_msg)

            output_temp_files = [
                f"temp_subset_{c_counter}_{f_counter}.vcf.gz" for f_counter, _ in enumerate(current_inputs)
            ]

            command_details = {
                "command": cmd,
                "counter": c_counter,
                "input_files": current_inputs,
                "sample_subset": sample_and_filename_subset,
                "output_temp_files": output_temp_files,
            }

            commands_config_structure.append(command_details)

            current_inputs = output_temp_files

        if not commands_config_structure:
            logger.error("No valid commands provided in input string")

        return commands_config_structure

    def process_bcftools_commands(self, commands_config: List[Dict[str, Any]]) -> str:
        """
        Process a JSON-like list of bcftools command configurations.
        Interprets the JSON configuration, processes and runs each command, ensures that files are indexed, and merges the results.
        """

        logger.info(f"Loaded configuration with {len(commands_config)} commands in the pipe")

        all_output_temp_files = []
        final_output_temp_files = None

        for cmd_config in commands_config:
            logger.info(f"Processing command #{cmd_config['counter'] + 1}: {cmd_config['command']}")
            output_temp_files = self.run_current_command(cmd_config)
            final_output_temp_files = output_temp_files
            all_output_temp_files.extend(output_temp_files)

        output_file = self.merge_bcftools_temp_files(final_output_temp_files)

        self.cleanup_temp_files(all_output_temp_files)

        logger.info("bcftools processing completed successfully")

        return output_file

    def run_current_command(self, cmd_config: Dict[str, Any]) -> List[str]:
        """Process a single command configuration."""
        command = cmd_config["command"]
        input_files = cmd_config["input_files"]
        output_temp_files = cmd_config["output_temp_files"]
        sample_subset = cmd_config["sample_subset"]

        for f_counter, file in enumerate(input_files):
            temp_file = output_temp_files[f_counter]

            samples_in_file = []
            for record in sample_subset:
                if record["Filename"] == file:
                    samples_in_file.append(record["Sample_ID"])

            samples_in_file_bcftools_formatted = ",".join(samples_in_file)

            cmd_with_samples = command.strip().replace("SAMPLES", samples_in_file_bcftools_formatted)
            formatted_cmd = f"{cmd_with_samples} {file} -Oz -o {temp_file}"
            self.run_bcftools(command=formatted_cmd)
            self.ensure_csi_index(temp_file)
        return output_temp_files

    def run_bcftools(self, command: str) -> None:
        """
        Run a bcftools command inside the Docker container.
        The method hancles celery sync and async tasks differently:

        sync tasks (submitted with .apply()) skips the queue and tries to find the container id of the running container
        and runs the command inside it using `docker exec` and `subprocess.run`. I.e. the host machine executes the command
        inside the container.

        async tasks (submitted with .apply_async()) are queued, picked up by the celery worker container and then execeuted
        inside the containter with `subprocess.run`. I.e. the container itself executes the command inside itself.

        To identify if the job is running async, check if it is process is running inside a docker container by checking
        for the existence of the /.dockerenv file.
        """
        logger.info(f"Running: bcftools {command}")

        in_docker = os.path.exists("/.dockerenv")

        if in_docker:
            logger.info("Running inside Docker container, executing bcftools directly")
            try:
                subprocess.run(["bcftools"] + command.split(), check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run bcftools directly: {e}")
                raise
        else:
            try:
                container_id = self.get_container_id(self.IMAGE_NAME)
                if container_id:
                    logger.debug(f"Excuting command in container with id{container_id}")
                    docker_cmd = ["docker", "exec", container_id, "bcftools"] + command.split()
                    subprocess.run(docker_cmd, check=True)
                else:
                    logger.warning("No Docker container found, trying to run bcftools locally")
                    subprocess.run(["bcftools"] + command.split(), check=True)
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to run bcftools: {e}")
                raise

    def ensure_csi_index(self, file: str) -> None:
        """
        Ensure that the given VCF file has a .csi index. If not, create it using bcftools.

        bcftools can often work on VCF files that lack an index file, but for consistency
        it is better to create an index file for all VCF files.
        """
        index_file = f"{file}.csi"
        if not os.path.exists(index_file):
            index_command = f"index -f {file}"
            self.run_bcftools(command=index_command)

    def merge_bcftools_temp_files(self, output_temp_files: List[str]) -> str:
        """
        Merge all temporary files produced by pipe_query_command into a single output file.
        """
        # TODO error handling for when temp files are missing or has not been cleaned up since last run (e.g. if the last run aborted)
        # TODO handle naming of output file better, e.g. by using a timestamp or a unique identifier

        output_file = "merged.vcf.gz"

        if len(output_temp_files) > 1:
            merge_command = f"merge --force-samples -Oz -o {output_file} {' '.join(output_temp_files)}"
            self.run_bcftools(command=merge_command)
            logger.info(f"Merged all temporary files into '{output_file}'.")

        return output_file

    def cleanup_temp_files(self, output_temp_files: List[str]) -> None:
        """
        Delete all temporary files and their associated .csi index files
        generated during the pipe_query_command execution.
        """
        for temp_file in output_temp_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
                index_file = f"{temp_file}.csi"
                if os.path.exists(index_file):
                    os.remove(index_file)
            except Exception as e:
                logger.error(f"Failed to delete temporary file or index {temp_file}: {e}")

    def ensure_docker_image(self) -> None:
        """
        Ensure that the image containing bcftools is available; if not, build it.
        """
        image_name = self.IMAGE_NAME
        image_exists = self.check_bcftools_docker_image(image_name=image_name)

        if not image_exists:
            logger.info(f"Docker image '{image_name}' not found. Building it now...")
            self.build_bcftools_docker_image(image_name=image_name)

    def check_bcftools_docker_image(self, image_name: str) -> bool:
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
            check=True,
        )
        image_exists = bool(result.stdout.strip())
        container_id = self.get_container_id(image_name)

        if not container_id:
            logger.warning(f"Starting {image_name}...")
            subprocess.run(["docker", "compose", "-f", "docker/docker-compose.yaml", "up", "-d"], check=True)
            container_id = self.get_container_id(image_name)

        logger.info(f"Using existing docker-compose container: {container_id}")

        return image_exists

    def get_container_id(self, image_name: str) -> str:
        """
        Helper function to get the container ID of the running bcftools Docker container.
        """
        ##debug to list all running containers
        # list_result = subprocess.run(
        #     ["docker", "ps", "--format", "{{.Names}}"], capture_output=True, text=True, check=False
        # )
        # logger.warning(f"Listing all running containers:{list_result}")

        result = subprocess.run(
            ["docker", "ps", "--filter", f"name={image_name}", "--format", "{{.ID}}"],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    def build_bcftools_docker_image(self, image_name: str) -> None:
        """
        Build the Docker image for bcftools if it does not exist.
        """
        dockerfile_path = "./docker/worker.dockerfile"
        subprocess.run(["docker", "build", "-f", dockerfile_path, "-t", image_name, "."], check=True)


class SidecarQueryManager:
    """Class to manage queries on sidecar metadata files."""

    def __init__(self, file: Path):
        self.file = file
        self.filter_string = None
        self.df = None
        self.query_result = None
        self.query_message = ""
        self.load_file()

    def load_file(self) -> "SidecarQueryManager":
        """
        Load the TSV file into a pandas DataFrame.
        Assumes that the first row is a header row and that the file is tab-separated.
        Also removes any leading '#' characters from the column names
        """
        self.df = pd.read_csv(self.file, sep="\t")
        self.df.columns = self.df.columns.str.lstrip("#")
        return self

    def run_query(self, filter_string: str = None) -> "SidecarQueryManager":
        """
        Run a query against the loaded data
        """

        if self.df is None:
            raise ValueError("No data loaded. Call load_file() first")

        if filter_string is not None:
            self.filter_string = filter_string

        if self.filter_string == "":
            logger.warning("Empty filter provided - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "ALL RECORDS (no filter)"
            return self

        if self.filter_string is None:
            raise ValueError("Filter cannot be None. Use an empty string ('') if you want all records.")

        key_values = self.filter_string.split(";")
        filter_conditions = []

        for key_value in key_values:
            if not key_value.strip():
                continue
            try:
                key, values = key_value.split(":", 1)
                values_list = values.split(",")

                if key in self.df.columns:
                    condition = self.df[key].isin(values_list)
                    if not condition.any():
                        logger.warning(f"None of the values {values_list} were found in column '{key}'")
                    filter_conditions.append(condition)
                else:
                    logger.warning(f"Column '{key}' not found in the TSV file. Skipping this filter condition.")
            except ValueError:
                logger.warning(f"Invalid filter format: '{key_value}'. Expected format 'key:value1,value2'")

        if filter_conditions:
            combined_condition = pd.Series(True, index=self.df.index)
            for condition in filter_conditions:
                combined_condition = combined_condition & condition

            self.query_result = self.df[combined_condition].copy()
            self.query_message = self.filter_string
        else:
            logger.warning("Invalid filter conditions found - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "Invalid filter conditions - returning ALL records"

        return self

    def get_unique_values(self, column: str) -> list:
        """Get unique values from a column in the query result"""
        if self.query_result is None:
            raise ValueError("No query result available. Run run_query() first.")

        if column in self.query_result.columns:
            return self.query_result[column].unique().tolist()
        else:
            raise ValueError(f"Column '{column}' not found in query result")


def fetch_query_files_from_bucket(
    bucket_name: str | None, config_path: Path, files: list[str], download_dir: Path = None, bucket_version=None
) -> None:
    """
    Helper function to fetch files needed for queries from the bucket if they do not exist locally.
    """
    if not download_dir:
        download_dir = Path.cwd()

    bucket_name = resolve_bucket_name(bucket_name, config_path)
    download_files_command(
        bucket_name=bucket_name,
        all_files=files,
        download_dir=download_dir,
        bucket_version=bucket_version,
    )
