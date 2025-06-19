import logging
import os
import subprocess
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

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
    IMAGE_NAME = "docker/worker"

    def execute_pipe(self, command: str, bcftools_inputs: dict, run_local_docker: bool = False) -> None:
        """
        Main entrypoint for executing executing divbase queries that require bcftools.
        """
        commands_config_structure = self.build_commands_config(command, bcftools_inputs)

        if run_local_docker:
            self.ensure_docker_image()

        self.process_bcftools_commands(commands_config_structure)

    def build_commands_config(self, command: str, bcftools_inputs: Dict[str, Any]) -> List[Dict[str, Any]]:
        filenames = bcftools_inputs.get("filenames")
        sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")
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

    def process_bcftools_commands(self, commands_config: List[Dict[str, Any]]) -> None:
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

        self.merge_bcftools_temp_files(final_output_temp_files)

        self.cleanup_temp_files(all_output_temp_files)

        logger.info("bcftools processing completed successfully")

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
        logger.info(f"Running: bcftools {command}")
        subprocess.run(["bcftools"] + command.split(), check=True)

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

    def merge_bcftools_temp_files(self, output_temp_files: List[str]) -> None:
        """
        Merge all temporary files produced by pipe_query_command into a single output file.
        """
        # TODO error handling for when temp files are missing or has not been cleaned up since last run (e.g. if the last run aborted)
        if len(output_temp_files) > 1:
            merge_command = f"merge --force-samples -Oz -o merged.vcf.gz {' '.join(output_temp_files)}"
            self.run_bcftools(command=merge_command)
            logger.info("Merged all temporary files into 'merged.vcf.gz'.")

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


def tsv_query_command(file: Path, filter: str = None) -> tuple[pd.DataFrame, str]:
    """
    Query a TSV file for specific key-value pairs. Assumes that the TSV file has a header row
    and a column named Filename. The user can use any column header as a key, and the values
    can be any string. Returns the unique filenames that match the query.
    """
    # TODO if value not in df[key], there is no warning or error. only if both values are missing
    # is there is message saying that an empty df is returned
    # TODO what if the user wants to make queries on the Filename column? It should work, but might result in wierd edge-cases?

    df = pd.read_csv(file, sep="\t")

    if any(col.startswith("#") for col in df.columns):
        df.columns = df.columns.str.lstrip("#")

    if filter == "":
        logger.warning("Empty filter provided - returning ALL records. This may be a large result set.")
        return df, "ALL RECORDS (no filter)"

    if filter is None:
        raise ValueError("Filter cannot be None. Use an empty string ('') if you want all records.")

    key_values = filter.split(";")
    filter_conditions = []

    for key_value in key_values:
        if not key_value.strip():
            continue
        try:
            key, values = key_value.split(":", 1)
            values_list = values.split(",")

            if key in df.columns:
                condition = df[key].isin(values_list)
                if not condition.any():
                    logger.warning(f"None of the values {values_list} were found in column '{key}'")
                filter_conditions.append(condition)
            else:
                logger.warning(f"Column '{key}' not found in the TSV file. Skipping this filter condition.")
        except ValueError:
            logger.warning(f"Invalid filter format: '{key_value}'. Expected format 'key:value1,value2'")

    if filter_conditions:
        combined_condition = pd.Series(True, index=df.index)
        for condition in filter_conditions:
            combined_condition = combined_condition & condition

        filtered_df = df[combined_condition].copy()
        return filtered_df, filter

    logger.warning("Invalid filter conditions found - returning ALL records. This may be a large result set.")
    return df, "Invalid filter conditions - returning ALL records"
