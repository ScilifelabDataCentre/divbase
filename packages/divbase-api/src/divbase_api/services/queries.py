"""
TODO - consider split metadata query and bcftools query into two separate modules.
"""

import contextlib
import datetime
import logging
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import psutil

from divbase_lib.api_schemas.queries import SampleMetadataQueryTaskResult
from divbase_lib.exceptions import (
    BcftoolsCommandError,
    BcftoolsEnvironmentError,
    BcftoolsPipeEmptyCommandError,
    BcftoolsPipeUnsupportedCommandError,
    SidecarColumnNotFoundError,
    SidecarInvalidFilterError,
    SidecarNoDataLoadedError,
)

logger = logging.getLogger(__name__)


@dataclass
class BCFToolsInput:
    """
    Contains the inputs required to run a bcftools query.
    """

    sample_and_filename_subset: list[dict[str, str]]
    sampleIDs: list[str]
    filenames: list[str]


def run_sidecar_metadata_query(
    file: Path,
    filter_string: str = None,
    project_id: int = None,
    vcf_dimensions_data: dict = None,
) -> SampleMetadataQueryTaskResult:
    """
    Run a query on a sidecar metadata TSV file and map samples to VCF files.

    Takes vcf_dimensions_data fetched by an API call in the task layer

    """

    sidecar_manager = SidecarQueryManager(file=file).run_query(filter_string=filter_string)
    query_message = sidecar_manager.query_message
    warnings = sidecar_manager.warnings
    unique_sample_ids = sidecar_manager.get_unique_values("Sample_ID")

    logger.info(f"Metadata query returned {len(unique_sample_ids)} unique sample IDs")

    if not vcf_dimensions_data or not vcf_dimensions_data.get("vcf_files"):
        error_msg = f"No VCF dimensions data provided for project {project_id}. "
        error_msg += "Please run 'divbase-cli dimensions update' first."
        raise ValueError(error_msg)

    sample_and_filename_subset = []
    unique_filenames = set()

    for vcf_entry in vcf_dimensions_data.get("vcf_files", []):
        filename = vcf_entry["vcf_file_s3_key"]
        sample_names = vcf_entry.get("samples", [])

        for sample_id in sample_names:
            if sample_id in unique_sample_ids:
                sample_and_filename_subset.append({"Sample_ID": sample_id, "Filename": filename})
                unique_filenames.add(filename)

    return SampleMetadataQueryTaskResult(
        sample_and_filename_subset=sample_and_filename_subset,
        unique_sample_ids=list(unique_sample_ids),
        unique_filenames=list(unique_filenames),
        query_message=query_message,
        warnings=warnings,
    )


class BcftoolsQueryManager:
    """
    A class that manages the execution of querys that require bcftools.

    # TODO - support different file paths for input files.

    Intended for use with a Celery architechture to run the queries as synchronous or asynchronous jobs.
    The bottom layer of the class - self.run_bcftools() - is designed for either being run inside a Celery worker container
    upon receiving a task from the queue (async job), or to be run inside the same Docker container as the Celery worker with
    `docker exec` instead of the queue (synchronous job). Either way, the class expects that the worker container with the
    name defined in CONTAINER_NAME is running.

    NOTE! The current implementation does not handle starting the container.

    Users are expected to interact with this class through the CLI. The CLI layer handles Celery task management for this class.

    However, it is possible to run a query without Celery using the class directly, e.g.:
    `output_file = BcftoolsQueryManager().execute_pipe(command, bcftools_inputs)`

    The main entry point is the `execute_pipe` method, which takes a semicolon-separated string of
    bcftools commands and a dictionary of inputs. It builds a configuration structure for the commands,
    ensure that the relevant Docker container is available (for synchronous runs), and processes the commands in
    a structured manner.

    Current implementation is a "merge-last" strategy. In short, all subsetting and filtering commands
    are executed on each input files, temp files are created, and the results are merged at the end.
    The benefit of this approach is that merge operations can be memory intensive, and by subsetting first,
    the row-columns of the input files for the eventual merge operation are smaller.
    Downside is that the logic is more complex than in a "merge-first" strategy.

    The "merge-first" logic in based on nesting of two loops:
    - Outer loop: iterates over the input bcftools commands in the pipe, and passes each command to `run_current_command()`.
    - Inner loop: iterates over the input files for each command , and runs the command on each file using `run_bcftools()`.
    The results of the inner loop are temporary files that are passed to the next command in the outer loop.
    At the end of the outer loop, the temporary files are merged into a single output file using `merge_or_concat_bcftools_temp_files()`.
    The temporary files are cleaned up after the processing is done using `cleanup_temp_files()`.
    The class also provides a context manager for temporary file management to ensure that temporary files are cleaned
    up even if the processing fails or exits unexpectedly.

    """

    VALID_BCFTOOLS_COMMANDS = ["view"]  # white-list of valid bcftools commands to run in the pipe.
    CONTAINER_NAME = "divbase-worker-quick-1"  # for synchronous tasks, use this container name to find the container ID
    ENABLE_SUBPROCESS_MONITORING = (
        os.environ.get("ENABLE_WORKER_METRICS_PER_TASK", "true").lower() == "true"
    )  # Monitoring bcftools subprocesses is controled with an environment variable. Comes with some overhead, hence optional.

    def execute_pipe(self, command: str, bcftools_inputs: dict, job_id: int) -> tuple[str, dict[str, float]]:
        """
        Main entrypoint for executing executing divbase queries that require bcftools.
        First calls on a method to build a structure of input parameters for bcftools, and then
        passes that on to another method that process the commands according to the "merge-last" strategy.

        Returns a tuple of (output_file, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses.
        """

        walltime_start = time.time()

        in_docker = os.path.exists("/.dockerenv")
        in_k8s = self._is_in_kubernetes()
        if not in_docker and not in_k8s:
            logger.info("Running outside Docker container, ensuring Docker container is available")
            try:
                get_container_id = self.get_container_id(self.CONTAINER_NAME)
                logger.info(f"Found the required {self.CONTAINER_NAME} container running with ID: {get_container_id}")
            except BcftoolsEnvironmentError:
                raise

        identifier = job_id if job_id else datetime.datetime.now().timestamp()

        commands_config_structure = self.build_commands_config(command, bcftools_inputs, identifier)
        output_file, metrics = self.process_bcftools_commands(commands_config_structure, identifier)

        walltime_end = time.time()
        metrics["walltime_seconds"] = walltime_end - walltime_start
        logger.info(f"Total bcftools pipeline walltime: {metrics['walltime_seconds']:.2f}s")

        return output_file, metrics

    @contextlib.contextmanager
    def temp_file_management(self):
        """Context manager to handle temporary file cleanup, even if processing fails/exits."""
        self.temp_files = []
        try:
            yield self
        finally:
            if self.temp_files:
                logger.info(f"Cleaning up {len(self.temp_files)} temporary files")
                self.cleanup_temp_files(self.temp_files)

    def build_commands_config(
        self, command: str, bcftools_inputs: dict[str, Any], identifier: str = None
    ) -> list[dict[str, Any]]:
        """
        Method that builds a configuration structure for the bcftools commands based on the provided command string and inputs.
        The command string is expected to be a semicolon-separated list of bcftools commands.
        Each command is processed to create a list of dictionaries containing the command details,
        input files, sample subsets, and output temporary files. Returns a list of dictionaries
        where each dictionary represents a command configuration with the bcftools command, input files and temporary
        output files.
        """
        filenames = bcftools_inputs.get("filenames")
        sample_and_filename_subset = bcftools_inputs.get("sample_and_filename_subset")

        if not command or command.strip() == ";" or command.strip() == "":
            raise BcftoolsPipeEmptyCommandError()
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
                raise BcftoolsPipeUnsupportedCommandError(
                    command=cmd_name, position=c_counter + 1, valid_commands=self.VALID_BCFTOOLS_COMMANDS
                )

            output_temp_files = [
                f"temp_subset_{identifier}_{c_counter}_{f_counter}.bcf" for f_counter, _ in enumerate(current_inputs)
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

    def process_bcftools_commands(
        self, commands_config: list[dict[str, Any]], identifier: str
    ) -> tuple[str, dict[str, float]]:
        """
        Method that handles the outer loop of the merge-last strategy: it loops over each of commands in
        the input and passes them to the command runner run_current_command() (which in turn handles the inner loop:
        looping over VCF files). Once the outer loop is done, it calls a method to merge all temporary results files
        into a single output VCF file.

        Returns a tuple of (output_file, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses.
        """
        with self.temp_file_management() as temp_file_manager:
            logger.info(f"Loaded configuration with {len(commands_config)} commands in the pipe")

            final_output_temp_files = None

            # Accumulate metrics across all commands
            total_cpu_seconds = 0.0
            max_peak_memory = 0
            all_memory_samples = []

            for cmd_config in commands_config:
                logger.info(f"Processing command #{cmd_config['counter'] + 1}: {cmd_config['command']}")
                output_temp_files, cmd_metrics = self.run_current_command(cmd_config)
                temp_file_manager.temp_files.extend(output_temp_files)
                final_output_temp_files = output_temp_files

                total_cpu_seconds += cmd_metrics.get("cpu_seconds", 0)
                max_peak_memory = max(max_peak_memory, cmd_metrics.get("peak_memory_bytes", 0))
                if cmd_metrics.get("avg_memory_bytes", 0) > 0:
                    all_memory_samples.append(cmd_metrics["avg_memory_bytes"])

            output_file = self.merge_or_concat_bcftools_temp_files(final_output_temp_files, identifier)

            logger.info("bcftools processing completed successfully")

            avg_memory = sum(all_memory_samples) / len(all_memory_samples) if all_memory_samples else 0

            metrics = {
                "cpu_seconds": total_cpu_seconds,
                "peak_memory_bytes": max_peak_memory,
                "avg_memory_bytes": avg_memory,
            }

            return output_file, metrics

    def run_current_command(self, cmd_config: dict[str, Any]) -> tuple[list[str], dict[str, float]]:
        """
        Method that handles the inner loop of the merge-last strategy: for each command pass from the outer loop,
        it processes all given input VCF files individually by running the command on each file using run_bcftools().
        For each processed file, a temporary output file is created, which is then used as input for the next command
        in the outer loop. Each temporary output file is indexed with a .csi index file using ensure_csi_index().

        Returns a tuple of (output_temp_files, metrics) where metrics contains accumulated CPU and memory stats
        for all bcftools subprocesses executed in this command.
        """
        command = cmd_config["command"]
        input_files = cmd_config["input_files"]
        output_temp_files = cmd_config["output_temp_files"]
        sample_subset = cmd_config["sample_subset"]

        # Accumulate metrics across all bcftools subprocess calls
        total_cpu_seconds = 0.0
        peak_memory = 0
        memory_samples = []

        for f_counter, file in enumerate(input_files):
            temp_file = output_temp_files[f_counter]

            samples_in_file = []
            for record in sample_subset:
                if record["Filename"] == file:
                    samples_in_file.append(record["Sample_ID"])

            samples_in_file_bcftools_formatted = ",".join(samples_in_file)

            cmd_with_samples = command.strip().replace("SAMPLES", samples_in_file_bcftools_formatted)
            formatted_cmd = f"{cmd_with_samples} {file} -Ou -o {temp_file}"

            # Run bcftools and optionally monitor the subprocess
            proc = self.run_bcftools(command=formatted_cmd)

            if self.ENABLE_SUBPROCESS_MONITORING:
                current_cpu = 0.0  # Initialize to track CPU usage
                loop_iterations = 0

                try:
                    # Monitor the bcftools subprocess
                    bcftools_proc = psutil.Process(proc.pid)
                    cpu_start = bcftools_proc.cpu_times()
                    logger.info(
                        f"Started monitoring bcftools PID {proc.pid}, initial CPU: user={cpu_start.user}, system={cpu_start.system}"
                    )

                    # Take initial memory sample immediately
                    try:
                        mem = bcftools_proc.memory_info().rss
                        memory_samples.append(mem)
                        peak_memory = max(peak_memory, mem)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass

                    # Sample memory and CPU while process is running
                    while proc.poll() is None:
                        loop_iterations += 1
                        try:
                            mem = bcftools_proc.memory_info().rss
                            memory_samples.append(mem)
                            peak_memory = max(peak_memory, mem)

                            cpu_current = bcftools_proc.cpu_times()
                            current_cpu = (cpu_current.user + cpu_current.system) - (cpu_start.user + cpu_start.system)

                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            break
                        # Set interval between measurements
                        time.sleep(0.01)
                        # TODO sleep time could be longer to reduce overhead, but then fast processes might be missed

                    logger.info(
                        f"Monitoring loop completed after {loop_iterations} iterations, current_cpu={current_cpu:.6f}s"
                    )

                    # Wait for process to complete and check return code
                    returncode = proc.wait()
                    if returncode != 0:
                        raise BcftoolsCommandError(
                            command=formatted_cmd, error_details=f"Process exited with code {returncode}"
                        ) from None

                    # Add CPU to total (will be 0.0 if process was too fast to measure)
                    total_cpu_seconds += current_cpu
                    logger.info(
                        f"Bcftools subprocess finished: CPU={current_cpu:.6f}s, Peak memory={peak_memory} bytes, Samples={len(memory_samples)}"
                    )

                except psutil.NoSuchProcess:
                    returncode = proc.wait()
                    if returncode != 0:
                        raise BcftoolsCommandError(
                            command=formatted_cmd, error_details=f"Process exited with code {returncode}"
                        ) from None
                    logger.info("Bcftools process finished too quickly to monitor - caught NoSuchProcess")
            else:
                # Monitoring disabled, just wait for process to complete
                returncode = proc.wait()
                if returncode != 0:
                    raise BcftoolsCommandError(
                        command=formatted_cmd, error_details=f"Process exited with code {returncode}"
                    ) from None
                logger.info("Bcftools subprocess finished (monitoring disabled)")

            self.ensure_csi_index(temp_file)
            self._log_file_size(temp_file)

        # Calculate average memory
        avg_memory = sum(memory_samples) / len(memory_samples) if memory_samples else 0

        metrics = {
            "cpu_seconds": total_cpu_seconds,
            "peak_memory_bytes": peak_memory,
            "avg_memory_bytes": avg_memory,
        }

        return output_temp_files, metrics

    def run_bcftools(self, command: str) -> subprocess.Popen:
        """
        Methid to run a bcftools command inside a Docker container that has bcftools installed.
        This method is specifically designed to with a Celery manager in an upper layer of the architechture.
        In short, the CLI layer handles the task management and submission to Celery, which can be either synchronous or asynchronous.

        Synchronous tasks (submitted with .apply()) skips the queue and tries to find the container id of the running container
        and runs the command inside it using `docker exec` and `subprocess.Popen`. I.e. the host machine executes the command
        inside the container.

        Asynchronous tasks (submitted with .apply_async()) are queued, picked up by the celery worker container and then execeuted
        inside the containter with `subprocess.Popen`. I.e. the container itself executes the command inside itself.

        To identify if the job is running async, the method evaluates if the current process is running inside a docker container by checking
        for the existence of the /.dockerenv file. If the file exists, it assumes that the command is run asynchronously inside the Celery worker
        container. If the file does not exist, it assumes that the command is run synchronously and tries to find the container ID of the running
        bcftools container using the get_container_id() method.

        Returns the subprocess.Popen object so the caller can monitor resource usage.
        """
        logger.info(f"Running: bcftools {command}")

        in_docker = os.path.exists("/.dockerenv")
        in_k8s = self._is_in_kubernetes()

        if in_docker or in_k8s:
            logger.debug("Running inside Celery worker Docker container, executing bcftools directly")
            try:
                proc = subprocess.Popen(["bcftools"] + command.split())
                return proc
            except Exception as e:
                logger.error(f"Failed to run bcftools directly: {e}")
                raise BcftoolsCommandError(command=command, error_details=str(e)) from e
        else:
            try:
                container_id = self.get_container_id(self.CONTAINER_NAME)
                logger.debug(f"Executing command in container with ID: {container_id}")
                docker_cmd = ["docker", "exec", container_id, "bcftools"] + command.split()
                proc = subprocess.Popen(docker_cmd)
                return proc
            except BcftoolsEnvironmentError:
                raise
            except Exception as e:
                logger.error(f"Failed to run bcftools in container: {e}")
                raise BcftoolsCommandError(command=command, error_details=str(e)) from e

    def ensure_csi_index(self, file: str) -> None:
        """
        Helper method that ensures that the given VCF file has a .csi index. If not, create it using bcftools.

        (bcftools can sometimes handle VCF files that lack an index file, but for consistency
        it is better to create an index file for all VCF files that are created.)
        """
        index_file = f"{file}.csi"
        if not os.path.exists(index_file):
            index_command = f"index -f {file}"
            proc = self.run_bcftools(command=index_command)
            proc.wait()

    def merge_or_concat_bcftools_temp_files(self, output_temp_files: list[str], identifier: str) -> str:
        """
        Helper method that merges the final temporary files produced by pipe_query_command into a single output file.

        This method adds some temp files list used by temp_file_management() for cleanup after processing. Note that
        the results file is not included in this list. Clean up of the results file is handled by tasks.py: only
        after successfull upload to S3 will the results file be deleted.

        # for all sets that have len(files) > 1, perform concat, save the temp filename to a new list
        # for all sets that have len(files) == 1, save the temp filename to a new list
        # for all the temp filenames in the new list, perform merge


        """
        # TODO handle naming of output file better, e.g. by using a timestamp or a unique identifier

        unsorted_output_file = f"merged_unsorted_{identifier}.bcf"
        annotated_unsorted_output_file = f"merged_annotated_unsorted_{identifier}.bcf"
        divbase_header_for_vcf = "divbase_header.txt"
        self.temp_files.append(unsorted_output_file)
        self.temp_files.append(annotated_unsorted_output_file)
        self.temp_files.append(divbase_header_for_vcf)

        output_file = f"result_of_job_{identifier}.vcf.gz"
        logger.info("Trying to determine if sample names overlap between temp files...")

        sample_names_per_VCF = self._get_all_sample_names_from_vcf_files(output_temp_files)
        sample_set_to_files = self._group_vcfs_by_sample_set(sample_names_per_VCF)
        non_overlapping_sample_names = self._check_non_overlapping_sample_names(sample_set_to_files)

        if len(output_temp_files) > 1:
            if non_overlapping_sample_names:
                logger.info("Sample names do not overlap between temp files, will continue with 'bcftools merge'")
                merge_command = f"merge --force-samples -Ou -o {unsorted_output_file} {' '.join(output_temp_files)}"
                # TODO double check if this should use output_temp_files or if that is an old remnant. the code below uses sample_set_to_files but that is perhaps to decide between concat and merge
                proc = self.run_bcftools(command=merge_command)
                proc.wait()
                logger.info(f"Merged all temporary files into '{unsorted_output_file}'.")
                self._log_file_size(unsorted_output_file)
            else:
                logger.info(
                    "Sample names overlap between some temp files, will concat overlapping sets, then merge if needed and possible."
                )
                temp_concat_files = []
                for sample_set, files in sample_set_to_files.items():
                    logger.debug(f"Processing sample set with {len(files)} files ({files}) and samples: {sample_set}")
                    if len(files) > 1:
                        logger.debug("Sample set occurs in multiple files, will concat these files.")
                        concat_temp = f"concat_{identifier}_{hash(sample_set)}.bcf"
                        concat_command = f"concat -Ou -o {concat_temp} {' '.join(files)}"
                        proc = self.run_bcftools(command=concat_command)
                        proc.wait()
                        temp_concat_files.append(concat_temp)
                        self.temp_files.append(concat_temp)
                        self.ensure_csi_index(concat_temp)
                        self._log_file_size(concat_temp)
                    elif len(files) == 1:
                        logger.debug(
                            "Sample set only occurs in a single file, will use this file as is for merging in a downstream step."
                        )
                        temp_concat_files.append(files[0])
                if len(temp_concat_files) > 1:
                    merge_command = f"merge --force-samples -Ou -o {unsorted_output_file} {' '.join(temp_concat_files)}"
                    proc = self.run_bcftools(command=merge_command)
                    proc.wait()
                    logger.info(f"Merged all files (including concatenated files) into '{unsorted_output_file}'.")
                    self._log_file_size(unsorted_output_file)
                elif len(temp_concat_files) == 1:
                    self._log_file_size(temp_concat_files[0])
                    os.rename(temp_concat_files[0], unsorted_output_file)
                    logger.info(
                        f"Only one file remained after concatenation, renamed this file to '{unsorted_output_file}'."
                    )
                    self._log_file_size(unsorted_output_file)
        elif len(output_temp_files) == 1:
            logger.info(f"Only one file was produced by the query, renamed this file to '{unsorted_output_file}'.")
            os.rename(output_temp_files[0], unsorted_output_file)
            self._log_file_size(unsorted_output_file)

        self._prepare_txt_with_divbase_header_for_vcf(header_filename=divbase_header_for_vcf)
        annotate_command = (
            f"annotate -h {divbase_header_for_vcf} -Ou -o {annotated_unsorted_output_file} {unsorted_output_file}"
        )
        proc = self.run_bcftools(command=annotate_command)
        proc.wait()
        self._log_file_size(annotated_unsorted_output_file)

        sort_command = f"sort -Oz -o {output_file} {annotated_unsorted_output_file}"
        proc = self.run_bcftools(command=sort_command)
        proc.wait()
        self._log_file_size(output_file)
        logger.info(
            f"Sorting the results file to ensure proper order of variants. Final results are in '{output_file}'."
        )

        return output_file

    def cleanup_temp_files(self, output_temp_files: list[str]) -> None:
        """
        Helper method that handles deletion of all temporary files and their associated .csi index files
        that were generated during the pipe_query_command execution.
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

    def get_container_id(self, container_name: str) -> str:
        """
        Helper method to get the container ID of the running bcftools Docker container.
        Raises BcftoolsEnvironmentError if the container is not found or if the command fails.
        The error can the be re-raised be the methods that call this method, e.g. run_bcftools() and execute_pipe().
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

    def _get_all_sample_names_from_vcf_files(self, output_temp_files: list[str]) -> dict[str, list[str]]:
        """
        Helper method that is used to determine if there are any sample names that recur across the temp files.
        If they do, bcftools concat is needed instead of bcftools merge.

        Use bcftools query to extract sample names from BCF/VCF files since it is designed for such operations.
        """
        sample_names_per_VCF = {}
        for vcf_file in output_temp_files:
            try:
                result = subprocess.run(
                    ["bcftools", "query", "-l", vcf_file],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                sample_names = result.stdout.strip().split("\n")
                sample_names_per_VCF[vcf_file] = sample_names if sample_names != [""] else []
            except subprocess.SubprocessError as e:
                logger.error(f"Failed to get sample names from {vcf_file}: {e}")
                raise

        return sample_names_per_VCF

    def _group_vcfs_by_sample_set(self, sample_names_per_VCF: dict[str, list[str]]) -> dict[tuple, list[str]]:
        """
        Helper method that groups VCF files by their sample sets. VCF files that contain the same sample set
        (=completely overlapping samples) need to be combined using bcftools concat instead of bcftools merge.
        Here, frozenset is used to create an immutable set of sample names for each VCF file.
        sample_set_to_files then stores all files that contain the same sample set.
        """
        sample_set_to_files = {}
        for vcf_file, sample_list in sample_names_per_VCF.items():
            sample_set = tuple(sample_list)
            sample_set_to_files.setdefault(sample_set, []).append(vcf_file)
        return sample_set_to_files

    def _check_non_overlapping_sample_names(self, sample_set_to_files: dict[frozenset, list[str]]) -> bool:
        """
        Helper method that looks at a mapping of sample set to VCF files and checks for non-overlapping sample names.
        Simply put, if any sample set in the input dict has more than one file, samples overlap between files
        """
        return not any(len(files) > 1 for files in sample_set_to_files.values())

    def _is_in_kubernetes(self) -> bool:
        """
        Check if k8s environment variable is set in containter.
        """
        return "KUBERNETES_SERVICE_HOST" in os.environ

    def _prepare_txt_with_divbase_header_for_vcf(self, header_filename: str) -> None:
        """
        The command 'bcftools annotate -h' can be used to append a custom line to the header of a VCF file.
        The command requires a text file. This method checks if such a file exists, and if not, creates it.
        The 'annotate' command expects the pattern '##key=value' in the text file.
        """
        try:
            with open(header_filename, "w", encoding="utf-8") as fh:
                timestamp_in_format_used_by_bcftools = datetime.datetime.now().strftime("%a %b %d %H:%M:%S %Y")
                fh.write(
                    f'##DivBase_created="This is a results file created by a DivBase query; Date={timestamp_in_format_used_by_bcftools}"\n'
                )
                logger.info("Added header to VCF saying that this results file was created with DivBase.")
        except Exception as e:
            logger.warning(f"Could not write {header_filename}: {e}")

    def _log_file_size(self, file_path: str):
        """
        Log the size of the a given file in both GB and GiB.
        """

        # TODO consider changing to logger debug later in the dev process
        try:
            size_bytes = os.path.getsize(file_path)
            size_gb = size_bytes / (1000 * 1000 * 1000)
            size_gi = size_bytes / (1024 * 1024 * 1024)
            logger.info(f"File '{file_path}' size: {size_gb:.2f} GB, {size_gi:.2f} Gi")
        except Exception as e:
            logger.warning(f"Could not determine size of file '{file_path}': {e}")


class SidecarQueryManager:
    """
    A class that manages the execution queries on sidecar metadata files.

    Expects a TSV file with a header row and tab-separated values.
    The class provides methods to load the TSV file into a pandas DataFrame, run queries against the data,
    and retrieve unique values from specific columns.

    TODO - consider seperation of concerns: this class currently handles both loading the TSV file and running queries on it.
    TODO - some of the __init__ params are perhaps better as properties?
    """

    def __init__(self, file: Path):
        self.file = file
        self.filter_string = None
        self.df = None
        self.query_result = None
        self.query_message: str = ""
        self.warnings: list[str] = []
        self.load_file()

    def load_file(self) -> "SidecarQueryManager":
        """
        Method that loads the TSV file into a pandas DataFrame. Assumes that the first row is a header row, and that the file is tab-separated.
        Also removes any leading '#' characters from the column names

        Strip empty filenames if there e.g. are typos with trailing commas
        """
        # TODO: pandas will likely read all plain files to df, so perhaps there should be a check that the file is a TSV file? or at least has properly formatted tabular columns and rows?
        try:
            logger.info(f"Loading sidecar metadata file: {self.file}")
            self.df = pd.read_csv(
                self.file, sep="\t"
            )  # Pandas has Type Inference and will detect numeric and string columns automatically
            self.df.columns = self.df.columns.str.lstrip("#")
            if "Sample_ID" not in self.df.columns:
                raise SidecarColumnNotFoundError("The 'Sample_ID' column is required in the metadata file.")

        except Exception as e:
            raise SidecarNoDataLoadedError(file_path=self.file, submethod="load_file") from e
        return self

    def get_unique_values(self, column: str) -> list:
        """
        Method to fetch unique values from a specific column in the query result. Intended to be invoked on a SidecarQueryManager
        instance after a query has been run with run_query().
        """
        if self.query_result is None:
            raise SidecarColumnNotFoundError("No query result available. Run run_query() first.")

        if column in self.query_result.columns:
            return self.query_result[column].unique().tolist()
        else:
            raise SidecarColumnNotFoundError(f"Column '{column}' not found in query result")

    def run_query(self, filter_string: str = None) -> "SidecarQueryManager":
        """
        Method to run a query against the loaded data. The filter_string should be a semicolon-separated list of key:value pairs,
        where key is a column name and value is a comma-separated list of filter values.
        For example: "key1:value1,value2;key2:value3,value4".

        The TSV that is loaded into the pandas DataFrame can have both string and numeric columns.
        - String columns are matched to filter string values with OR logic: if ANY value in a cell matches ANY filter value, the row matches.
        - Numeric columns support:
            - Inequalities: ">25", "<=40" (checks if any cell value satisfies the condition)
            - Ranges: "20-40" (checks if any cell value is within the range)
            - Discrete values: "25,30,50" (checks if any cell value matches any filter value)
            - All are combined with OR logic

        Filtering using the ! (NOT) operator:
        - "!" must prefix the filter value, e.g. "key:!value" means that rows with "value" in the "key" column should be excluded.
        - Numeric examples: "Population:!2" (exclude 2), "Age:<30,!25" (less than 30 but not 25), "Weight:!20-40" (exclude range 20-40)
        - String examples: "Area:!North" (exclude North), "Region:East,West,!South" (East or West but not South)
        - NOT conditions are applied with AND logic after positive conditions have been applied: rows must NOT match any negated value

        Filter string values in the query vs. cell values in the TSV:
        - Filter strings are handled per semicolon-separated key-value pair: in "key1:value1,value2;key2:value3,value4"
          "key1:value1,value2" is handled separately from "key2:value3,value4".
        - Filter string values can be comma-separated, e.g. "value1,value2" in "key1:value1,value2" and each filter string value is handled separately.
        - Cell values can be semicolon-separated, e.g. "25;30;35" in a TSV cell
        - Matching of filter string to cell values uses OR logic: if ANY value in a cell matches ANY filter value, the row matches.
          E.g. "key2:value3,value4" means that TSV cells in the "key2" column that contain "value3" will match, but also cells that contain "value3;value4" or "value4;value3" will match.

        Summary of how different input filter values are handled:
        - If the filter_string is empty, all records are returned.
        - If the filter_string is None, an error is raised.
        - If the filter_string is not empty, the method filters the dataframe based on the provided filter_string.
        - If any of the keys in the filter_string are not found in the dataframe columns, a warning is logged and those conditions are skipped.
        - If none of the filter string values match any cell values in the dataframe, a warning is logged and all records are returned.
        - If the filter_string is invalid, a SidecarInvalidFilterError is raised.

        The method returns the SidecarQueryManager instance with the query_result and query_message. The former is the filtered DataFrame results,
        and the latter is filter_string used for the query.
        """

        if self.df is None:
            raise SidecarNoDataLoadedError(
                file_path=self.file, submethod="run_query", error_details="No data loaded. Call load_file() first."
            )

        if filter_string is not None:
            self.filter_string = filter_string

        if self.filter_string == "":
            logger.warning("Empty filter provided - returning ALL records. This may be a large result set.")
            self.query_result = self.df
            self.query_message = "ALL RECORDS (no filter)"
            return self

        if self.filter_string is None:
            raise SidecarInvalidFilterError("Filter cannot be None. Use an empty string ('') if you want all records.")

        key_values = self.filter_string.split(";")
        filter_conditions = []

        # 1. Parse the input filter string and build a list of boolean conditions to apply to the dataframe
        for key_value in key_values:
            if not key_value.strip():
                continue
            try:
                key, filter_string_values = key_value.split(":", 1)
                key = key.strip()
                filter_string_values = filter_string_values.strip()

                if key not in self.df.columns:
                    warning_msg = f"Column '{key}' not found in the TSV file. Skipping this filter condition."
                    logger.warning(warning_msg)
                    self.warnings.append(warning_msg)
                    continue

                is_numeric = pd.api.types.is_numeric_dtype(self.df[key])
                is_semicolon_numeric = self._is_semicolon_separated_numeric_column(key) if not is_numeric else False

                # Handle numeric filtering: inequalities, ranges, and discrete values (all with OR logic)
                # e.g., "Weight:>25,<30,50" or "Weight:20-40,50,>100"
                # Supports filtering on semicolon-separated values in cells in the TSV: e.g. "25;30;35"
                # Also handles columns that pandas infers as strings but contain numeric values with semicolons (e.g., "1;2;3")
                # Also supports NOT operator with ! prefix: e.g., "Weight:!25" or "Weight:<4,!2"
                if is_numeric or is_semicolon_numeric:
                    filter_string_values_list = filter_string_values.split(",")

                    # Negated values are those that start with "!" in the filter string
                    positive_values, negated_values = self._separate_positive_and_negated_values(
                        filter_values=filter_string_values_list
                    )

                    filter_context = {
                        "key": key,
                        "filter_string_values": filter_string_values,
                        "is_negated": False,
                    }

                    inequality_conditions, range_conditions, discrete_values = self._parse_numeric_filter_values(
                        values_to_process=positive_values,
                        context=filter_context,
                    )

                    filter_context["is_negated"] = True
                    negated_inequality_conditions, negated_range_conditions, negated_discrete_values = (
                        self._parse_numeric_filter_values(
                            values_to_process=negated_values,
                            context=filter_context,
                        )
                    )

                    # Combine multiple conditions (inequality, range, discrete values) with OR logic
                    conditions = self._build_condition_list(
                        inequality_conditions=inequality_conditions,
                        range_conditions=range_conditions,
                        discrete_values=discrete_values,
                        key=key,
                    )

                    negated_conditions = self._build_condition_list(
                        inequality_conditions=negated_inequality_conditions,
                        range_conditions=negated_range_conditions,
                        discrete_values=negated_discrete_values,
                        key=key,
                    )

                    if conditions or negated_conditions:
                        # First combine all positive conditions with OR logic. Can be None if there are no positive conditions, e.g. if the filter string only contains negated conditions like "Weight:!20-40"
                        base_condition = self._combine_conditions_with_or(conditions=conditions) if conditions else None
                        # Then apply negated conditions with AND logic: rows must NOT match any negated condition.
                        combined = self._apply_not_conditions(
                            base_condition=base_condition, negated_conditions=negated_conditions
                        )

                        if not combined.any():
                            warning_msg = f"No values in column '{key}' match the filter: {filter_string_values}"
                            logger.warning(warning_msg)
                            self.warnings.append(warning_msg)
                        filter_conditions.append(combined)
                        logger.info("filter_conditions: " + str(filter_conditions))  # debug
                    else:
                        warning_msg = f"No valid numeric values, ranges, or inequalities provided for column '{key}'. Filter condition will not match any rows."
                        logger.warning(warning_msg)
                        self.warnings.append(warning_msg)
                else:
                    # Non-numeric column: handle as discrete string values
                    # Supports NOT operator with ! prefix: e.g., "Area:!North" or "Area:North,!South"
                    filter_string_values_list = filter_string_values.split(",")
                    self._validate_no_commas_in_column(key)

                    positive_values, negated_values = self._separate_positive_and_negated_values(
                        filter_values=filter_string_values_list
                    )

                    # Build condition
                    if positive_values or negated_values:
                        base_condition = (
                            self._create_string_condition(key=key, target_values=positive_values)
                            if positive_values
                            else None
                        )
                        negated_conditions = (
                            [self._create_string_condition(key=key, target_values=negated_values)]
                            if negated_values
                            else []
                        )

                        condition = self._apply_not_conditions(
                            base_condition=base_condition, negated_conditions=negated_conditions
                        )

                        if not condition.any():
                            warning_msg = f"None of the values {filter_string_values_list} were found in column '{key}'"
                            logger.warning(warning_msg)
                            self.warnings.append(warning_msg)
                        filter_conditions.append(condition)
            except SidecarInvalidFilterError:
                # Allow specific validation errors (like "contains commas") to propagate unchanged.
                # This preserves detailed error messages for user-facing exceptions.
                raise
            except Exception as e:
                raise SidecarInvalidFilterError(
                    f"Invalid filter format: '{key_value}'. Expected format 'key:value1,value2' or 'key:min-max' for numeric ranges"
                ) from e

        # 2. Apply the final boolean filters on the dataframe
        if filter_conditions:
            combined_condition = pd.Series(True, index=self.df.index)
            # Iteratively combine each condition in filter_conditions to create a final combined condition where each row must satisfy all filter conditions to be included.
            for condition in filter_conditions:
                # The ampersand (&) is pandas syntax for element-wise AND between boolean Series.
                combined_condition = combined_condition & condition

            self.query_result = self.df[combined_condition].copy()
            self.query_message = self.filter_string
        else:
            warning_msg = "Invalid filter conditions: none of the filters matched any records. Returning ALL records. This may be a large result set. Please check your filter keys, value spelling, and syntax."
            logger.warning(warning_msg)
            self.warnings.append(warning_msg)
            self.query_result = self.df
            self.query_message = f"Invalid filter conditions ({self.filter_string}) - returning ALL records"

        return self

    def _validate_no_commas_in_column(self, key: str) -> None:
        """
        Helper method to validate that column values in the imported TSV does not contain commas.
        Raises SidecarInvalidFilterError if any comma is found in the column values.
        """
        for row_index, cell_value in enumerate(self.df[key].dropna()):
            cell_str = str(cell_value).strip()
            if cell_str and "," in cell_str:
                raise SidecarInvalidFilterError(
                    f"Column '{key}' contains commas in value '{cell_str}' at row {row_index}. "
                    f"Commas are not allowed in DivBase metadata files. Use semicolons (;) to separate multiple values."
                )

    def _is_semicolon_separated_numeric_column(self, key: str) -> bool:
        """
        Helper method to detect if a column contains semicolon-separated numeric values.
        Pandas correctly infers type from single-value columns. But for columns with semicolon-separated values
        (e.g.: "1;2;3"), it infers them as strings (object dtype). This is an issue since numeric operations
        (inequalities, ranges) cannot be performed on string values.

        This helper method checks ALL non-null values in the column to determine if they can be parsed as numeric
        after splitting by semicolon. Note that it only detects if a column value is semicolon-separated numeric, it does not
        convert the column to numeric type. It also validates that the colum value is not of mixed string-numerical type, which
        is invalid input for the query system. The actual parsing and handling of the semicolon-separated numeric values is
        done in the lamda functions in the run_query() method.

        Returns True if the column contains ONLY numeric values (with or without semicolons).
        Returns False if the column contains ONLY non-numeric values (regular string column), or if the column is empty.
        Raises SidecarInvalidFilterError if mixed types detected (e.g., "1;2" and "abc;def" in same column) or if commas are found.
        """
        if key not in self.df.columns:
            return False

        non_null_values = self.df[key].dropna()
        if len(non_null_values) == 0:
            return False

        self._validate_no_commas_in_column(key)

        has_numeric_type = False
        has_non_numeric_type = False

        for row_index, cell_value in enumerate(non_null_values):
            cell_str = str(cell_value).strip()
            if not cell_str:
                continue

            parts = cell_str.split(";")
            for part in parts:
                part = part.strip()
                if not part:
                    continue

                # Check if the value contains a hyphen and looks like it could be numeric (e.g., "1-2", "3-4")
                if "-" in part and any(p.isdigit() for p in part):
                    raise SidecarInvalidFilterError(
                        f"Column '{key}' contains value '{part}' with a hyphen at row {row_index}. "
                        f"Hyphens are not allowed in numeric column values (only in string columns). "
                        f"If this is meant to be a string column, all values should be non-numeric strings."
                    )

                try:
                    float(part)
                    has_numeric_type = True
                except ValueError:
                    has_non_numeric_type = True

                if has_numeric_type and has_non_numeric_type:
                    raise SidecarInvalidFilterError(
                        f"Column '{key}' in the metadata file contains mixed types. Value '{cell_str}' at row {row_index} "
                        f"has both numeric and non-numeric parts. All values in a column must be consistently "
                        f"numeric or string for filtering to work correctly."
                    )

        return has_numeric_type and not has_non_numeric_type

    def _split_cell_values(self, cell_value: Any) -> list[str]:
        """
        Helper method to split cell value by semicolon and return list of non-empty values.
        If the cell contains a single value without semicolon, it will return a list with that single value.
        If the cell is empty or NaN, it will return an empty list.
        """
        if pd.isna(cell_value):
            return []
        return [val.strip() for val in str(cell_value).split(";") if val.strip()]

    def _parse_numeric_value(self, value_str: str) -> float | int:
        """Helper method to parse a string value to int or float. To be used when other checks have already confirmed that the value can be parsed as numeric."""
        return float(value_str) if "." in value_str else int(value_str)

    def _create_inequality_condition(self, key: str, operator: str, threshold: float) -> pd.Series:
        """
        Helper method to create a condition for inequality filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to defined the
        logic that will be applied to the Pandas dataframe.
        """

        def check_inequality(cell_value):
            if pd.isna(cell_value):
                return False
            cell_values = self._split_cell_values(cell_value)
            for val_str in cell_values:
                try:
                    val_num = self._parse_numeric_value(val_str)
                    if (
                        operator == ">"
                        and val_num > threshold
                        or operator == ">="
                        and val_num >= threshold
                        or operator == "<"
                        and val_num < threshold
                        or operator == "<="
                        and val_num <= threshold
                    ):
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_inequality)

    def _create_range_condition(self, key: str, min_val: float, max_val: float) -> pd.Series:
        """
        Helper method to create a condition for range filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_range(cell_value):
            if pd.isna(cell_value):
                return False
            cell_values = self._split_cell_values(cell_value)
            for val_str in cell_values:
                try:
                    val_num = self._parse_numeric_value(val_str)
                    if min_val <= val_num <= max_val:
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_range)

    def _create_discrete_numeric_condition(self, key: str, target_values: list[float | int]) -> pd.Series:
        """
        Helper method to create a condition for discrete numeric value filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_discrete(cell_value):
            if pd.isna(cell_value):
                return False
            cell_values = self._split_cell_values(cell_value)
            for val_str in cell_values:
                try:
                    val_num = self._parse_numeric_value(val_str)
                    if val_num in target_values:
                        return True
                except ValueError:
                    continue
            return False

        return self.df[key].apply(check_discrete)

    def _create_string_condition(self, key: str, target_values: list[str]) -> pd.Series:
        """
        Helper method to create a condition for string value filtering on a column.
        Uses a named nested function instead of a lambda to improve readability to define the
        logic that will be applied to the Pandas dataframe.
        """

        def check_string(cell_value):
            if pd.isna(cell_value):
                return False
            cell_values = self._split_cell_values(cell_value)
            return any(val in target_values for val in cell_values)

        return self.df[key].apply(check_string)

    def _combine_conditions_with_or(self, conditions: list[pd.Series]) -> pd.Series:
        """
        Helper method to combine multiple Pandas boolean Series with OR logic.
        Returns a single boolean Series that is True if any of the input conditions is True for each row.

        The resulting Series is used at the end of self.run_query() to filter the DataFrame values.
        """
        if not conditions:
            return pd.Series(False, index=self.df.index)
        combined = conditions[0]
        for cond in conditions[1:]:
            # The bar (|) is pandas syntax for element-wise OR between boolean Series.
            combined = combined | cond
        return combined

    def _build_condition_list(
        self,
        inequality_conditions: list[pd.Series],
        range_conditions: list[pd.Series],
        discrete_values: list[float | int],
        key: str,
    ) -> list[pd.Series]:
        """
        Helper method to build a list of conditions from inequality, range, and discrete value filters.
        """
        conditions = []

        if inequality_conditions:
            conditions.append(self._combine_conditions_with_or(conditions=inequality_conditions))

        if range_conditions:
            conditions.append(self._combine_conditions_with_or(conditions=range_conditions))

        if discrete_values:
            discrete_condition = self._create_discrete_numeric_condition(key=key, target_values=discrete_values)
            conditions.append(discrete_condition)

        return conditions

    def _parse_numeric_filter_values(
        self, values_to_process: list[str], context: dict[str, str | bool]
    ) -> tuple[list[pd.Series], list[pd.Series], list[float | int]]:
        """
        Helper method to identify if a numeric filter values is an inequality, range, or discrete value and process it accordingly.

        The context dict is intended to keep the kwargs manageable when passing positive and negative values back-to-back:
            - key: Column name being filtered
            - filter_string_values: Original filter string (for error messages)
            - is_negated: Whether these are negated (NOT) conditions
        """
        key = context["key"]
        filter_string_values = context["filter_string_values"]
        is_negated = context["is_negated"]

        inequality_conditions = []
        range_conditions = []
        discrete_values = []

        for filter_string_value in values_to_process:
            # Check for common mistakes: =< or => instead of <= or >=
            if re.match(r"^=<\d+\.?\d*$", filter_string_value) or re.match(r"^=>\d+\.?\d*$", filter_string_value):
                raise SidecarInvalidFilterError(
                    f"Invalid operator format '{filter_string_value[:2]}' in filter '{key}:{filter_string_values}'."
                    f"Use standard operators: '<=' (not '=<') or '>=' (not '=>')"
                )

            # Check if it's an inequality (e.g., ">25", "<=40")
            inequality_match = re.match(r"^(>=|<=|>|<)(\d+\.?\d*)$", filter_string_value)
            if inequality_match:
                operator = inequality_match.group(1)
                threshold = float(inequality_match.group(2))
                condition = self._create_inequality_condition(key, operator, threshold)
                inequality_conditions.append(condition)
                prefix = "NOT " if is_negated else ""
                logger.debug(
                    f"Applied {'negated ' if is_negated else ''}inequality filter on '{key}': {prefix}{operator} {threshold}"
                )
                continue

            # Check if it's a range (e.g., "20-40")
            range_match = re.match(r"^(\d+\.?\d*)-(\d+\.?\d*)$", filter_string_value)
            if range_match:
                min_val = float(range_match.group(1))
                max_val = float(range_match.group(2))
                condition = self._create_range_condition(key, min_val, max_val)
                range_conditions.append(condition)
                prefix = "NOT " if is_negated else ""
                logger.debug(
                    f"Applied {'negated ' if is_negated else ''}range filter on '{key}': {prefix}{min_val} to {max_val}"
                )
                continue

            # Otherwise, treat as discrete value
            try:
                numeric_value = float(filter_string_value) if "." in filter_string_value else int(filter_string_value)
                discrete_values.append(numeric_value)
            except ValueError:
                logger.warning(
                    f"Cannot convert '{filter_string_value}' to numeric for column '{key}'. Skipping this value."
                )

        return inequality_conditions, range_conditions, discrete_values

    def _separate_positive_and_negated_values(self, filter_values: list[str]) -> tuple[list[str], list[str]]:
        """
        Helper method to separate filter values into positive and negated lists.
        Values prefixed with '!' are negated (NOT conditions).
        """
        positive_values = []
        negated_values = []

        for value in filter_values:
            value = value.strip()
            if value.startswith("!"):
                negated_values.append(value[1:].strip())
            else:
                positive_values.append(value)

        return positive_values, negated_values

    def _apply_not_conditions(self, base_condition: pd.Series | None, negated_conditions: list[pd.Series]) -> pd.Series:
        """
        Helper method to apply NOT conditions to a base condition. The base condition contains positive filters combined with OR, or None if there were only negations
        in the input filter string from the CLI. Returns a combined condition where rows must match base_condition AND NOT match any negated condition
        """
        if base_condition is None:
            # If only negated conditions (no positive conditions), start with all True
            combined = pd.Series(True, index=self.df.index)
        else:
            combined = base_condition

        # Apply negated conditions (must NOT match any negated condition)
        for negated_condition in negated_conditions:
            combined = combined & ~negated_condition

        return combined
