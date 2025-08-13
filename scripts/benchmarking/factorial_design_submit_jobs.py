"""
A script to submit jobs for a factorial design experiment to benchmark how the dimensions of a VCF file affect processing time.
This script assumes that there is a docker image that contains 'fake-vcf' (https://github.com/endast/fake-vcf).

The script generates a factorial design and then uses 'fake-vcf' to generate the VCF files. Mock sample metadata is generated for each VCF file.
These are then uploaded to a local minIO bucket and a bcftools-queries are submitted to the job system to subset the VCF files based on the metadata.

Currently, the script has code for two different types of experimental designs: full factorial and Latin hypercube.
Latin hypercube is random, but if the same number of samples AND same random seed is used, the results can be reproduced.

Since tasks can take long to complete for larger VCF files, this script does not poll for task completion. Instead it writes the task IDs to a json file.
The json file can later be read by the companion script 'factorial_design_analyze_results.py', which fetches runtimes from the completed tasks and analyses the results.

Usage:
    python scripts/benchmarking/factorial_design_submit_jobs.py
"""

import itertools
import json
import os
import re
import shlex
import subprocess
from pathlib import Path

import numpy as np
from _benchmarking_shared_utils import (
    LOCAL_ENV,
    ensure_project_exists,
    ensure_project_in_config,
)
from scipy.stats import qmc


def generate_full_factorial_design(factor1_levels: np.ndarray, factor2_levels: np.ndarray) -> list[tuple]:
    """Generate all combinations of factor1_levels and factor2_levels (full factorial design)."""
    return list(itertools.product(factor1_levels.astype(int), factor2_levels.astype(int)))


def generate_latin_hypercube_design(
    n_samples: int, factor1_levels: np.ndarray, factor2_levels: np.ndarray, random_seed: int = None
) -> list[tuple]:
    """
    Run a latin hypercube.
    Note! The nature of the latin hypercube is that it will produce different results when repeated. They should still cover the same general design space.
    """
    sampler = qmc.LatinHypercube(d=2, seed=random_seed)
    samples = sampler.random(n=n_samples)

    indices = (samples * (len(factor1_levels) - 1)).round().astype(int)
    factor1 = factor1_levels[indices[:, 0]]
    factor2 = factor2_levels[indices[:, 1]]

    hypercube_design = []
    for i in range(n_samples):
        hypercube_design.append((int(factor1[i]), int(factor2[i])))

    return hypercube_design


def ensure_required_files_in_bucket(project_name: str, vcf_filenames: list) -> None:
    """Check if the required files are in the bucket. Otherwise: download the vcf file from the URL from the public repository;
    generate the mock metadata; upload the files to the bucket so that the worker container can access them from there during the query."""

    print("\nChecking if files are already in bucket...")
    cmd = f"divbase-cli files list --project {project_name}"
    env = dict(os.environ)
    env["DIVBASE_ENV"] = "local"

    result = subprocess.run(shlex.split(cmd), env=env, capture_output=True, text=True, check=True)
    output = result.stdout

    for vcf_filename in vcf_filenames:
        if vcf_filename in output:
            print(f"{vcf_filename} is already present in the bucket.")
        else:
            print(f"{vcf_filename} is NOT present in the bucket. Checking to see if it exists locally...")
            if not os.path.isfile(vcf_filename):
                print(f"{vcf_filename} not found locally. Please generate it first.")
                continue
            print(f"{vcf_filename} found locally. Uploading to bucket...")
            cmd_upload_vcf = shlex.split(f"divbase-cli files upload --project {project_name} {vcf_filename}")
            subprocess.run(cmd_upload_vcf, check=True, env=LOCAL_ENV)

        metadata_filename = vcf_filename.replace("mock_vcf_file", "mock_metadata_file").replace(".vcf.gz", ".tsv")
        if metadata_filename in output:
            print(f"{metadata_filename} is already present in the bucket.")
        else:
            print(f"{metadata_filename} is NOT present in the bucket. Checking to see if it exists locally...")
            if not os.path.isfile(metadata_filename):
                print(f"{metadata_filename} not found locally. Generating...")
                cmd_generate_metadata = (
                    f"python scripts/generate_mock_sample_metadata.py --vcf {vcf_filename} --output {metadata_filename}"
                )
                subprocess.run(shlex.split(cmd_generate_metadata), check=True)
                print("Mock metadata generated.")
            print(f"Uploading {metadata_filename} to bucket...")
            cmd_upload_metadata = shlex.split(f"divbase-cli files upload --project {project_name} {metadata_filename}")
            subprocess.run(cmd_upload_metadata, check=True, env=LOCAL_ENV)

    # TODO use batch upload to the bucket!


def generate_mock_vcfs(design: list[tuple], random_seed: int) -> list[str]:
    """Generate mock VCF files based on the factorial design and random seed. Uses fake-vcf to do so, which is installed in a separate docker image"""

    # TODO check if docker image exists, otherwise build it.
    # TODO keep docker image spun up until all files are generated
    # TODO first check if file in bucket, then run docker command to generate it

    output_files = []
    for _, (s, r) in enumerate(design):
        s_int = int(round(s))
        r_int = int(round(r))
        output_file = f"mock_vcf_file_s_{s_int}_r_{r_int}.vcf.gz"
        output_files.append(output_file)
        output_path = Path(output_file)
        if not output_path.exists():
            docker_cmd = (
                f"docker run --rm -v {Path.cwd()}:/data benchmarking "
                f"fake-vcf generate -s {s_int} -r {r_int} -o /data/{output_file} --seed {random_seed}"
            )
            subprocess.run(shlex.split(docker_cmd), check=True)
        else:
            print(f"{output_file} already exists for random seed {random_seed}, skipping.")
    return output_files


def submit_jobs_to_queue(output_files: list[str], project_name: str, n_replicates: int = 1) -> None:
    task_records = []
    for replicate in range(n_replicates):
        for filename in output_files:
            metadata_filename = filename.replace("mock_vcf_file", "mock_metadata_file").replace(".vcf.gz", ".tsv")
            cmd_query = f"divbase-cli query bcftools-pipe --tsv-filter 'Area:North,East' --command 'view -s SAMPLES' --metadata-tsv-name {metadata_filename} --project {project_name}"
            result = subprocess.run(shlex.split(cmd_query), check=True, env=LOCAL_ENV, capture_output=True, text=True)
            output = result.stdout

            match = re.search(r"task id:\s*([a-f0-9\-]+)", output, re.IGNORECASE)
            task_id = match.group(1) if match else None

            sr_match = re.search(r"_s_(\d+)_r_(\d+)", filename)
            if sr_match:
                samples = int(sr_match.group(1))
                variants = int(sr_match.group(2))
            else:
                samples = None
                variants = None
                print("Pattern not found")

            task_records.append(
                {
                    "filename": filename,
                    "metadata_filename": metadata_filename,
                    "cmd_query": cmd_query,
                    "task_id": task_id,
                    "number_of_samples": samples,
                    "number_of_variants": variants,
                    "replicate": replicate + 1,
                }
            )
    return task_records

    # TODO should be able to send single or nested bcftools view query


def main():
    n_samples = 30  # only used for  latin hypercube experiments, the number of samples should be at least the same as the number of levels in each factor
    factor1_levels = np.linspace(10, 1000, num=10)
    factor2_levels = np.logspace(np.log10(10), np.log10(1_000_000), num=10)
    random_seed: int = 12345
    n_replicates = 3
    project_name = "factorial-design"

    design = generate_full_factorial_design(factor1_levels, factor2_levels)

    # TODO allow running latin hyoercube instead with argparse
    # design = generate_latin_hypercube_design(
    #     n_samples=n_samples, factor1_levels=factor1_levels, factor2_levels=factor2_levels, random_seed=random_seed
    # )

    output_files = generate_mock_vcfs(design, random_seed)

    ensure_project_exists(project_name)

    ensure_project_in_config(project_name)

    ensure_required_files_in_bucket(
        project_name=project_name,
        vcf_filenames=output_files,
    )

    task_records = submit_jobs_to_queue(output_files, project_name, n_replicates)

    with open("task_records.json", "w") as f:
        json.dump(task_records, f, indent=2)

    # TODO add an option to specify the name of the task_records file

    print("The tasks have been submitted to the job queue. Wait until all tasks are finished, and then run:")
    print("python scripts/benchmarking/factorial_design_analyze_results.py")
    print("to analyze the results.")


if __name__ == "__main__":
    main()
