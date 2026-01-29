"""
Helper script to fetch per-task metrics from Prometheus for specified job IDs and dates (configurable in an external YAML file).

Assumes that a prometheus server is running locally on port 9090 and has scraped the relevant metrics.
To scrape per-task metrics, the environment variables ENABLE_WORKER_METRICS and ENABLE_WORKER_METRICS_PER_TASK need
to have been set to True in the worker at the time the job ran.

To be able to use this script for a variety of historical jobs in the Prometheus data store, a YAML file is used to specify
the job IDs and the dates when the jobs ran. The date is used to specify the time range for the Prometheus query.

The 9090/api/v1/query_range endpoint is used to get historical data in the Prometheus data store, and it requires a start and end time.
Therefore the YAML needs to contain the job IDs and the dates when the jobs ran.

YAML example:
jobs:
  - job_id: 23
    date: 2026-01-28
  - job_id: 28
    date: 2026-01-29


Usage:

python scripts/benchmarking/fetch_per_task_metrics_from_prometheus.py --yaml scripts/benchmarking/task_metrics.yaml

"""

import argparse
import json
import subprocess

import yaml

PROM_URL = "http://localhost:9090/api/v1/query_range"
METRICS = [
    "celery_task_cpu_seconds_total",
    "celery_task_python_overhead_cpu_seconds",
    "celery_task_memory_peak_bytes",
    "celery_task_memory_avg_bytes",
    "celery_task_bcftools_cpu_seconds_total",
    "celery_task_bcftools_memory_peak_bytes",
    "celery_task_bcftools_memory_avg_bytes",
    "celery_task_bcftools_walltime_seconds",
    "celery_task_vcf_download_walltime_seconds",
    "celery_task_vcf_download_cpu_seconds",
    "celery_task_vcf_download_memory_peak_bytes",
    "celery_task_vcf_download_memory_avg_bytes",
    "celery_task_walltime_seconds",
]


def load_jobs(yaml_path):
    """
    Load job IDs and dates from an external YAML file.
    """
    with open(yaml_path) as f:
        return yaml.safe_load(f)["jobs"]


def fetch_metric(metric, job_id, start, end):
    """
    Curl the Prometheus API to fetch the specified metric for the given job ID and time range.
    In order to get historical data in the Prometheus data store, the 9090/api/v1/query_range
    endpoint is used, and it requires a start and end time.
    """
    url = f"{PROM_URL}?query={metric}%7Bjob_id%3D%22{job_id}%22%7D&start={start}&end={end}&step=60"
    result = subprocess.run(["curl", "-s", url], capture_output=True, text=True)
    return result.stdout


def extract_single_value(data):
    """
    The worker metrics endpoint only exposes per-task metrics after a task is completed.
    However, the metrics are kept in a TTL cache for some time so Prometheus might store
    duplicates in the data store. Since the reported metric values will be static for a
    given task, it is enough to extract the last value from the time series to deduplicate '
    the results.
    """
    try:
        parsed = json.loads(data)
        results = parsed.get("data", {}).get("result", [])
        if results and "values" in results[0] and results[0]["values"]:
            # Take the last value (timestamp, value)
            _, value = results[0]["values"][-1]
            return value
        else:
            return "No data"
    except Exception as e:
        return f"Error: {e}"


def main():
    parser = argparse.ArgumentParser(description="Fetch per-task metrics from Prometheus for specified job IDs.")
    parser.add_argument("--yaml", "-y", required=True, help="Path to the YAML file with job IDs and dates.")
    args = parser.parse_args()

    jobs = load_jobs(args.yaml)
    for job in jobs:
        job_id = job["job_id"]
        date = job["date"]
        start = f"{date}T06:00:00Z"
        end = f"{date}T18:00:00Z"
        for metric in METRICS:
            data = fetch_metric(metric, job_id, start, end)
            value = extract_single_value(data)
            print(f"Job {job_id}, Metric {metric}: {value}")


if __name__ == "__main__":
    main()
