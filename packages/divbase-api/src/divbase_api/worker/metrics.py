import logging
import os
import socket
import threading
import time
from collections import defaultdict
from datetime import datetime, timedelta

import psutil
from prometheus_client import Gauge, Info, start_http_server

logger = logging.getLogger(__name__)

WORKER_NAME = socket.gethostname()

ENABLE_WORKER_METRICS = os.environ.get("ENABLE_WORKER_METRICS", "true").lower() == "true"


class MemoryMonitor:
    """Monitor memory usage of a process, tracking peak and average usage."""

    def __init__(self, process: psutil.Process, sample_interval: float = 0.5, baseline_memory: int = 0):
        self.process = process
        self.sample_interval = sample_interval
        self.baseline_memory = baseline_memory
        self.samples = []
        self.peak_memory = 0
        self.monitoring = False
        self._thread = None

    def start(self):
        """Start monitoring memory in a background thread."""
        self.monitoring = True
        self._thread = threading.Thread(target=self._monitor, daemon=True)
        self._thread.start()

    def stop(self):
        """Stop monitoring and return statistics."""
        self.monitoring = False
        if self._thread:
            self._thread.join(timeout=2)

        if not self.samples:
            return {"peak_bytes": 0, "avg_bytes": 0}

        peak_delta = max(0, self.peak_memory - self.baseline_memory)
        avg_value = sum(self.samples) / len(self.samples)
        avg_delta = max(0, avg_value - self.baseline_memory)

        return {"peak_bytes": peak_delta, "avg_bytes": avg_delta}

    def _monitor(self):
        """Background monitoring loop."""
        while self.monitoring:
            try:
                if self.process.is_running():
                    mem = self.process.memory_info().rss
                    self.samples.append(mem)
                    self.peak_memory = max(self.peak_memory, mem)
                else:
                    break
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            time.sleep(self.sample_interval)

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()


## Define Prometheus metrics. Use the same prefix to simplify discovery in Prometheus and Grafana UIs.
# System metrics
worker_cpu_percent = Gauge("celery_worker_prom_client_cpu_percent", "CPU usage percent", ["worker_name", "pid"])
worker_memory_mb = Gauge("celery_worker_prom_client_memory_mb", "Memory usage in MB", ["worker_name", "pid"])
worker_memory_percent = Gauge(
    "celery_worker_prom_client_memory_percent", "Memory usage percent", ["worker_name", "pid"]
)
worker_info = Info("celery_worker_prom_client_info", "Worker information")
worker_cpu_time_total = Gauge(
    "celery_worker_prom_client_cpu_time_total", "Total CPU time (user + system) in seconds", ["worker_name", "pid"]
)
worker_memory_bytes = Gauge("celery_worker_prom_client_memory_bytes", "Memory usage in bytes", ["worker_name", "pid"])
worker_num_threads = Gauge("celery_worker_prom_client_num_threads", "Number of threads", ["worker_name", "pid"])
worker_open_fds = Gauge("celery_worker_prom_client_open_fds", "Number of open file descriptors", ["worker_name", "pid"])

# Per-task metrics - These track resource usage for individual Celery tasks
# Each task is identified by job_id and task_name labels
#
# TASK METRICS: Measure the TOTAL resources for the entire task execution.
# Task CPU = Python worker process CPU + bcftools subprocess CPU (artificially summed in tasks.py)
#   - This is NOT a real measurement, but a constructed sum for visualization purposes
#   - Use python_overhead_cpu for actual measured Python worker CPU
# Task Memory = Peak/avg memory (RSS) of Python worker process ONLY (separate from bcftools subprocess)
#   - This is the worker process memory, NOT including bcftools subprocess memory (separate process).
#   - Memory used by bcftools subprocesses is measured separately and is not included in the Python worker’s memory metrics.
#   - Memory overhead cannot be calculated by subtraction (separate process memory spaces).
#   - Memory metrics for the Python worker and bcftools subprocesses are measured independently. Due to separate process memory spaces, these values should not be summed or directly compared as parts of a whole.
task_cpu_seconds = Gauge(
    "celery_task_cpu_seconds_total",
    "ARTIFICIAL SUM: Python overhead CPU + bcftools subprocess CPU. Not a real measurement. Use for stacked visualizations only.",
    ["job_id", "task_name"],
)
task_python_overhead_cpu_seconds = Gauge(
    "celery_task_python_overhead_cpu_seconds",
    "ACTUAL MEASURED: CPU seconds used by Python worker process (overhead). Does NOT include bcftools subprocess.",
    ["job_id", "task_name"],
)
task_memory_peak_bytes = Gauge(
    "celery_task_memory_peak_bytes",
    "Peak memory (RSS) of Python worker process only. Does NOT include bcftools subprocess memory (separate process). Includes VCF download. Sampled every 0.5s.",
    ["job_id", "task_name"],
)
task_memory_avg_bytes = Gauge(
    "celery_task_memory_avg_bytes",
    "Average memory (RSS) of Python worker process only. Does NOT include bcftools subprocess memory (separate process). Includes VCF download. Sampled every 0.5s.",
    ["job_id", "task_name"],
)

# BCFTOOLS SUBPROCESS METRICS: Measure ONLY the bcftools subprocess portion.
# For tasks with multiple bcftools commands, CPU is summed and memory peak is the max across all subprocesses.
# To calculate overhead: python_overhead_cpu = task_cpu - bcftools_cpu
bcftools_monitoring_config = Info(
    "celery_task_bcftools_monitoring_config",
    "Configuration for bcftools subprocess monitoring. Check 'enabled' field to see if monitoring is active.",
)
task_bcftools_step_only_cpu_seconds = Gauge(
    "celery_task_bcftools_cpu_seconds_total",
    "CPU seconds used ONLY by bcftools subprocesses (sum across all bcftools calls). Sampled every 0.01s. "
    "Value of 0.0 may indicate: (1) process was too fast to measure, or (2) monitoring is disabled (check celery_task_bcftools_monitoring_config). "
    "Subtract from task_cpu to get Python overhead.",
    ["job_id", "task_name"],
)
task_bcftools_memory_peak_bytes = Gauge(
    "celery_task_bcftools_memory_peak_bytes",
    "Peak memory (RSS) of bcftools subprocess(es) in bytes. Max peak observed in any single bcftools subprocess (not summed). Sampled every 0.01s. "
    "Value of 0.0 may indicate: (1) process was too fast to measure, or (2) monitoring is disabled (check celery_task_bcftools_monitoring_config).",
    ["job_id", "task_name"],
)
task_bcftools_memory_avg_bytes = Gauge(
    "celery_task_bcftools_memory_avg_bytes",
    "Average memory (RSS) of bcftools subprocess(es) in bytes. Mean of all samples from all bcftools calls (not an average of averages). Sampled every 0.01s. "
    "Value of 0.0 may indicate: (1) process was too fast to measure, or (2) monitoring is disabled (check celery_task_bcftools_monitoring_config).",
    ["job_id", "task_name"],
)
task_bcftools_walltime_seconds = Gauge(
    "celery_task_bcftools_walltime_seconds",
    "Walltime (elapsed real time) for all bcftools subprocess execution in seconds. Includes all bcftools commands in the pipeline.",
    ["job_id", "task_name"],
)
task_vcf_download_walltime_seconds = Gauge(
    "celery_task_vcf_download_walltime_seconds",
    "Walltime (elapsed real time) for downloading VCF files from S3 in seconds.",
    ["job_id", "task_name"],
)
task_vcf_download_cpu_seconds = Gauge(
    "celery_task_vcf_download_cpu_seconds",
    "CPU seconds used for downloading VCF files from S3 (boto3 operations run in worker process).",
    ["job_id", "task_name"],
)
task_vcf_download_memory_peak_bytes = Gauge(
    "celery_task_vcf_download_memory_peak_bytes",
    "Peak INCREMENTAL memory used during VCF download from S3 (delta from baseline). Measures actual memory increase, not a separate allocation. Do not add this to Python worker memory; it is already included during download. Sampled every 0.5s.",
    ["job_id", "task_name"],
)
task_vcf_download_memory_avg_bytes = Gauge(
    "celery_task_vcf_download_memory_avg_bytes",
    "Average INCREMENTAL memory used during VCF download from S3 (delta from baseline). Measures actual memory increase, not a separate allocation. Do not add this to Python worker memory; it is already included during download. Sampled every 0.5s.",
    ["job_id", "task_name"],
)

task_walltime_seconds = Gauge(
    "celery_task_walltime_seconds",
    "Walltime (elapsed real time) for the entire Celery task execution in seconds.",
    ["job_id", "task_name"],
)

# Metrics cache for per-task metrics to persist across multiple task executions
# Structure: {metric_name: {(job_id, task_name): {"value": float, "timestamp": datetime}}}
task_metrics_cache = defaultdict(dict)
metrics_cache_lock = threading.Lock()  # lock threads to avoid race conditions since multiple threads can access the cache (main thread, celery worker threads, purge thread)

# Prometheus scrapes every 15 seconds in DivBase setup. A TLL of 5 min means it is available for 20 scrapes. Once Prometheus has scraped it, it will store the data in its own volume for its retention time (default 15d).
TASK_METRICS_CACHE_TTL_MINUTES = int(os.environ.get("TASK_METRICS_CACHE_TTL_MINUTES", "5"))


def collect_system_metrics():
    """Collect system metrics for a specific worker process (by PID)."""
    process = psutil.Process()
    pid = str(os.getpid())
    while True:
        try:
            cpu = process.cpu_percent(interval=1)
            worker_cpu_percent.labels(worker_name=WORKER_NAME, pid=pid).set(cpu)
            mem_info = process.memory_info()
            mem_mb = mem_info.rss / 1024 / 1024
            worker_memory_mb.labels(worker_name=WORKER_NAME, pid=pid).set(mem_mb)
            mem_percent = process.memory_percent()
            worker_memory_percent.labels(worker_name=WORKER_NAME, pid=pid).set(mem_percent)

            cpu_times = process.cpu_times()
            worker_cpu_time_total.labels(worker_name=WORKER_NAME, pid=pid).set(cpu_times.user + cpu_times.system)
            worker_memory_bytes.labels(worker_name=WORKER_NAME, pid=pid).set(mem_info.rss)
            worker_num_threads.labels(worker_name=WORKER_NAME, pid=pid).set(process.num_threads())
            # num_fds() is only available on Unix. if open fds increase it might be indicative of a a resource leak
            if hasattr(process, "num_fds"):
                worker_open_fds.labels(worker_name=WORKER_NAME, pid=pid).set(process.num_fds())
        except Exception as e:
            logger.error(f"Error collecting metrics: {e}")
        time.sleep(5)


def store_task_metric(metric_name: str, job_id: int, task_name: str, value: float):
    """Store a task metric in the cache with current timestamp."""
    with metrics_cache_lock:
        task_metrics_cache[metric_name][(job_id, task_name)] = {
            "value": value,
            "timestamp": datetime.now(),
        }


def purge_old_metrics():
    """
    Remove metrics older than TASK_METRICS_CACHE_TTL_MINUTES from cache.
    Each metric is timestamped when stored and can thus be purged individually based on age.
    """
    cutoff_time = datetime.now() - timedelta(minutes=TASK_METRICS_CACHE_TTL_MINUTES)
    with metrics_cache_lock:
        for metric_name, gauge in [
            ("task_cpu_seconds", task_cpu_seconds),
            ("task_python_overhead_cpu_seconds", task_python_overhead_cpu_seconds),
            ("task_bcftools_cpu_seconds", task_bcftools_step_only_cpu_seconds),
        ]:
            tasks_to_remove = [
                key for key, data in task_metrics_cache[metric_name].items() if data["timestamp"] < cutoff_time
            ]
            for key in tasks_to_remove:
                del task_metrics_cache[metric_name][key]
                gauge.remove(*key)
                logger.debug(f"Purged old metric: {metric_name} for job_id={key[0]} from Prometheus client memory")

        for metric_name, gauge in [
            ("task_memory_peak_bytes", task_memory_peak_bytes),
            ("task_memory_avg_bytes", task_memory_avg_bytes),
            ("task_bcftools_memory_peak_bytes", task_bcftools_memory_peak_bytes),
            ("task_bcftools_memory_avg_bytes", task_bcftools_memory_avg_bytes),
            ("task_bcftools_walltime_seconds", task_bcftools_walltime_seconds),
            ("task_vcf_download_walltime_seconds", task_vcf_download_walltime_seconds),
            ("task_vcf_download_cpu_seconds", task_vcf_download_cpu_seconds),
            ("task_vcf_download_memory_peak_bytes", task_vcf_download_memory_peak_bytes),
            ("task_vcf_download_memory_avg_bytes", task_vcf_download_memory_avg_bytes),
            ("task_walltime_seconds", task_walltime_seconds),
        ]:
            tasks_to_remove = [
                key for key, data in task_metrics_cache[metric_name].items() if data["timestamp"] < cutoff_time
            ]
            for key in tasks_to_remove:
                del task_metrics_cache[metric_name][key]
                gauge.remove(*key)
                logger.debug(f"Purged old metric: {metric_name} for job_id={key[0]} from Prometheus client memory")


def get_all_cached_metrics():
    """Get all cached metrics for exposure via Prometheus."""
    with metrics_cache_lock:
        return {
            metric_name: {key: data["value"] for key, data in tasks.items()}
            for metric_name, tasks in task_metrics_cache.items()
        }


def update_prometheus_gauges_from_cache(
    task_cpu_gauge,
    python_overhead_cpu_gauge,
    bcftools_cpu_gauge,
    task_mem_peak_gauge,
    task_mem_avg_gauge,
    bcftools_mem_peak_gauge,
    bcftools_mem_avg_gauge,
    bcftools_walltime_gauge,
    vcf_download_walltime_gauge,
    vcf_download_cpu_gauge,
    vcf_download_mem_peak_gauge,
    vcf_download_mem_avg_gauge,
    task_walltime_gauge,
):
    """Update all Prometheus Gauges with values from the cache."""
    cached = get_all_cached_metrics()

    for (job_id, task_name), value in cached.get("task_cpu_seconds", {}).items():
        task_cpu_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_python_overhead_cpu_seconds", {}).items():
        python_overhead_cpu_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_bcftools_cpu_seconds", {}).items():
        bcftools_cpu_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_memory_peak_bytes", {}).items():
        task_mem_peak_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_memory_avg_bytes", {}).items():
        task_mem_avg_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_bcftools_memory_peak_bytes", {}).items():
        bcftools_mem_peak_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_bcftools_memory_avg_bytes", {}).items():
        bcftools_mem_avg_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_walltime_seconds", {}).items():
        task_walltime_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_bcftools_walltime_seconds", {}).items():
        bcftools_walltime_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_vcf_download_walltime_seconds", {}).items():
        vcf_download_walltime_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_vcf_download_cpu_seconds", {}).items():
        vcf_download_cpu_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_vcf_download_memory_peak_bytes", {}).items():
        vcf_download_mem_peak_gauge.labels(job_id=job_id, task_name=task_name).set(value)

    for (job_id, task_name), value in cached.get("task_vcf_download_memory_avg_bytes", {}).items():
        vcf_download_mem_avg_gauge.labels(job_id=job_id, task_name=task_name).set(value)


def metrics_purge_loop():
    """Background thread that periodically purges old metrics."""
    while True:
        time.sleep(60)
        purge_old_metrics()


_metrics_server_started = False
_metrics_lock = threading.Lock()


def start_metrics_server(port=8101):
    """
    Start the Prometheus-client metrics server and metrics collection. Designed to handle celery prefork concurrency >=1.

    Only the first worker process to call this will start the HTTP server.
    All worker processes will start their own metrics collection threads.

    Subsequent processes will gracefully skip the HTTP server initialization but
    continue collecting metrics. All metrics are aggregated and served via the
    single HTTP endpoint.

    Celery prefork concurrency=1 is needed both for docker compose and k8s deployments.
    For local dev with docker compose, it is needed to get the per-task CPU and RAM metrics to work correctly.

    For k8s deployment, concurrency must be set to 1 (1 worker process per pod). Scaling will be handled by increasing the number of pods.
    Example with concurrency=1: 1 worker process, 1 collection thread, 1 HTTP server.
    """
    if not ENABLE_WORKER_METRICS:
        logger.info("Metrics collection disabled via ENABLE_WORKER_METRICS environment variable")
        return
    global _metrics_server_started

    pid = os.getpid()
    worker_info.info({"worker_name": WORKER_NAME, "pid": str(pid)})

    # Import here to avoid circular dependency
    from divbase_api.services.queries import BcftoolsQueryManager

    bcftools_monitoring_config.info(
        {"enabled": str(BcftoolsQueryManager.ENABLE_SUBPROCESS_MONITORING), "sample_interval": "0.01s"}
    )
    logger.info(f"Bcftools subprocess monitoring: {BcftoolsQueryManager.ENABLE_SUBPROCESS_MONITORING}")

    thread = threading.Thread(target=collect_system_metrics, daemon=True)
    thread.start()
    logger.info(f"Metrics collection started for PID {pid}")

    with _metrics_lock:
        if not _metrics_server_started:
            try:
                start_http_server(port)
                _metrics_server_started = True
                logger.info(f"Metrics HTTP server started on port {port}")
                logger.info(f"Metrics cache TTL: {TASK_METRICS_CACHE_TTL_MINUTES} minutes")

                # Start background thread for purging old metrics
                purge_thread = threading.Thread(target=metrics_purge_loop, daemon=True)
                purge_thread.start()
                logger.info("Metrics purge thread started")
            except OSError as e:
                if e.errno == 98:  # Address already in use
                    logger.info(f"Metrics server already running on port {port} (started by another process)")
                else:
                    raise
