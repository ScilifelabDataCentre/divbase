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
# Per-task metrics
task_cpu_seconds = Gauge("celery_task_cpu_seconds_total", "CPU seconds used by task", ["task_id", "task_name"])
task_memory_bytes = Gauge("celery_task_memory_bytes", "RAM used by task", ["task_id", "task_name"])
task_bcftools_step_only_cpu_seconds = Gauge(
    "celery_task_bcftools_cpu_seconds_total",
    "CPU seconds used by bcftools subprocess in task",
    ["task_id", "task_name"],
)
task_bcftools_step_only_memory_bytes = Gauge(
    "celery_task_bcftools_memory_bytes", "RAM used by bcftools subprocess in task", ["task_id", "task_name"]
)

# Metrics cache for per-task metrics to persist across multiple task executions
# Structure: {metric_name: {(task_id, task_name): {"value": float, "timestamp": datetime}}}
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


def store_task_metric(metric_name: str, task_id: str, task_name: str, value: float):
    """Store a task metric in the cache with current timestamp."""
    with metrics_cache_lock:
        task_metrics_cache[metric_name][(task_id, task_name)] = {
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
            ("task_memory_bytes", task_memory_bytes),
            ("task_bcftools_cpu_seconds", task_bcftools_step_only_cpu_seconds),
            ("task_bcftools_memory_bytes", task_bcftools_step_only_memory_bytes),
        ]:
            tasks_to_remove = [
                key for key, data in task_metrics_cache[metric_name].items() if data["timestamp"] < cutoff_time
            ]
            for key in tasks_to_remove:
                del task_metrics_cache[metric_name][key]
                gauge.remove(*key)
                logger.debug(f"Purged old metric: {metric_name} for task_id={key[0]} from Prometheus client memory")


def get_all_cached_metrics():
    """Get all cached metrics for exposure via Prometheus."""
    with metrics_cache_lock:
        return {
            metric_name: {key: data["value"] for key, data in tasks.items()}
            for metric_name, tasks in task_metrics_cache.items()
        }


def update_prometheus_gauges_from_cache(task_cpu_gauge, task_mem_gauge, bcftools_cpu_gauge, bcftools_mem_gauge):
    """Update all Prometheus Gauges with values from the cache."""
    cached = get_all_cached_metrics()

    for (task_id, task_name), value in cached.get("task_cpu_seconds", {}).items():
        task_cpu_gauge.labels(task_id=task_id, task_name=task_name).set(value)

    for (task_id, task_name), value in cached.get("task_memory_bytes", {}).items():
        task_mem_gauge.labels(task_id=task_id, task_name=task_name).set(value)

    for (task_id, task_name), value in cached.get("task_bcftools_cpu_seconds", {}).items():
        bcftools_cpu_gauge.labels(task_id=task_id, task_name=task_name).set(value)

    for (task_id, task_name), value in cached.get("task_bcftools_memory_bytes", {}).items():
        bcftools_mem_gauge.labels(task_id=task_id, task_name=task_name).set(value)


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
    global _metrics_server_started

    pid = os.getpid()
    worker_info.info({"worker_name": WORKER_NAME, "pid": str(pid)})
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
