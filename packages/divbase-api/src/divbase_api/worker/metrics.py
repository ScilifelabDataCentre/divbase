import os
import socket
import threading
import time

import psutil
from prometheus_client import Gauge, Info, start_http_server

# Define Prometheus metrics. Use the same prefix to simplify discovery in Prometheus and Grafana UIs.
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

WORKER_NAME = socket.gethostname()


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
            print(f"Error collecting metrics: {e}")
        time.sleep(5)


_metrics_server_started = False
_metrics_lock = threading.Lock()


def start_metrics_server(port=8001):
    """
    Start the Prometheus-client metrics server and metrics collection. Designed to handle celery prefork concurrency >=1.

    Only the first worker process to call this will start the HTTP server.
    All worker processes will start their own metrics collection threads.

    Subsequent processes will gracefully skip the HTTP server initialization but
    continue collecting metrics. All metrics are aggregated and served via the
    single HTTP endpoint.

    For local dev with docker compose, any celery concurrency value can be used.
    Example with concurrency=4: 4 worker processes, 4 collection threads, 1 HTTP server.

    For k8s deployment, concurrency must be set to 1 (1 worker process per pod). Scaling will be handled by increasing the number of pods.
    Example with concurrency=1: 1 worker process, 1 collection thread, 1 HTTP server.
    """
    global _metrics_server_started

    pid = os.getpid()
    worker_info.info({"worker_name": WORKER_NAME, "pid": str(pid)})
    thread = threading.Thread(target=collect_system_metrics, daemon=True)
    thread.start()
    print(f"Metrics collection started for PID {pid}")

    with _metrics_lock:
        if not _metrics_server_started:
            try:
                start_http_server(port)
                _metrics_server_started = True
                print(f"Metrics HTTP server started on port {port}")
            except OSError as e:
                if e.errno == 98:  # Address already in use
                    print(f"Metrics server already running on port {port} (started by another process)")
                else:
                    raise
