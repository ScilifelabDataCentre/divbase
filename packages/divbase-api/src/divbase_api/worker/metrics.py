import os
import socket
import threading
import time

import psutil
from prometheus_client import Gauge, Info, start_http_server

# Define Prometheus metrics. Use the same prefix to simplify discovery in Prometheus and Grafana UIs.
worker_cpu_percent = Gauge("celery_worker_prom_client_cpu_percent", "CPU usage percent", ["worker_name"])
worker_memory_mb = Gauge("celery_worker_prom_client_memory_mb", "Memory usage in MB", ["worker_name"])
worker_memory_percent = Gauge("celery_worker_prom_client_memory_percent", "Memory usage percent", ["worker_name"])
worker_info = Info("celery_worker_prom_client_info", "Worker information")
worker_cpu_time_total = Gauge(
    "celery_worker_prom_client_cpu_time_total", "Total CPU time (user + system) in seconds", ["worker_name"]
)
worker_memory_bytes = Gauge("celery_worker_prom_client_memory_bytes", "Memory usage in bytes", ["worker_name"])
worker_num_threads = Gauge("celery_worker_prom_client_num_threads", "Number of threads", ["worker_name"])
worker_open_fds = Gauge("celery_worker_prom_client_open_fds", "Number of open file descriptors", ["worker_name"])

WORKER_NAME = socket.gethostname()


def collect_system_metrics():
    process = psutil.Process()
    while True:
        try:
            cpu = process.cpu_percent(interval=1)
            worker_cpu_percent.labels(worker_name=WORKER_NAME).set(cpu)
            mem_info = process.memory_info()
            mem_mb = mem_info.rss / 1024 / 1024
            worker_memory_mb.labels(worker_name=WORKER_NAME).set(mem_mb)
            mem_percent = process.memory_percent()
            worker_memory_percent.labels(worker_name=WORKER_NAME).set(mem_percent)

            cpu_times = process.cpu_times()
            worker_cpu_time_total.labels(worker_name=WORKER_NAME).set(cpu_times.user + cpu_times.system)
            worker_memory_bytes.labels(worker_name=WORKER_NAME).set(mem_info.rss)
            worker_num_threads.labels(worker_name=WORKER_NAME).set(process.num_threads())
            # num_fds() is only available on Unix. if open fds increase it might be indicative of a a resource leak
            if hasattr(process, "num_fds"):
                worker_open_fds.labels(worker_name=WORKER_NAME).set(process.num_fds())
        except Exception as e:
            print(f"Error collecting metrics: {e}")
        time.sleep(5)


def start_metrics_server(port=8001):
    worker_info.info({"worker_name": WORKER_NAME, "pid": str(os.getpid())})
    start_http_server(port)
    print(f"Metrics server started on port {port}")
    thread = threading.Thread(target=collect_system_metrics, daemon=True)
    thread.start()
