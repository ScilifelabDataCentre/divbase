# Monitoring and Observability in DivBase

Monitoring and observability are related concepts. To paraphrase [this guide](https://www.ibm.com/think/topics/observability-vs-monitoring): monitoring tells you when something is wrong; observability helps you understand why it’s wrong and how to fix it. In short, monitoring is about assessing system health and performance, and observability is about gaining understanding of the internal state of a system based on its outputs. To do this, monitoring and observability both rely on three types of data: metrics, logs, and traces. This document describes how these concepts are implemented in DivBase, and to which degree.

## Implemented aspects of monitoring and observability

At the time of writing, a full monitoring/observability stack has not been implemented for DivBase. This section outlines the different subsystems and topics that are currently implemented.

### CPU and RAM

CPU and RAM metrics can be collected for the Celery workers, as is descrcibed in detail in [Monitoring: Celery Worker Metrics](worker_metrics.md). There are two levels: **General (container-level) metrics** — CPU utilisation, RSS memory, thread count, and open file descriptors — are continuously sampled every ~6 seconds per worker process regardless of whether a task is running. **Per-task metrics** provide finer granularity: wall time, CPU time, and peak/average RSS are measured for the the S3 VCF download step, the `bcftools` subprocess calls, and the remaining Python worker process overhead during the task. These metrics are stored in a short-lived in-memory cache (default TTL: 5 minutes) exposed by the worker container in `/metrics` endpoint for Prometheus to scrape. Data can be fetched by PromQL queries, as explained in [Section 3.3 in Monitoring: Celery Worker Metrics](worker_metrics.md#33-fetching-the-scraped-data-from-the-prometheus-database).

### RabbitMQ queue

RabbitMQ metrics can be scraped using the built-in Prometheus plugin included with the `management-alpine` image. The plugin must be explicitly enabled via two environment variables: `RABBITMQ_SERVER_ADDITIONAL_ERL_ARGS=-rabbitmq_prometheus true` (enables the plugin) and `RABBITMQ_PROMETHEUS_RETURN_PER_OBJECT_METRICS=true` (enables per-queue metrics at the `/metrics/per-object` path). For the local Docker Compose stack, Prometheus scrapes that endpoint at port `15692` using the non-default path `/metrics/per-object` (as configured in `docker/prometheus.yml`).

In the future, this metric could be of interest for setting up alerting based on the queue depth and status.

### Grafana dashboards

Custom DivBase dashboards are configured at `docker/grafana` and are made available at `http://localhost:3000` when the local monitoring stack is running (see [Monitoring: Celery Worker Metrics — Section 3.1](worker_metrics.md#31-local-monitoring-stack-docker-compose) for how to start it). The follwing dashboards have been configured:

| Dashboard | Source | What it shows |
|---|---|---|
| **Celery Worker Prometheus Client Metrics** | `celery_worker_prom_client_*` gauges (port 8101) | CPU utilisation (%), memory usage (bytes and %), open file descriptors, thread count, and cumulative CPU time — one time series per worker process |
| **Container CPU & Memory Usage (cAdvisor)** | cAdvisor (port 8080) | CPU cores/s and memory working-set bytes for all running containers — useful for container-level comparison independent of the custom metrics server |
| **RabbitMQ Overview** | RabbitMQ Prometheus plugin (port 15692) | Messages published per minute; current queue depth and consumer count for the `celery`, `quick`, and `long` queues |
