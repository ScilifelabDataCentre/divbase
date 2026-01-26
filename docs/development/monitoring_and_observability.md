# Monitoring and Observability in DivBase

Monitoring and observability are related concepts (for an overview, see e.g. [this guide](https://www.ibm.com/think/topics/observability-vs-monitoring)). In short, monitoring is about assessing system health and performance, and observeability is about gaining understanding of the internal state of a system based on its outputs. To do this, monitoring and observability both rely on three types of data: metrics, logs, and traces. This document describes how these concepts are implemented in DivBase, and to which degree.

## Metrics

### CPU and RAM

A per-task resource monitoring system for the Celery workers is implemented using a custom `prometheus-client` metrics server and `psutils`. For details on this see [Monitoring: Celery Worker Metrics](worker_metrics.md).

### RabbitMQ queue

Scraped using the build-it Prometheus plugin that comes with RabbitMQ `management-alpine` image

TODO add more details on available metrics, PromQL queries, and Grafana dashboards.
