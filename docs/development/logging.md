# Logging

DivBase API and worker use [structlog](https://www.structlog.org/en/stable/). divbase-lib and divbase-cli use the standard library logging module.

Structlog used so we can output logs as JSON in deployed environments for easier to parse logs with e.g. grafana alloy. In `local_dev` and `test` environments, we render logs not as JSON, so easier to read.

We don't use structlog in the cli or lib modules as would add unnecessary bloat.

## Configuration

Shared logging configuration lives in logging_config.py. The config can handle both:

- structlog logger calls.
- non-structlog stdlib logging calls (for example from third-party libraries like uvicorn).

## Special context added to logs

1. Request context (API)

    - We use middleware to create and bind (aka add) a uuid to all logging events associated with a specific request. This uuid is then present in all logs covering the entire request-response cycle - so easy to filter.
    - The uuid is included as a response header ("X-Request-ID") and divbase-cli includes the response header in the user facing error message if an error occurs, so we can use this for debugging.

2. Task context (worker)

    - celery signals used to bind (aka add) the task id to all logging events associated with a specific task, - again same idea as above, nice for filtering.

## Logging to file (by default turned off)

File logging is optional (for api + worker) and controlled by setting `LOG_TO_FILE=1` and is by default turned off.

- If on, it will write logs to /logs/ in deployed environments (as well as stdout).
- If on, it will write logs to `{REPO_ROOT}/docker/logs/` in local dev (as well as stdout).

There is a celery cron job in place to clean up old logs files (older than 30 days).

At time of writing, in deployed environments, logs are written to a shared PVC and stdout. In local dev logs are only written to stdout. At some point it would be nice to setup loki, grafana alloy and grafana for log aggregation and visualisation in the deployed environments. If that happens then the logging to file approch can most likely be deleted.
