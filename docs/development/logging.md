---
search:
  exclude: true
---

# Logging

DivBase API and worker use [structlog](https://www.structlog.org/en/stable/). divbase-lib and divbase-cli use the standard library logging module.

Structlog used so we can output logs as JSON in deployed environments for easier to parse logs with e.g. grafana alloy. In `local_dev` and `test` environments, we render logs as text and not as JSON, so easier to read.

We don't use structlog in the cli or lib modules as would add unnecessary bloat.

## Configuration

Shared logging configuration lives in `logging_config.py`. The config can handle both:

- structlog logger calls.
- non-structlog stdlib logging calls (for example from third-party libraries like uvicorn).

All logs go to stdout only. Successful (200) healthcheck requests are filtered out of uvicorn's access logs to reduce noise. A failed healthcheck is still logged.

## Special context added to logs

1. Request context (API)

    - We use middleware to create and bind (aka add) a uuid to all logging events associated with a specific request. This uuid is then present in all logs covering the entire request-response cycle - so easy to filter.
    - The uuid is included as a response header ("X-Request-ID") and divbase-cli includes the response header in the user facing error message if an error occurs, so we can use this for debugging.

2. Task context (worker)

    - celery signals used to bind (aka add) the task id to all logging events associated with a specific task, - again same idea as above, nice for filtering.
