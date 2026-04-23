## Stage 1: Build 
FROM ghcr.io/astral-sh/uv:python3.13-alpine3.23 AS builder


# UV_COMPILE_BYTECODE=1 Compile to bytecode during install so faster startup time
# UV_LINK_MODE=copy recommended for docker 
# UV_PYTHON_DOWNLOADS=0 use images python install 
ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=0 

WORKDIR /app

# This installs dependencies but not divbase packages so layer only invalidated when dependencies change, not source code
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=/app/uv.lock \
    --mount=type=bind,source=pyproject.toml,target=/app/pyproject.toml \
    --mount=type=bind,source=packages/divbase-lib/pyproject.toml,target=/app/packages/divbase-lib/pyproject.toml \
    --mount=type=bind,source=packages/divbase-api/pyproject.toml,target=/app/packages/divbase-api/pyproject.toml \
    uv sync --frozen --no-dev --no-install-workspace --package divbase-api

# Copy source code and install workspace packages
COPY README.md pyproject.toml uv.lock ./
COPY packages/divbase-lib/ ./packages/divbase-lib/
COPY packages/divbase-api/ ./packages/divbase-api/
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-editable --package divbase-api

## Stage 2: Final image (without uv installed)
FROM python:3.13-alpine3.23

WORKDIR /app

# curl needed for healthcheck
RUN apk add --no-cache curl

# Setup a non-root user
RUN addgroup -g 1000 appuser && \
    adduser -u 1000 -G appuser -s /bin/sh -D appuser && \
    chown -R appuser:appuser /app

COPY --from=builder --chown=appuser:appuser /app/.venv /app/.venv

# PYTHONUNBUFFERED=1 allows for log messages to be sent immediately rather than buffered, can lose log message on application crash otherwise
ENV PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1

USER appuser

CMD ["fastapi", "run", "--host", "0.0.0.0", "/app/.venv/lib/python3.13/site-packages/divbase_api/divbase_api.py"]

