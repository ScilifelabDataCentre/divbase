## Stage 1: Build
FROM ghcr.io/astral-sh/uv:python3.13-alpine3.23 AS builder

ARG BCFTOOLS_VERSION="1.22"

ENV UV_COMPILE_BYTECODE=1 \
    UV_LINK_MODE=copy \
    UV_PYTHON_DOWNLOADS=0

WORKDIR /app

RUN apk add --no-cache \
    gcc \
    musl-dev \
    python3-dev \
    build-base \
    ca-certificates \
    curl \
    curl-dev \
    libffi-dev \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    openssl-dev \
    perl-dev

RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    | tar -C /tmp -xjf- \
    && cd /tmp/bcftools-${BCFTOOLS_VERSION} \
    && make \
    && make install \
    && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}

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

# Only install runtime dependencies
RUN apk add --no-cache \
    ca-certificates \
    curl \
    libffi \
    zlib \
    bzip2 \
    xz \
    openssl \
    perl

# Setup a non-root user. Write access to /app is needed for bcftools.
RUN addgroup -g 1000 appuser && \
    adduser -u 1000 -G appuser -s /bin/sh -D appuser && \
    chown -R appuser:appuser /app

COPY --from=builder /usr/local/bin/bcftools /usr/local/bin/bcftools
COPY --from=builder --chown=appuser:appuser /app/.venv /app/.venv

# PYTHONUNBUFFERED=1 allows for log messages to be sent immediately rather than buffered, can lose log message on application crash otherwise
ENV PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1

USER appuser

ENTRYPOINT ["celery", "-A", "divbase_api.worker.tasks", "worker", "--loglevel=info"]