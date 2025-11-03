## Stage 1: Build
FROM python:3.12.11-alpine3.22 AS builder

WORKDIR /app

ARG BCFTOOLS_VERSION="1.22"

RUN apk update && \
    apk add --no-cache \
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
    && cd - && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}  \
    && pip install --upgrade pip

# Pip will complain if the readme is not copied over, since it is referenced in pyproject.toml   
COPY README.md ./

# Copy all package sources and install in dependency order
COPY packages/divbase-lib/ ./packages/divbase-lib/
RUN pip install ./packages/divbase-lib/
COPY packages/divbase-api/ ./packages/divbase-api/
RUN pip install ./packages/divbase-api/


## Stage 2: Final image
FROM python:3.12.11-alpine3.22

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

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# Create a proper user and group to avoid Celery warnings. Write access to /app is needed for bcftools.
RUN addgroup -g 1000 appuser && \
    adduser -u 1000 -G appuser -s /bin/sh -D appuser && \
    chown -R appuser:appuser /app

USER appuser

ENTRYPOINT ["celery", "-A", "divbase_api.worker.tasks", "worker", "--loglevel=info"]