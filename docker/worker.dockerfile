FROM python:3.12-alpine

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
    && cd - && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}

RUN pip install --upgrade pip 

COPY README.md ./

# Copy all package sources and install in dependency order
COPY packages/divbase-lib/ ./packages/divbase-lib/
RUN pip install -e ./packages/divbase-lib/
COPY packages/divbase-worker/ ./packages/divbase-worker/
RUN pip install -e ./packages/divbase-worker/

# Create a proper user and group to avoid Celery warnings. Write access to /app is needed for bcftools.
RUN addgroup -g 1000 appuser && \
    adduser -u 1000 -G appuser -s /bin/sh -D appuser && \
    chown -R appuser:appuser /app

USER appuser

ENTRYPOINT ["celery", "-A", "divbase_worker.tasks", "worker", "--loglevel=info"]