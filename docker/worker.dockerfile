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

# Copy all package sources
COPY packages/divbase-lib/ ./packages/divbase-lib/
COPY packages/divbase-worker/ ./packages/divbase-worker/

# Install packages in dependency order
RUN pip install -e ./packages/divbase-lib/
RUN pip install -e ./packages/divbase-worker/

ENTRYPOINT ["celery", "-A", "divbase_worker.tasks", "worker", "--loglevel=info"]