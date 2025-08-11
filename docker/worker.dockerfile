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

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./

RUN pip install --upgrade pip 

COPY src/ ./src/
RUN pip install -e .

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks", "worker", "--loglevel=info"]