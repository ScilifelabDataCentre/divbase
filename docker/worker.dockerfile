FROM ghcr.io/astral-sh/uv:python3.13-alpine

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


ENV UV_COMPILE_BYTECODE=1
ENV UV_LINK_MODE=copy

RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev

# watch requires that the container's USER can write to the target path so it can update files. 
COPY --chown=app:app . /app

RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --locked --no-dev

ENV PATH="/app/.venv/bin:$PATH"

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks", "worker", "--loglevel=info"]