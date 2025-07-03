FROM python:3.12-alpine

WORKDIR /app

# Install Alpine-specific packages
RUN apk update && \
    apk add --no-cache \
    gcc \
    musl-dev \
    python3-dev \
    build-base \
    docker \
    ca-certificates \
    curl \
    libffi-dev

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./
COPY src/ ./src/

RUN pip install --upgrade pip && pip install -e .

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks"]