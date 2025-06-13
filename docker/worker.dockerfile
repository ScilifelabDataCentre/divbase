FROM python:3.12-slim

WORKDIR /app

# Install system updates and dependencies
RUN apt-get update && apt-get upgrade -y && apt-get install -y --no-install-recommends gcc build-essential && rm -rf /var/lib/apt/lists/*

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./
COPY src/ ./src/

RUN pip install --upgrade pip && pip install -e .

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks", "worker", "--loglevel=info"]