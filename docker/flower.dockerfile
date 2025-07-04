FROM python:3.12-alpine

WORKDIR /app

# Install curl for healthcheck
RUN apk add --no-cache curl

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./
COPY src/ ./src/

RUN pip install --upgrade pip && pip install -e .

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks"]