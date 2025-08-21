FROM python:3.12-alpine

WORKDIR /app

# Install curl for healthcheck
RUN apk add --no-cache curl gcc g++ musl-dev

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./

RUN pip install --upgrade pip 

COPY src/ ./src/
RUN pip install -e .

ENTRYPOINT ["celery", "-A", "divbase_tools.tasks"]