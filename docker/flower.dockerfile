FROM python:3.12-alpine

WORKDIR /app

# Install curl for healthcheck
RUN apk add --no-cache curl gcc g++ musl-dev


RUN pip install --upgrade pip 

# Copy workspace files, TODO - do I need outer pyproject.toml?
COPY pyproject.toml ./
COPY README.md ./

# Copy all package sources
COPY packages/divbase-lib/ ./packages/divbase-lib/
COPY packages/divbase-worker/ ./packages/divbase-worker/

# Install packages in dependency order
RUN pip install -e ./packages/divbase-lib/
RUN pip install -e ./packages/divbase-worker/

ENTRYPOINT ["celery", "-A", "divbase_worker.tasks"]