FROM python:3.12-alpine

WORKDIR /app

# Install curl for healthcheck
RUN apk add --no-cache curl gcc g++ musl-dev


RUN pip install --upgrade pip 

COPY README.md ./

# Copy all package sources and install in dependency order
COPY packages/divbase-lib/ ./packages/divbase-lib/
RUN pip install -e ./packages/divbase-lib/
COPY packages/divbase-worker/ ./packages/divbase-worker/
RUN pip install -e ./packages/divbase-worker/

ENTRYPOINT ["celery", "-A", "divbase_worker.tasks"]