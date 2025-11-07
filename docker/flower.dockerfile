## Stage 1: Build 
FROM python:3.12.11-alpine3.22  AS builder

WORKDIR /app

RUN apk add --no-cache curl gcc g++ musl-dev && \
    pip install --upgrade pip

# Pip will complain if the readme is not copied over, since it is referenced in pyproject.toml 
COPY README.md ./

# Copy all package sources and install in dependency order
COPY packages/divbase-lib/ ./packages/divbase-lib/
RUN pip install ./packages/divbase-lib/
COPY packages/divbase-api/ ./packages/divbase-api/
RUN pip install ./packages/divbase-api/

## Stage 2: Final stage
FROM python:3.12.11-alpine3.22 

WORKDIR /app

# curl is needed for healthchecks
RUN apk add --no-cache curl

COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

ENTRYPOINT ["celery", "-A", "divbase_api.worker.tasks"]