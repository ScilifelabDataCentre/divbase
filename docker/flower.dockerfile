## Stage 1: Build 
FROM python:3.12.11-alpine3.22  AS builder

WORKDIR /app

RUN apk add --no-cache curl gcc g++ musl-dev && \
    pip install --upgrade pip

# Pip will complain if the readme is not copied over, since it is referenced in pyproject.toml  
COPY README.md ./

# Copy all package sources and install in dependency order
COPY packages/divbase-lib/ ./packages/divbase-lib/
RUN pip install -e ./packages/divbase-lib/
COPY packages/divbase-worker/ ./packages/divbase-worker/
RUN pip install -e ./packages/divbase-worker/

## Stage 2: Final stage
FROM python:3.12.11-alpine3.22 

WORKDIR /app

# curl is needed for healthchecks
RUN apk add --no-cache curl

# for use together with pip editable install in the builder stage
COPY --from=builder /app/packages/divbase-lib/ ./packages/divbase-lib/
COPY --from=builder /app/packages/divbase-worker/ ./packages/divbase-worker/

# If pip editable installs are used, the .pth file is copied over in this step. 
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin


ENTRYPOINT ["celery", "-A", "divbase_worker.tasks"]