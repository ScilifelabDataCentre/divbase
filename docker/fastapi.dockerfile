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
COPY packages/divbase-api/ ./packages/divbase-api/
RUN pip install -e ./packages/divbase-api/

## Stage 2: Final stage
FROM python:3.12.11-alpine3.22 

WORKDIR /app

# curl is needed for healthchecks
RUN apk add --no-cache curl

# for use together with pip editable install in the builder stage
COPY --from=builder /app/packages/divbase-lib/ ./packages/divbase-lib/
COPY --from=builder /app/packages/divbase-worker/ ./packages/divbase-worker/
COPY --from=builder /app/packages/divbase-api/ ./packages/divbase-api/

# If pip editable installs are used, the .pth file is copied over in this step. 
COPY --from=builder /usr/local/lib/python3.12/site-packages /usr/local/lib/python3.12/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# host needs to be set to 0.0.0.0 to be accessible from outside the container
CMD ["fastapi", "run", "--host", "0.0.0.0", "/app/packages/divbase-api/src/divbase_api/divbase_api.py"]