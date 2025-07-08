FROM python:3.12-alpine

WORKDIR /app

# Install curl for healthcheck
RUN apk add --no-cache curl

# copy readme to avoid pip complaining about missing files
COPY pyproject.toml README.md ./

RUN pip install --upgrade pip 

COPY src/ ./src/
RUN pip install -e .

# host needs to be set to 0.0.0.0 to be accessible from outside the container
CMD ["uvicorn", "divbase_tools.divbase_api:app", "--host", "0.0.0.0", "--port", "8000", "--reload", "--reload-dir", "src/divbase_tools"]