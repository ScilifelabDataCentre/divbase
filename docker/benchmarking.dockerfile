FROM python:3.12-slim

RUN apt-get update && apt-get install -y git make && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/endast/fake-vcf.git /opt/fake-vcf
WORKDIR /opt/fake-vcf
RUN pip install poetry biopython
RUN pip install .

RUN make poetry-download \
    && poetry self add poetry-plugin-export \
    && make install-all

WORKDIR /app
COPY . /app
