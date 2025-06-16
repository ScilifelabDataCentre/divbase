FROM python:3.12-alpine

ARG BCFTOOLS_VERSION="1.22"

RUN apk add --no-cache \
    build-base \
    curl \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    perl-dev 

RUN curl -fsSL https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    | tar -C /tmp -xjf- \
    && cd /tmp/bcftools-${BCFTOOLS_VERSION} \
    && make \
    && make install \
    && cd - && rm -rf /tmp/bcftools-${BCFTOOLS_VERSION}

WORKDIR /app

CMD ["tail", "-f", "/dev/null"]