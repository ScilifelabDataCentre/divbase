#!/bin/bash
#
# Wrapper script that generates a mock VCF file using the fake-vcf package in a Docker container, and copies the resulting file to the host machine. 
# The number of samples and variants to be generated in the mock file can be specified as arguments.
# The random seed is fixed to ensure that the same mock file is generated each time for the same input parameters.
# 
# Usage: 
# bash scripts/benchmarking/generate_mock_vcf.sh -s <number_of_samples_to_generate> -r <number_of_variants_to_generate>
#
# Example:
# bash scripts/benchmarking/generate_mock_vcf.sh -s 1000 -r 50

set -e


SAMPLES=5000
VARIANTS=100

while [[ $# -gt 0 ]]; do
  case $1 in
    -s)
      SAMPLES="$2"
      shift 2
      ;;
    -r)
      VARIANTS="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      echo "Usage: $0 -s <number_of_samples_to_generate> -r <number_of_variants_to_generate>"
      exit 1
      ;;
  esac
done

IMAGE="benchmarking-fake-vcf"
CONTAINER_NAME="fake-vcf-gen"

# Build the Docker image
docker build -f docker/benchmarking.dockerfile -t $IMAGE .

docker rm -f $CONTAINER_NAME 2>/dev/null || true

docker run --name $CONTAINER_NAME -d $IMAGE tail -f /dev/null

OUTFILE="/tmp/mock_vcf_${SAMPLES}s_${VARIANTS}r.vcf.gz"
HOST_OUTFILE="mock_vcf_${SAMPLES}s_${VARIANTS}r.vcf.gz"

docker exec $CONTAINER_NAME bash -c "cd /opt/fake-vcf && poetry run fake-vcf generate -s $SAMPLES -r $VARIANTS --seed 12345 | bgzip > $OUTFILE"

# Copy the file from the container to the host
docker cp $CONTAINER_NAME:$OUTFILE ./$HOST_OUTFILE

docker rm -f $CONTAINER_NAME