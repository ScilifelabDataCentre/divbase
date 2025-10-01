#!/bin/sh
# Note: using getopts here could be a awkward as DOCKER_TAG sometimes exists, sometimes doesn't.
# DOCKER_TAG occurs if triggered by workflow_dispatch and tag added.

EVENT_NAME=$1
REF_NAME=$2
SHA=$3
DOCKER_TAG=$4


if [ "$EVENT_NAME" = "release" ]; then
    tag="$REF_NAME"
else
    tag="${DOCKER_TAG:-$SHA}"
fi
echo "tag=$tag"