# Deployment of DivBase Server

- DivBase is deployed on a kubernetes cluster and its deployment is managed in our [private repository, argocd-divbase](https://github.com/ScilifelabDataCentre/argocd-divbase), see there for details.

- In this repository the main thing that needs to be done is to build and push the new docker images for the API and worker.

- We have a GitHub Actions workflow that does this (`.github/workflows/publish_images.yaml`). You can trigger a workflow_dispatch to manually create the docker images at a given commit.

- Head to the actions tag on GitHub and to the action "Build and push both docker images to the GitHub Container Registry". From there click run manual workflow. You can choose to specify the name of the image tag if you want. Otherwise leave the input blank and it will be tagged with the full commit hash.

## Deployment of DivBase CLI

TODO - when the time comes add the info on how to deploy/publish new versions of the CLI package to PyPI.
