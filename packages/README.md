DivBase is structured as a monorepo and contains 2 packages for end use, and 1 library module:

1. divbase-cli: A command-line interface for users to interact with the a DivBase instance
2. divbase-api: An API deployed on the server which handles all communication between the user, job system and S3. 
    - divbase-worker (submodule called "worker") A celery worker deployed on the server that performs long running tasks.
3. divbase-lib: A library that allows us to share code between the CLI and API packages. One example use case is the API schemas used for request and response models between the API and CLI. 

The API package has 2 entrypoints (fastapi and celery worker) and the CLI package has a single entrypoint (the command line interface).

