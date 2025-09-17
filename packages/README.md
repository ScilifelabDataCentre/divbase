DivBase contains 3 packages for end use, and 1 library module:

1. divbase-cli - A command-line interface for users to interact with the a DivBase instance
2. divbase-api - An API deployed on the server which handles all communication between the user, job system and S3. 
3. divbase-worker - A celery worker deployed on the server that performs long running tasks.
4. divbase-lib - A library that provides core functionalities and utilities shared between 2 or more of the divbase packages. 

