# DivBase API

!!! attention "For End Users: Use divbase-cli instead"
    **Users are strongly encouraged to use divbase-cli rather than calling the API endpoints directly.**

    divbase-cli offers several advantages over direct API usage, including:
    - Handles authentication logic, including automatic token refresh
    - Simplifies commands and workflows (e.g. file uploads/downloads require working with presigned URLs).
    - Better error messages and user experience

    If there is something you cannot do with divbase-cli that you think should be possible, please let us know.
---

## Overview

- The API is written with fastapi and is used by divbase-cli. It is not designed for direct end user usage but we should still try to make the endpoints user-friendly/understandable.

- The API and frontend are served from the same origin, so to differentiate we serve all API routes under the `/api/v1/` path.

- You can see the 2 versions of the automatically generated api docs by going to the `api/v1/docs` or `api/v1/redoc` routes. The API docs are publically available.
