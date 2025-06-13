"""
The API server for DivBase.
"""

import uvicorn
from fastapi import FastAPI

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello World"}


def main():
    uvicorn.run("divbase_tools.divbase_api:app", host="127.0.0.1", port=8000, reload=True)


if __name__ == "__main__":
    main()
