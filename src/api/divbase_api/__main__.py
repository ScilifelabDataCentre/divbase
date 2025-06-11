import uvicorn
from fastapi import FastAPI

app = FastAPI()


@app.get("/")
async def root():
    return {"message": "Hello World"}


def main():
    uvicorn.run("divbase_api.__main__:app", host="0.0.0.0", port=8000, reload=True)


if __name__ == "__main__":
    main()
