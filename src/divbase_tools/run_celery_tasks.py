import time

from divbase_tools.tasks import add

print("Sending task...")
result = add.delay(2, 3)
print(f"Task ID: {result.id}")

print("Waiting for task completion...")

while not result.ready():
    print("Task still processing...")
    time.sleep(1)
print(f"Result: {result.get()}")
