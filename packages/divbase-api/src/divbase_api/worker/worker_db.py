"""
Handles connection between celery workers and the postgresql db.
"""

import logging
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

logger = logging.getLogger(__name__)

# if block is used as API and worker are in same package and
# API imports modules from the worker (e.g. tasks.py).
# Without the if:
# - Would have to add the WORKER_DATABASE_URL environment variable to the API and
# - create the SyncSessionLocal object that would not used by the API.
# TODO - consider creating SyncSessionLocal object using https://docs.celeryq.dev/en/stable/userguide/signals.html#worker-init instead?
WORKER_DATABASE_URL = os.getenv("WORKER_DATABASE_URL")
if WORKER_DATABASE_URL:
    sync_engine = create_engine(WORKER_DATABASE_URL, echo=False, pool_pre_ping=True)
    SyncSessionLocal = sessionmaker(bind=sync_engine)
else:
    SyncSessionLocal = None
