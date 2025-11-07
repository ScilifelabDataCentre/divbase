"""
Handles connection between celery workers and the postgresql db.
"""

import logging
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

logger = logging.getLogger(__name__)

WORKER_DATABASE_URL = os.getenv("WORKER_DATABASE_URL", "NOT_SET")
sync_engine = create_engine(WORKER_DATABASE_URL, echo=False)
SyncSessionLocal = sessionmaker(bind=sync_engine)
