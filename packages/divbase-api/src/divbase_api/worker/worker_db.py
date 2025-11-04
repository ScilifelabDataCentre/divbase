"""
Handles connection between celery workers and the postgresql db.
"""

import logging

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

logger = logging.getLogger(__name__)

WORKER_DATABASE_URL = "postgresql+psycopg://divbase_user:badpassword@postgres:5432/divbase_db"

sync_engine = create_engine(WORKER_DATABASE_URL, echo=True)
SyncSessionLocal = sessionmaker(bind=sync_engine)
