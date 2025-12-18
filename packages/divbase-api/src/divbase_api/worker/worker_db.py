"""
Handles connection between celery workers and the postgresql db.
"""

import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# NOTE: As SyncSessionLocal is created at module load time,
# both the API and WORKER applications will import an instance of the SyncSessionLocal object.
# Both already need the enviroment variable SYNC_DATABASE_URL set.
# And as the overhead of creating the SyncSessionLocal object is low, we can just create it for both,
# even if the API will not actually use it.
SYNC_DATABASE_URL = os.environ["SYNC_DATABASE_URL"]
sync_engine = create_engine(SYNC_DATABASE_URL, echo=False, pool_pre_ping=True)
SyncSessionLocal = sessionmaker(bind=sync_engine)
