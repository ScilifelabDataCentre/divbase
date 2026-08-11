"""
Handles connection between celery workers and the postgresql db.
"""

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from divbase_api.worker.worker_config import worker_settings

# NOTE: As SyncSessionLocal is created at module load time,
# both the API and WORKER applications will import an instance of the SyncSessionLocal object.
# Both already need the environment variable SYNC_DATABASE_URL set.
# And as the overhead of creating the SyncSessionLocal object is low, we can just create it for both,
# even if the API will not actually use it.
sync_engine = create_engine(worker_settings.general.sync_url.get_secret_value(), echo=False, pool_pre_ping=True)
SyncSessionLocal = sessionmaker(bind=sync_engine)
