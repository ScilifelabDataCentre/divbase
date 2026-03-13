"""
Constants used by multiple CRUD modules, to avoid circular imports.
Can be especially useful for constants used by both API and worker CRUD operations, since these are in the same package and more likely to have circular import issues.
"""

PROVISIONAL_ENTRY_TTL_SECONDS = 120

CELERYTASKMETA_ENTRY_GAP_TTL_SECONDS = 3600  # 1 hour in the queue without being picked up by a worker

ACTIVE_CELERY_STATUSES = {"PENDING", "STARTED", "RETRY"}
