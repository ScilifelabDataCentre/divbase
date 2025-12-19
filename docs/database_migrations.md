# Database Migrations with Alembic

The directory `packages/divbase-api/src/divbase_api/migrations` contains database migrations for the DivBase postgres db instance. We use [Alembic](https://alembic.sqlalchemy.org/) to manage database schema changes.

## Overview

- **Migration files/scripts** are stored in `versions/` with descriptive names: `YYYY-MM-DD_description.py`
- **Migrations are generated** by comparing SQLAlchemy models with the current database instance deployed using the local docker compose stack.
- **Migrations are applied automatically** in local development when the docker compose stack starts up using an init-container (db-migrator).
- **In production/deployed environments**, migrations should be applied as part of the deployment process - multi step and not automatic, see below.
- **In deployed environments with actual user data**, database migrations are a dangerous operation that can lead to data loss if not done carefully. A backup of the database must be taken before applying migrations in production.
- **FastAPI's lifespan event automatically checks** if migrations are up to date during startup. - will raise an error if not.
- **`pytest test/migrations`** can run tests to ensure all migrations can be applied cleanly to a fresh database.
- **Celery managed tables** (`celery_taskmeta` and `celery_groupmeta`) are excluded from alembic (see the `migrations/env.py`), as Celery handles their creation and updates. Keep this in mind when creating migrations.

### Creating a New Migration

Note: The FastAPI lifespan event checks (but does not apply) that all migration scripts have been applied on startup - aka at the latest version.

This means when you create a new migration script in watch mode, the app will error on restart until the migration is applied.
If you instead make your changes using compose up, then make sure to restart the stack before generating the migration script (step 2) so the changes to your db models are actually included.

The fastapi lifespan event is a little awkward yes, but it is primarily for protecting the deployed instances.

#### 1. Make changes to the SQLAlchemy models

If you add a new model, include it in the `packages/divbase-api/src/divbase_api/models/__init__.py` so Alembic can detect it.

Make sure the divbase stack is running:

```bash
docker compose -f docker/divbase_compose.yaml down && docker compose -f docker/divbase_compose.yaml watch
```

#### 2. Generate the migration script by entering the FastAPI container

```bash
# Enter the running FastAPI container
docker compose -f docker/divbase_compose.yaml exec -it fastapi sh

# Generate migration (use descriptive names)
alembic revision --autogenerate -m "write_your_useful_slug_here"

# Exit container - if needed - see NOTE below
exit
```

**NOTE: If you've been using `docker compose watch` up to this point, it is expected that now the FastAPI service will be in a crashed state due to pending migrations -> the lifespan event will fail when the app restarts.**
**NOTE2** The commands can be run as a one-liner with:
`docker compose -f docker/divbase_compose.yaml exec -it fastapi alembic revision --autogenerate -m "write_your_useful_slug_here"`

#### 3. Review the generated migration file

- The migration script created in the container is automatically synced to your host machine (`packages/divbase-api/src/divbase_api/migrations/versions/`)
- Open the created migration script and check that the auto generated changes correctly represent your intended changes
- Be extra careful if you rename a table, or if you are dropping columns or tables.
- Modify the `upgrade()` and `downgrade()` functions if needed to actually achieve the desired schema changes.
- If you have created a new Enum you need to append a drop statement into the downgrade function for the enum.

```python
# rest of downgrade function
# where projectroles is the lower case name of the enum you created
op.execute("DROP TYPE IF EXISTS projectroles")
```

See the bottom of the first migration script (2025-12-04_initial_migration.py) for an example of this.

#### 4. Test out the migration

- In local development, migrations are applied automatically when the stack starts up using the `db-migrator` init-container.
- The `db-migrator` ensures that migrations are applied before starting FastAPI and worker services (via `condition: service_completed_successfully`).
- Do not down the volumes (e.g. `docker compose -f docker/divbase_compose.yaml down -v`), this is a good initial test to ensure the migration works as expected.
- To test the migration:

```bash
# interrupt if you have a process already running.
docker compose -f docker/divbase_compose.yaml down && docker compose -f docker/divbase_compose.yaml watch
```

This will make the db-migrator init-container run the migrations, and then start fastapi and workers if successful.
If issues check the db-migrator logs.

You can furthermore check that an update to a table schema has been applied as intended by running the postgreSQL `\d <TABLE_NAME>` command to inspect how the table looks like in the database engine after the restart:
`docker exec -it divbase-postgres-1 psql -U divbase_user -d divbase_db -c '\d "<TABLE_NAME>";'`

To check the current status of applied migrations in the database, you can run:

```bash
docker exec -it divbase-postgres-1 psql -U divbase_user -d divbase_db -c 'SELECT * FROM "alembic_version";'
```

You can also run the `pytest-alembic` tests to further validate the newly created migration script.

```bash
pytest tests/migrations
```

## Production Deployment

TODO

## Troubleshooting

### "No changes in schema detected"

- Ensure your models are properly imported in `models/base.py`
- Check that your model changes are actually different from the database.

### Migration fails to apply

- Check the migration file for syntax errors
- Ensure the database is in the expected state
- Use `alembic current` to check current migration status or inspect the `alembic_version` table in the database (see section 4).

### Starlette admin/admin panel not showing the models

- You need the models in your src to match the postgres schema. So if you have pending changes (that you have or have not created migrations for) they won't display until you've actually done the migration.
