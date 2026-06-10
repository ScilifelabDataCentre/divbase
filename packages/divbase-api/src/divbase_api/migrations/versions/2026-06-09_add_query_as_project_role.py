"""add_query_as_project_role

Revision ID: 9885301b9073
Revises: 427a022ebccf
Create Date: 2026-06-09 14:10:27.943160

This migration adds the new project role QUERY.

Enums in Postgres are a little special, hence the below steps.
"""

from typing import Sequence, Union

from alembic import op

# revision identifiers, used by Alembic.
revision: str = "9885301b9073"
down_revision: Union[str, Sequence[str], None] = "427a022ebccf"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Upgrade schema."""
    # manually created steps to add the QUERY role to the projectroles enum
    # Convert the colum to string first, so we can drop and then recreate the new enum type safely.
    op.execute("ALTER TABLE project_membership ALTER COLUMN role TYPE VARCHAR(50)")
    op.execute("DROP TYPE projectroles")
    op.execute("CREATE TYPE projectroles AS ENUM ('READ', 'QUERY', 'EDIT', 'MANAGE')")
    op.execute("ALTER TABLE project_membership ALTER COLUMN role TYPE projectroles USING role::projectroles")


def downgrade() -> None:
    """Downgrade schema."""
    # Like upgrade reversed, but we just reset any QUERY roles to READ first
    op.execute("UPDATE project_membership SET role = 'READ' WHERE role = 'QUERY'")
    op.execute("ALTER TABLE project_membership ALTER COLUMN role TYPE VARCHAR(50)")
    op.execute("DROP TYPE projectroles")
    op.execute("CREATE TYPE projectroles AS ENUM ('READ', 'EDIT', 'MANAGE')")
    op.execute("ALTER TABLE project_membership ALTER COLUMN role TYPE projectroles USING role::projectroles")
