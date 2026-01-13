"""Make task_id nullable in task_history

Revision ID: 3e168ddc857e
Revises: 5b559e3008f4
Create Date: 2026-01-13 09:12:26.855801

"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "3e168ddc857e"
down_revision: Union[str, Sequence[str], None] = "5b559e3008f4"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Upgrade schema."""
    op.drop_constraint("task_history_pkey", "task_history", type_="primary")
    op.create_primary_key("task_history_pkey", "task_history", ["id"])
    op.alter_column("task_history", "task_id", existing_type=sa.VARCHAR(), nullable=True)


def downgrade() -> None:
    """Downgrade schema."""
    op.alter_column("task_history", "task_id", existing_type=sa.VARCHAR(), nullable=False)
    op.drop_constraint("task_history_pkey", "task_history", type_="primary")
    op.create_primary_key("task_history_pkey", "task_history", ["task_id", "id"])
