"""task_name_to_task_history

Revision ID: 1d9d677c1eee
Revises: 427a022ebccf
Create Date: 2026-05-22 00:00:00.000000

"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "1d9d677c1eee"
down_revision: Union[str, Sequence[str], None] = "427a022ebccf"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("task_history", sa.Column("task_name", sa.String(), nullable=True))
    op.create_index(op.f("ix_task_history_task_name"), "task_history", ["task_name"], unique=False)
    op.create_index(
        "ix_task_history_project_id_task_name",
        "task_history",
        ["project_id", "task_name"],
        unique=False,
        postgresql_where=sa.text("task_name IS NOT NULL"),
    )


def downgrade() -> None:
    """Downgrade schema."""
    op.drop_index("ix_task_history_project_id_task_name", table_name="task_history", if_exists=True)
    op.drop_index(op.f("ix_task_history_task_name"), table_name="task_history")
    op.drop_column("task_history", "task_name")
