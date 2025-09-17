"""
User DB Model.
"""

from sqlalchemy.orm import Mapped, mapped_column

from divbase_api.models.base import BaseDBModel


class UserDB(BaseDBModel):
    """
    DB Model for a user.

    id, created_at and updated_at are inherited from BaseDBModel.
    """

    __tablename__ = "user"

    name: Mapped[str] = mapped_column(index=True, unique=True)
    email: Mapped[str] = mapped_column(index=True, unique=True)
    hashed_password: Mapped[str]

    is_admin: Mapped[bool] = mapped_column(default=False)
    is_active: Mapped[bool] = mapped_column(default=True)
