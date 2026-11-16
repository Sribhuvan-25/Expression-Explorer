"""contrasts table (pre-specified design contrasts, METHODS.md 11.1)

Revision ID: b3c1f0a2d4e5
Revises: 97262eb77552
Create Date: 2026-11-15 18:00:00

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa

revision: str = "b3c1f0a2d4e5"
down_revision: Union[str, Sequence[str], None] = "97262eb77552"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "contrasts",
        sa.Column("id", sa.Integer(), autoincrement=True, nullable=False),
        sa.Column("dataset_id", sa.String(length=64), nullable=False),
        sa.Column("contrast_id", sa.String(length=64), nullable=False),
        sa.Column("spec", sa.JSON(), nullable=False),
        sa.Column("reference_samples", sa.JSON(), nullable=False),
        sa.Column("treatment_samples", sa.JSON(), nullable=False),
        sa.Column("details", sa.JSON(), nullable=False),
        sa.Column("results_uri", sa.String(length=1024), nullable=False),
        sa.ForeignKeyConstraint(["dataset_id"], ["datasets.dataset_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint("dataset_id", "contrast_id", name="uq_contrast_per_dataset"),
    )
    op.create_index("ix_contrasts_dataset", "contrasts", ["dataset_id"], unique=False)


def downgrade() -> None:
    op.drop_index("ix_contrasts_dataset", table_name="contrasts")
    op.drop_table("contrasts")
