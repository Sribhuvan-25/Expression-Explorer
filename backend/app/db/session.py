"""
Database engine and session handling.

`DATABASE_URL` selects the backend. Postgres in production (Railway's
add-on injects this variable automatically); SQLite for local runs and
tests, which keeps `pytest` and `uvicorn` working with zero setup --
the same zero-setup property app/config.py already preserves for
CACHE_DIR and CORS_ORIGINS.

The schema is written to work on both. That is why `models.py` uses the
generic `JSON` type rather than Postgres-native `JSONB`: SQLAlchemy maps
it to JSONB on Postgres and TEXT-backed JSON on SQLite, so tests exercise
the same models the deployment uses. If JSONB-specific indexing is ever
needed, it becomes a Postgres-only migration rather than a schema fork.
"""
from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session, sessionmaker

from app.config import settings
from app.db.models import Base


def _database_url() -> str:
    url = os.environ.get("DATABASE_URL", "").strip()
    if url:
        # Railway (and Heroku-style providers) hand out "postgres://",
        # which SQLAlchemy 2 does not accept -- it wants an explicit
        # driver. Rewriting here rather than asking every deployment to
        # edit the variable by hand.
        if url.startswith("postgres://"):
            url = url.replace("postgres://", "postgresql+psycopg://", 1)
        elif url.startswith("postgresql://"):
            url = url.replace("postgresql://", "postgresql+psycopg://", 1)
        return url
    # Local default: a file beside the dataset cache, so a dev machine
    # needs no database service running.
    settings.cache_dir.mkdir(parents=True, exist_ok=True)
    return f"sqlite:///{Path(settings.cache_dir) / 'metadata.db'}"


@lru_cache(maxsize=1)
def get_engine() -> Engine:
    url = _database_url()
    kwargs: dict = {"future": True}
    if url.startswith("sqlite"):
        # The warm-up thread and request handlers touch the same SQLite
        # connection pool; SQLite refuses cross-thread use by default.
        kwargs["connect_args"] = {"check_same_thread": False}
    else:
        # A dropped connection otherwise surfaces as a request-time error
        # long after the network blip that caused it.
        kwargs["pool_pre_ping"] = True
    return create_engine(url, **kwargs)


@lru_cache(maxsize=1)
def _session_factory() -> sessionmaker[Session]:
    return sessionmaker(bind=get_engine(), expire_on_commit=False, future=True)


@contextmanager
def session_scope() -> Iterator[Session]:
    """Transactional scope -- commits on success, rolls back on error."""
    session = _session_factory()()
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def create_all() -> None:
    """Create any missing tables.

    Alembic owns schema *changes*; this exists so a fresh local database
    or a test run works without first running a migration.
    """
    Base.metadata.create_all(get_engine())


def reset_engine() -> None:
    """Drop the cached engine/session factory. For tests that point
    DATABASE_URL somewhere new mid-process."""
    get_engine.cache_clear()
    _session_factory.cache_clear()
