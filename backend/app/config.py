"""
Environment-driven settings. Nothing environment-specific belongs inline
in application code — CORS origins, cache location, and (later) any
credentials all come from here, sourced from env vars with sane local
defaults so `uvicorn app.api.main:app` still works with zero setup.
"""
from __future__ import annotations

import os
from pathlib import Path


def _split_csv(value: str) -> list[str]:
    return [v.strip() for v in value.split(",") if v.strip()]


class Settings:
    def __init__(self) -> None:
        self.cors_origins: list[str] = _split_csv(
            os.environ.get("CORS_ORIGINS", "http://localhost:5173,http://localhost:3000")
        )
        # Cache root for downloaded datasets. In a container this should be
        # a mounted volume; locally it defaults to <repo>/data/cache.
        default_cache = Path(__file__).resolve().parents[2] / "data" / "cache"
        self.cache_dir: Path = Path(os.environ.get("CACHE_DIR", str(default_cache)))
        self.environment: str = os.environ.get("APP_ENV", "development")
        # Warm every dataset's cache on a background thread at startup.
        # On by default so a deployed container populates its volume
        # without a user waiting on the first request; off for tests and
        # any local run that shouldn't pull ~500MB. Never blocks the
        # server from binding either way -- see app/services/warmup.py.
        self.warmup_on_startup: bool = os.environ.get("WARMUP_ON_STARTUP", "1") not in (
            "0",
            "false",
            "False",
        )


settings = Settings()
