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


settings = Settings()
