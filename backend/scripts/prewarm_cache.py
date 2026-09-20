"""
Pre-warms every registered dataset's on-disk cache, synchronously.

**No longer run at container startup.** It used to be the Dockerfile's
CMD prefix (`prewarm_cache.py && uvicorn ...`), which meant the platform
healthcheck could not pass until every download finished -- and once the
auxiliary layers landed, measured cold start reached ~848s against a
900s healthcheck budget. The app now warms itself on a background thread
while serving (app/services/warmup.py), so the deploy no longer races the
download.

This script is kept for the cases where blocking is what you actually
want: seeding a fresh volume on purpose, priming a dev machine before
going offline, or checking how long a cold load really takes. It shares
its per-dataset failure behaviour with the background warmer -- one
unreachable source is logged and skipped, not fatal.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app.registry import ensure_loaded, list_descriptors


def main() -> None:
    ensure_loaded()
    for descriptor in list_descriptors():
        start = time.monotonic()
        print(f"[prewarm] {descriptor.dataset_id}: loading...", flush=True)
        try:
            dataset = descriptor.loader()
        except Exception as exc:
            # A prewarm failure (e.g. GDC transiently unreachable at
            # container start) must not block the API from starting --
            # the lazy load-on-first-request path still exists as a
            # fallback, just without the timeout protection this script
            # normally provides. Log loudly and move on to the next
            # dataset rather than crashing container startup entirely.
            print(f"[prewarm] {descriptor.dataset_id}: FAILED ({exc}) -- will lazy-load on first request", flush=True)
            continue
        elapsed = time.monotonic() - start
        print(f"[prewarm] {descriptor.dataset_id}: ready ({dataset.matrix.shape[1]} samples, {elapsed:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
