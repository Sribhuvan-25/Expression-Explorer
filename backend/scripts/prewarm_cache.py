"""
Pre-warms every registered dataset's on-disk cache before the API starts
serving traffic. Run once at container startup (see Dockerfile).

Without this, the first real user request to a not-yet-cached dataset
triggers its full ingest pipeline inline -- for target_all_p2 that's a
530-file sequential download from GDC, several minutes long, which
exceeds most PaaS reverse-proxy request timeouts (measured: Railway
returns 502 at ~5 minutes) even once the loader itself is memory-safe.
Running the same loaders here, before uvicorn binds, means every real
request just reads an already-assembled parquet file.
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
