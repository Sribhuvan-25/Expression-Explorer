"""
Shared network-fetch helper for ingest modules pulling a single file over
HTTP (a clinical supplement, a supplementary table) — as opposed to the
hundreds-of-files case (GDC MAFs, STAR-Counts), which retries inline at
the call site because each file also needs its own gunzip/parse step.

Exists because a real production failure showed the gap this closes:
`target_mrd.py`'s single-shot `httpx.get(...).content` fetch of TARGET's
clinical supplement returned a 200 with a body that failed pandas'
"cannot be determined" sniff — not a connection drop `httpx.TransportError`
would catch, and not a 4xx/5xx `raise_for_status()` would catch either.
The proximate cause (a corrupted/truncated body on an otherwise-200
response) couldn't be pinned down further after the fact, but the module
had zero defense against it: no retry, and no check that the bytes it
got back actually looked like the file it asked for before handing them
to pandas. `load_mrd_status()` failing killed the entire 469-sample
target_all_p2 dataset in production (2026-09-21) — a bare call with no
retry and no isolation, for a single supplementary grouping column,
brought down the dataset's expression matrix and survival support along
with it.

`gdc_target.py`'s own `_download_star_counts()` already retries
`httpx.TransportError` with backoff, with a comment noting GDC drops
connections from cloud IPs "observed repeatedly from Railway, not
reproducible locally" — the same platform, the same symptom. This
generalizes that pattern to the single-file case and adds the piece that
pattern didn't have: verifying the response actually looks like the
expected file type before parsing it, so a bad body fails fast with a
clear message on attempt 1 instead of reaching pandas with garbage on
whichever attempt happens to return one.
"""
from __future__ import annotations

import time

import httpx

# Retryable at the transport layer (connection drops, timeouts) AND at
# the application layer (a 200 whose body isn't the file promised) --
# both observed from this exact platform, so both get the same backoff
# rather than only the one `httpx.TransportError` already covers.
class BadResponseBody(Exception):
    pass


def fetch_with_retry(
    url: str,
    *,
    timeout: float = 60.0,
    attempts: int = 4,
    validate: "callable[[bytes], None] | None" = None,
) -> bytes:
    """GET `url`, retrying transient failures with exponential backoff.

    `validate(content)` should raise `BadResponseBody` if the bytes don't
    look like what was expected (e.g. missing a file-format magic
    number) — that failure is retried exactly like a transport error,
    since both are "the network gave us something unusable" and a retry
    is cheap relative to failing the whole dataset over it.

    A real HTTP error status (4xx/5xx) still raises immediately via
    `raise_for_status()` — that's not transient, retrying it just
    delays a failure that won't change.
    """
    last_exc: Exception | None = None
    for attempt in range(attempts):
        try:
            resp = httpx.get(url, follow_redirects=True, timeout=timeout)
            resp.raise_for_status()
            content = resp.content
            if validate is not None:
                validate(content)
            return content
        except (httpx.TransportError, BadResponseBody) as exc:
            last_exc = exc
            if attempt < attempts - 1:
                time.sleep(2**attempt)
    assert last_exc is not None
    raise last_exc


# Excel 2007+ (.xlsx) is a zip archive; the legacy .xls binary format has
# its own distinct magic number. Checking both because GDC's clinical
# supplements and Liu et al.'s supplementary table use modern .xlsx, but
# validating on the narrower signature would be a trap for the next file
# added this way.
_XLSX_MAGIC = b"PK\x03\x04"
_XLS_MAGIC = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"


def validate_excel_bytes(content: bytes) -> None:
    if not content:
        raise BadResponseBody("response body is empty")
    if not (content.startswith(_XLSX_MAGIC) or content.startswith(_XLS_MAGIC)):
        raise BadResponseBody(
            f"response body ({len(content)} bytes) doesn't start with an Excel "
            f"file signature -- got {content[:16]!r}"
        )
