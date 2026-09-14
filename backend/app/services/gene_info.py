"""
Gene annotation lookup -- symbol/alias/Ensembl-ID/summary, via MyGene.info.

Free, public, no auth or registration required (verified directly against
the live API, not assumed from docs). This is the "cheap win" identified
against GEPIA3's Gene Card feature (see reference/gepia3-feature-review.md)
and doubles as the fix for METHODS.md 7.4 (gene-symbol aliases): looking a
symbol up here, rather than only checking a dataset's own exact-match
index, is what lets a user search "CIP2A" and find a gene a platform's own
annotation only knows as "KIAA1524".

Deliberately NOT wired into `_resolve_gene()` in api/main.py yet -- that
would change what genes a query can match, which needs its own decision
(silently falling back to an alias match changes analysis behavior, not
just the annotation panel). This module only powers the read-only
annotation lookup for now; alias-based gene *resolution* stays a separate
decision for whoever revisits METHODS.md 7.4 next.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import httpx

from app.config import settings

MYGENE_QUERY_URL = "https://mygene.info/v3/query"
CACHE_DIR = settings.cache_dir / "gene_info"
CACHE_TTL_SECONDS = 30 * 24 * 60 * 60  # 30 days -- gene annotations change rarely

FIELDS = "symbol,name,summary,alias,ensembl.gene,entrezgene"


def _cache_path(symbol: str) -> Path:
    # Gene symbols are a small, safe alphabet (letters/digits/hyphen) in
    # practice, but never trust external input as a filename verbatim.
    safe = "".join(c for c in symbol.upper() if c.isalnum() or c in "-_")
    return CACHE_DIR / f"{safe}.json"


def _read_cache(path: Path) -> dict | None:
    if not path.exists():
        return None
    if time.time() - path.stat().st_mtime > CACHE_TTL_SECONDS:
        return None
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def lookup_gene(symbol: str) -> dict | None:
    """Best-effort annotation for one gene symbol (or alias). None if
    MyGene.info has nothing under that name -- not every dataset's local
    feature_id/symbol will resolve externally, and that's not an error,
    just nothing to show in the annotation panel."""
    path = _cache_path(symbol)
    cached = _read_cache(path)
    if cached is not None:
        return cached or None  # cached {} means "looked up, found nothing"

    resp = httpx.get(
        MYGENE_QUERY_URL,
        params={"q": f"symbol:{symbol}", "species": "human", "fields": FIELDS, "size": 1},
        timeout=10.0,
    )
    resp.raise_for_status()
    hits = resp.json().get("hits", [])

    if not hits:
        # Fall back to a free-text query (not exact-symbol) -- this is
        # what lets an alias like KIAA1524 resolve to CIP2A (verified
        # directly: mygene.info/v3/query?q=KIAA1524 returns the CIP2A
        # record with KIAA1524 listed under `alias`).
        resp = httpx.get(
            MYGENE_QUERY_URL,
            params={"q": symbol, "species": "human", "fields": FIELDS, "size": 1},
            timeout=10.0,
        )
        resp.raise_for_status()
        hits = resp.json().get("hits", [])

    result = _shape(hits[0]) if hits else {}

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result))
    return result or None


def _ensembl_gene_id(hit: dict) -> str | None:
    """mygene.info returns `ensembl` as a dict for most genes, but as a
    LIST of dicts when a symbol maps to several Ensembl genes -- ZAP70 is
    one, which crashed this endpoint with
    `AttributeError: 'list' object has no attribute 'get'` (caught in QA
    by watching the server log, not the UI: the pane just showed no
    annotation, so nothing looked wrong from the browser).

    Where there are several, take the first stable `gene` id rather than
    inventing a merge -- callers want one id to link out with, and the
    alternatives here are transcript-level variants of the same locus.
    """
    ensembl = hit.get("ensembl")
    if isinstance(ensembl, list):
        for entry in ensembl:
            if isinstance(entry, dict) and entry.get("gene"):
                return entry["gene"]
        return None
    if isinstance(ensembl, dict):
        return ensembl.get("gene")
    return None


def _aliases(hit: dict) -> list[str]:
    """mygene.info returns `alias` as a list for genes with several
    aliases but as a BARE STRING for a gene with exactly one -- LYL1
    ("bHLHa18") is the case that mattered, since LYL1 ships in this app's
    own ETP-TF5 preset.

    The contract declares `aliases: list[str]`, and the frontend does
    `aliases.length > 0 && aliases.join(", ")`. A string passes `.length`
    (8, truthy) and then throws on `.join`, which white-screened the whole
    workspace rather than degrading one panel (caught in QA).

    Same polymorphism as `_ensembl_gene_id` handles one field over; this
    one was missed on the first pass, so normalise every shape here rather
    than trusting the upstream type.
    """
    alias = hit.get("alias")
    if alias is None:
        return []
    if isinstance(alias, str):
        return [alias]
    if isinstance(alias, list):
        return [str(a) for a in alias if a]
    return []


def _shape(hit: dict) -> dict:
    return {
        "symbol": hit.get("symbol"),
        "name": hit.get("name"),
        "summary": hit.get("summary"),
        "aliases": _aliases(hit),
        "ensembl_gene_id": _ensembl_gene_id(hit),
        "entrez_gene_id": hit.get("entrezgene"),
    }
