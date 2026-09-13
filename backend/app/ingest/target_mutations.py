"""
TARGET-ALL-P2 somatic mutation calls -- an auxiliary layer on the existing
`target_all_p2` dataset.

Open-access Masked Somatic Mutation MAFs from GDC (739 files, one per
aliquot; verified open-access via the GDC API, same status as the RNA-seq
this project already ingests). Attached as an `AuxLayer` on the existing
patient cohort rather than registered as its own dataset -- same patients,
a second measurement -- see METHODS.md 8.1.

**Partial coverage is the headline caveat.** 717 distinct cases have MAF
files, but only **279 of the 469 samples in our RNA-seq matrix** overlap
(the rest were exome-sequenced without matching RNA-seq, or vice versa).
Any analysis stratifying expression by mutation status therefore runs on
~60% of the cohort and MUST say so -- the exclusion-transparency pattern
in METHODS.md 6.1 exists for exactly this, and skipping it here would
repeat the n=466 confusion that the domain-expert review already caught
once.

Gene-level only for now: a sample is "mutated" in gene X if it carries at
least one non-silent variant in X. Amino-acid-level resolution (GEPIA3's
other mode) needs per-call HGVSp positions retained; this loader keeps the
columns that would allow it but doesn't build that index yet -- see
METHODS.md 8.2.
"""
from __future__ import annotations

import gzip
import io
import logging
import time

import httpx
import numpy as np
import pandas as pd

from app.config import settings
from app.models.contract import AuxLayer, AuxMatrix

GDC_API = "https://api.gdc.cancer.gov"
PROJECT_ID = "TARGET-ALL-P2"
CACHE_DIR = settings.cache_dir / "target_mutations"

log = logging.getLogger(__name__)

# Variant classes that change the protein product. Everything outside this
# set (Silent, 3'UTR, 5'UTR, Intron, RNA, IGR, ...) is excluded from the
# boolean mutation matrix: counting a synonymous or untranslated-region
# variant as "this gene is mutated" would inflate every mutation rate with
# calls that don't alter the protein. Verified against a real file: 7 of
# 21 calls in one sampled aliquot were Silent and 2 were 3'UTR, so this is
# not a hypothetical distinction.
#
# This is the standard non-silent set used in MAF analysis; it's written
# out explicitly rather than expressed as "not Silent" so the decision is
# auditable and a future reader can see exactly what was counted.
NON_SILENT = frozenset(
    {
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Splice_Site",
        "Translation_Start_Site",
    }
)


def _list_maf_files(client: httpx.Client) -> list[dict]:
    payload = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "in", "content": {"field": "cases.project.project_id", "value": [PROJECT_ID]}},
                {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
                {"op": "in", "content": {"field": "access", "value": ["open"]}},
            ],
        },
        "fields": "file_id,cases.submitter_id",
        "size": 2000,
    }
    resp = client.post("/files", json=payload)
    resp.raise_for_status()
    return resp.json()["data"]["hits"]


def _download_maf(client: httpx.Client, file_id: str) -> pd.DataFrame:
    """One MAF -> (Hugo_Symbol, Variant_Classification) rows.

    Retries transient transport failures the same way gdc_target.py does,
    and for the same reason: GDC drops connections from cloud IPs often
    enough that a single-attempt request kills a 700+ file run over one
    flaky file.
    """
    last_exc: Exception | None = None
    for attempt in range(4):
        try:
            resp = client.get(f"/data/{file_id}", follow_redirects=True)
            resp.raise_for_status()
            raw = resp.content
            break
        except httpx.TransportError as exc:
            last_exc = exc
            if attempt < 3:
                time.sleep(2**attempt)
    else:
        raise last_exc

    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    lines = [ln for ln in raw.decode(errors="replace").splitlines() if not ln.startswith("#")]
    if not lines:
        return pd.DataFrame(columns=["Hugo_Symbol", "Variant_Classification"])
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", low_memory=False)
    keep = [c for c in ("Hugo_Symbol", "Variant_Classification") if c in df.columns]
    return df[keep] if len(keep) == 2 else pd.DataFrame(columns=["Hugo_Symbol", "Variant_Classification"])


def load_mutations(use_cache: bool = True, limit: int | None = None) -> AuxMatrix | None:
    """Boolean mutation matrix (gene symbol x sample_id): True if that
    sample carries >= 1 non-silent variant in that gene.

    Returns None on any failure to reach GDC -- additive layer, must not
    take down the expression dataset it attaches to.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path = CACHE_DIR / "mutation_status.parquet"

    if use_cache and cache_path.exists():
        matrix = pd.read_parquet(cache_path)
    else:
        with httpx.Client(base_url=GDC_API, timeout=60.0) as client:
            files = _list_maf_files(client)
            if limit:
                files = files[:limit]

            # sample_id -> set of mutated gene symbols. Accumulating sets
            # keyed by sample keeps peak memory to the distinct symbols
            # actually seen, rather than materialising a 20k x 700 dense
            # boolean frame before it's needed.
            by_sample: dict[str, set[str]] = {}
            for f in files:
                case = (f.get("cases") or [{}])[0]
                sample_id = case.get("submitter_id")
                if not sample_id:
                    continue
                try:
                    calls = _download_maf(client, f["file_id"])
                except Exception as exc:  # noqa: BLE001 - one bad file shouldn't kill the run
                    log.warning("target_mutations: skipping %s (%s)", f["file_id"], exc)
                    continue
                non_silent = calls[calls["Variant_Classification"].isin(NON_SILENT)]
                mutated = set(non_silent["Hugo_Symbol"].dropna().astype(str))
                by_sample.setdefault(sample_id, set()).update(mutated)

        if not by_sample:
            return None

        all_genes = sorted({g for genes in by_sample.values() for g in genes})
        sample_ids = sorted(by_sample)
        # Build the boolean frame directly rather than assembling sparse
        # True-only Series and filling the gaps: `.fillna(False)` on an
        # object-dtype frame is deprecated in pandas (silent downcast to
        # bool, slated to change behaviour), and constructing the dense
        # array up front is both warning-free and unambiguous about dtype.
        gene_index = {g: i for i, g in enumerate(all_genes)}
        data = np.zeros((len(all_genes), len(sample_ids)), dtype=bool)
        for col, sid in enumerate(sample_ids):
            for gene in by_sample[sid]:
                data[gene_index[gene], col] = True
        matrix = pd.DataFrame(data, index=all_genes, columns=sample_ids)
        matrix.to_parquet(cache_path)

    return AuxMatrix(
        layer=AuxLayer.MUTATION_STATUS,
        matrix=matrix,
        value_label="Non-silent mutation present",
        value_description=(
            "True if the sample carries at least one non-silent somatic variant in that "
            "gene (missense, nonsense, frameshift, splice-site, in-frame indel, or "
            "translation-start). Silent and untranslated-region variants are excluded."
        ),
        source_note=(
            f"GDC open-access Masked Somatic Mutation MAFs for {PROJECT_ID} "
            "(whole-exome; covers a subset of the RNA-seq cohort)."
        ),
    )
