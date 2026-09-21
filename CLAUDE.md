# Expression Explorer — project rules

Gene expression analysis platform: group comparison, signature scoring,
survival analysis. FastAPI backend, React/Vite frontend, datasets behind a
shared contract.

> **Note (2026-09-20):** this file was accidentally deleted and
> reconstructed from the project's conventions. It was never committed, so
> no original exists to diff against. The rules below are accurate as far
> as they go, but if you remember a rule that used to be here and isn't,
> it was lost — please add it back.

## Non-negotiables

**METHODS.md is the record.** Every analytical decision that can change a
number a user reads gets an entry: what was decided, why, and where it
lives in the code. Update the entry in place when a decision changes and
add a line to the change log — never silently rewrite history. This is
the project's most important convention: the tool exists so researchers
can analyse without a biostatistician, which only works if every choice
is stated and defensible.

**Sample accounting is mandatory.** Any analysis running on fewer samples
than the dataset holds must report `n_dataset_total`, `n_excluded`, and a
plain-language reason, surfaced in the UI. Where samples are excluded for
genuinely different reasons, report them separately — a merged count
misattributes samples and is the failure mode this rule exists to
prevent. See METHODS.md §6.1, §8.6.

**Never pool values across datasets.** Different assay types and units
are not on a comparable scale. Results are shown side by side, each with
its own axis and statistics. See METHODS.md §7.2.

**Avoid circular analysis.** Probe selection, thresholds, and gene
choices must be independent of the group labels being compared.

## Adding a dataset

Write `backend/app/ingest/<name>.py` with a `load()` returning a
`Dataset`, then run `python scripts/ingest_dataset.py <name>`. The
registry and API pick it up from the metadata database — no API code
changes, no redeploy. See METHODS.md §10.

Auxiliary measurements on an existing cohort's samples (mutations, CRISPR,
drug response) are `AuxLayer`s on that dataset, not new datasets — same
samples, a second measurement. See METHODS.md §8.1.

## Backend

- Python 3.12, FastAPI. `app/api/main.py` is thin by design: it resolves a
  dataset through the registry and calls into `app/analysis/`. No
  dataset- or disease-specific logic belongs in it.
- Route ordering matters: static paths must be registered **before**
  dynamic `{param}` paths or they get shadowed.
- Normalise user input centrally (`_resolve_gene`, `_resolve_aux_feature`)
  rather than per-endpoint — inconsistent normalisation has bitten twice.
- Prefer vectorised scipy/pandas over per-gene Python loops; verify a
  rewrite is bit-identical before trusting the speedup.
- Dataset warm-up must never block the server binding or `/datasets`.
  See METHODS.md §7.9.

## Frontend

- React 19 + TypeScript + Vite, dockview panes. Node **22.12.0**
  (`nvm use 22.12.0`) — Node 18 fails the build.
- Keep both tabs mounted and toggle with CSS rather than conditionally
  rendering; a conditional render unmounts and destroys child state.
- Never let one bad API field take down the page — components that render
  inside every results page have no error boundary above them.

## Ports

Backend **8420** (`vite.config.ts` proxies `/api` there, and
`docker-compose.yml`/`Dockerfile` agree). Frontend dev server 5173.
Starting the backend anywhere else makes every frontend call 502 and every
pane render empty.

## Testing

- `pytest` in `backend/`, `npm test` in `frontend/`.
- Tests are for defects that actually shipped, not for a coverage number.
- **Validate a new test by reintroducing the bug it guards and confirming
  it fails.** A test that cannot fail is worse than no test, because it
  reads as coverage. This has caught a useless test at least once.
- Don't let tests hit the network — pass `with_aux=False` or point at
  fixtures.

## Verification

Claims about behaviour get checked against the running system, not
asserted from memory. Numbers in commit messages and METHODS.md should be
measured. Where something is unverified, say so explicitly rather than
implying it was tested.
