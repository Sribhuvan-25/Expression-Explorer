# Expression Explorer

A gene expression analysis platform covering group comparison, signature
scoring, and survival analysis. Datasets plug in through a shared
contract, so new data sources and disease areas extend the tool without
rewriting the analysis layer.

**[METHODS.md](METHODS.md)** records every analytical decision — which
statistical tests are used and why, how signatures are scored, how
survival groups are split, which samples get excluded from an analysis,
and what is still open. Any change that affects a reported number should
be recorded there.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r backend/requirements.txt
```

## Run the API

```bash
cd backend
uvicorn app.api.main:app --reload --port 8420
```
