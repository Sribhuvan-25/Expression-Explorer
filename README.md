# Expression Explorer

A gene expression analysis platform covering group comparison, signature
scoring, and survival analysis. Datasets plug in through a shared
contract, so new data sources and disease areas extend the tool without
rewriting the analysis layer.

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
