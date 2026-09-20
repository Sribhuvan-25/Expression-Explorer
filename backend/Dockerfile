FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app
COPY scripts ./scripts

# Downloaded dataset matrices live here — mount a volume in production so
# a container restart doesn't force a multi-hundred-MB re-download.
ENV CACHE_DIR=/data/cache
RUN mkdir -p /data/cache

EXPOSE 8420

# Railway (and most PaaS) inject $PORT at runtime rather than using a fixed
# port -- default to 8420 for local/docker-compose use where PORT isn't set.
#
# uvicorn starts IMMEDIATELY; the cache warms on a background thread that
# the app itself launches (app/services/warmup.py). This used to be
# `python scripts/prewarm_cache.py && uvicorn ...`, which made the
# healthcheck wait on the download and turned deploy success into a race
# against GDC: the aux layers pushed cold start to ~848s against a 900s
# healthcheck budget, from a home connection, with the MAF loader
# retrying up to 4x per flaky file. Binding first removes the race
# entirely rather than buying time against it -- a bigger dataset now
# costs a longer warm-up, not a failed deploy.
CMD uvicorn app.api.main:app --host 0.0.0.0 --port ${PORT:-8420}
