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
# Cache is pre-warmed before uvicorn binds, not lazily on first request:
# target_all_p2's first-ever load is a 530-file sequential GDC download,
# several minutes long, which exceeds Railway's ~5 minute reverse-proxy
# timeout if a live user request has to wait on it. Running the same
# loader here means the healthcheck (and therefore "deployment successful")
# doesn't fire until the data is actually ready, and every real request
# after that just reads the cached parquet file. Already-cached datasets
# (the common case after the first deploy, since the volume persists)
# make this a fast no-op.
CMD python scripts/prewarm_cache.py && uvicorn app.api.main:app --host 0.0.0.0 --port ${PORT:-8420}
