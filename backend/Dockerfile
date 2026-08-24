FROM python:3.12-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

# Downloaded dataset matrices live here — mount a volume in production so
# a container restart doesn't force a multi-hundred-MB re-download.
ENV CACHE_DIR=/data/cache
RUN mkdir -p /data/cache

EXPOSE 8420

# Railway (and most PaaS) inject $PORT at runtime rather than using a fixed
# port -- default to 8420 for local/docker-compose use where PORT isn't set.
CMD uvicorn app.api.main:app --host 0.0.0.0 --port ${PORT:-8420}
