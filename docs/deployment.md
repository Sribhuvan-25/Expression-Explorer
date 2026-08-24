# Deployment

Backend on Railway, frontend on Vercel, CI via GitHub Actions. Both platforms deploy
automatically on push to `main` once connected — no deploy step lives in CI itself.

## 1. Backend — Railway

1. **New Project → Deploy from GitHub repo** → select `Expression-Explorer`.
2. **Settings → Source → Root Directory**: set to `backend`. Railway's config file
   (`backend/railway.toml`) doesn't follow this setting itself, but the Dockerfile's
   relative `COPY` paths need it to resolve correctly.
3. **Settings → Volumes → New Volume**: mount path `/data/cache`. This is the fix for
   the real risk here — without a persistent volume, every redeploy/restart wipes the
   downloaded dataset cache and the next request re-triggers a ~10 minute sequential
   download from GDC (TARGET, ~530 files) or a large fetch from Figshare (DepMap).
4. **Variables**, set:
   - `CORS_ORIGINS` — the Vercel frontend's URL once deployed (step 2 below), e.g.
     `https://expression-explorer.vercel.app`. Comma-separate if there's more than one
     (e.g. a preview URL too).
   - `APP_ENV=production`
   - Railway sets `PORT` itself — don't override it; the Dockerfile already reads it.
5. Deploy. Note the public URL Railway assigns (Settings → Networking → Generate
   Domain if one isn't automatic).
6. **Pre-warm the cache** once the service is live, before sending the URL to anyone:
   ```bash
   curl https://<railway-url>/datasets/target_all_p2
   curl https://<railway-url>/datasets/depmap
   ```
   Both will hang for several minutes on this first call (that's the expected
   download) — let them finish. Every request after that is fast, cache-backed.

## 2. Frontend — Vercel

1. **Add New Project** → import the same GitHub repo.
2. Vercel should pick up `vercel.json` at the repo root automatically (build command,
   output directory already declared there). If it doesn't, set manually:
   - Build Command: `cd frontend && npm install && npm run build`
   - Output Directory: `frontend/dist`
3. **Environment Variables**, add:
   - `VITE_API_BASE_URL` = the Railway backend's public URL from step 1.5 above
     (no trailing slash, e.g. `https://expression-explorer-production.up.railway.app`).
4. Deploy. Vercel assigns a `*.vercel.app` URL.
5. **Go back to Railway** and set `CORS_ORIGINS` to this exact Vercel URL (step 1.4
   above) — the backend rejects cross-origin requests from anywhere not on that list.

## 3. CI — GitHub Actions

`.github/workflows/ci.yml` runs on every push/PR to `main`:
- Backend: installs `requirements.txt`, runs `pytest`.
- Frontend: installs deps, runs `tsc --noEmit`, runs `npm run build` (the real check —
  this caught a type error `tsc --noEmit` alone missed, see commit history).

This is a gate, not a deploy trigger — Railway and Vercel each watch the repo
independently via their own GitHub integration and redeploy on push once connected.

## Known limitations, honestly

- **First real query per dataset is slow** (see step 1.6) — this is inherent to the
  data sources (GDC's API serves one file per request, DepMap's file is just large),
  not something fixable by more caching. It only happens once per cache lifetime.
- **AALL0434-dependent results aren't available anywhere** — that cohort was never
  publicly obtainable (confirmed during this project, see `docs/plan.html`). Anything
  needing it will 400 with a clear message, not silently return wrong data.
- **The `data/cache` volume is not backed up.** If Railway's volume is ever deleted,
  the next request re-downloads from the original public sources — annoying, not
  catastrophic (no data is lost that can't be re-fetched).
