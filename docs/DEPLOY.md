# PAVS Deployment Guide

## Overview

The PAVS stack runs as five Docker services orchestrated by `docker-compose-sparql.yml`:

| Service | Image | Role |
|---|---|---|
| `virtuoso` | `openlink/virtuoso-opensource-7` | Virtuoso triple store (SPARQL database) |
| `init` | `Dockerfile.pavs` | Generates RDF from source data (runs once) |
| `loader` | `openlink/virtuoso-opensource-7` | Bulk-loads TTL files into Virtuoso via isql |
| `backend` | `Dockerfile.pavs` | FastAPI SPARQL search API |
| `frontend` | `website/frontend/Dockerfile` | React app served by nginx |

All source data (annotated TSV, HPO ontology, phenotype.hpoa, Arabic translations, literature phenopackets, pre-generated phenopackets) is **baked into the `Dockerfile.pavs` image** — no host data files are required for normal operation.

Service start order: `virtuoso` → `init` → `loader` → `backend` → `frontend`

---

## Quick Start (local)

```bash
# 1. Copy and review the environment file
cp .env.example .env

# 2. Start the full stack
docker compose -f docker-compose-sparql.yml up -d

# Backend startup takes ~15 s (loads IC values and caches)
# Health check:
curl http://localhost:8000/api/health

# Frontend:  http://localhost:3000
# API docs:  http://localhost:8000/docs
# Virtuoso:  http://localhost:8890/sparql
```

---

## Configuration (`.env`)

Copy `.env.example` to `.env` and adjust as needed. Docker Compose automatically loads `.env` from the same directory.

| Variable | Default | Description |
|---|---|---|
| `PAVS_PUBLIC_URL` | *(empty)* | Public base URL (scheme + host, no trailing slash). Leave blank to use relative paths — the nginx frontend proxies `/api/` to the backend automatically. Set to `https://pavs.phenomebrowser.net` for production. |
| `PAVS_FRONTEND_PORT` | `3000` | Host port for the nginx frontend (map to `80` in production) |
| `PAVS_BACKEND_PORT` | `8000` | Host port for the FastAPI backend |
| `PAVS_VIRTUOSO_HTTP_PORT` | `8890` | Host port for Virtuoso HTTP / SPARQL UI |
| `PAVS_VIRTUOSO_SQL_PORT` | `1111` | Host port for Virtuoso iSQL |
| `VIRT_DBA_PASSWORD` | `pavs_dba` | Virtuoso DBA password — **change before production** |

---

## Production Deployment (`pavs.phenomebrowser.net`)

### Server

The production server is **`cbontsr01.kaust.edu.sa`**, configured as SSH host `onto`:

```
Host onto
    Hostname cbontsr01.kaust.edu.sa
    User hohndor
```

The PAVS stack lives at `/data/pavs/` on that server (not a git repo — deploy by copying files).

**Port mapping** (Docker → host):
| Service | Container port | Host port |
|---|---|---|
| frontend (nginx) | 80 | 20000 |
| backend (FastAPI) | 8000 | 20001 |
| Virtuoso HTTP | 8890 | (internal) |
| Virtuoso iSQL | 1111 | (internal) |

All external traffic goes through **Imperva CDN** (`X-CDN: Imperva`) in front of port 80.

**CDN note**: Imperva aggressively caches static files. Dynamic content must be served through
`/api/` to bypass the cache. For this reason, `about.md` is served via `GET /api/about?lang=en|ar`
(backend reads `website/backend/about.md` / `about_ar.md`) rather than as a static file.

### Deploying code changes to production

```bash
# 1. Copy changed files to onto
scp website/backend/main.py onto:/data/pavs/website/backend/main.py
scp website/frontend/src/components/SomeComponent.tsx onto:/data/pavs/website/frontend/src/components/SomeComponent.tsx
# ... etc.

# 2. Rebuild and restart on onto
ssh onto 'cd /data/pavs && docker compose -f docker-compose-sparql.yml build backend frontend && docker compose -f docker-compose-sparql.yml up -d backend frontend'

# 3. Verify
curl https://pavs.phenomebrowser.net/api/health
```

### Updating about.md / about_ar.md

These files live in `website/backend/` (baked into the backend image) and are served via `/api/about`:

```bash
scp website/frontend/public/about.md onto:/data/pavs/website/backend/about.md
scp website/frontend/public/about_ar.md onto:/data/pavs/website/backend/about_ar.md
# Hot-update running container (no restart needed):
ssh onto 'docker cp /data/pavs/website/backend/about.md pavs-backend-1:/app/backend/about.md && docker cp /data/pavs/website/backend/about_ar.md pavs-backend-1:/app/backend/about_ar.md'
```

### Minimal `.env`

```bash
PAVS_PUBLIC_URL=https://pavs.phenomebrowser.net
PAVS_FRONTEND_PORT=80
PAVS_BACKEND_PORT=8000     # can be firewalled; not required externally
PAVS_VIRTUOSO_HTTP_PORT=8890  # can be firewalled; Virtuoso is proxied via /sparql
PAVS_VIRTUOSO_SQL_PORT=1111
VIRT_DBA_PASSWORD=<strong-password>
```

### Start

```bash
docker compose -f docker-compose-sparql.yml up -d
```

### Public endpoints (all served through nginx on port 80)

| URL | Description |
|---|---|
| `https://pavs.phenomebrowser.net/` | React frontend |
| `https://pavs.phenomebrowser.net/api/` | FastAPI backend (proxied by nginx) |
| `https://pavs.phenomebrowser.net/api/health` | Health check |
| `https://pavs.phenomebrowser.net/api/sparql` | PAVS SPARQL proxy (SELECT only, JSON/TSV/CSV) |
| `https://pavs.phenomebrowser.net/sparql` | Raw Virtuoso SPARQL endpoint (full protocol) |
| `https://pavs.phenomebrowser.net/api/docs` | FastAPI OpenAPI docs |

### Programmatic SPARQL access

```bash
# JSON (default)
curl 'https://pavs.phenomebrowser.net/api/sparql?format=json&query=SELECT+...'

# TSV download
curl 'https://pavs.phenomebrowser.net/api/sparql?format=tsv&query=SELECT+...' > results.tsv

# CSV download
curl 'https://pavs.phenomebrowser.net/api/sparql?format=csv&query=SELECT+...' > results.csv

# Raw Virtuoso endpoint (GET, full SPARQL protocol)
curl 'https://pavs.phenomebrowser.net/sparql?query=SELECT+...&format=application%2Fsparql-results%2Bjson'
```

---

## After Code Changes (local or production)

```bash
# Backend or frontend code changed:
docker compose -f docker-compose-sparql.yml build backend frontend
docker compose -f docker-compose-sparql.yml up -d backend frontend

# Data pipeline changed (intake/, normalization/):
docker compose -f docker-compose-sparql.yml build init backend
docker compose -f docker-compose-sparql.yml up -d

# Force full RDF regeneration (e.g. after changes to compute_hpo_ic.py):
docker compose -f docker-compose-sparql.yml run --rm init python intake/prepare_all.py --endpoint=http://virtuoso:8890 --skip-load --force
docker compose -f docker-compose-sparql.yml up loader
```

---

## Rebuilding the Knowledge Graph

The `init` service generates RDF from baked-in source data and writes it to the `pavs_rdf` named volume. It skips steps whose outputs already exist in the volume.

```bash
# Force full rebuild (regenerates all RDF and reloads into Virtuoso)
docker compose -f docker-compose-sparql.yml run --rm init --force

# Then restart the loader to reload into Virtuoso
docker compose -f docker-compose-sparql.yml up loader
```

To rebuild and start fresh (drops existing Virtuoso data):

```bash
docker compose -f docker-compose-sparql.yml down
docker volume rm pavs_pavs_rdf
docker compose -f docker-compose-sparql.yml up -d
```

---

## Volumes

| Volume | Contents | Shared between |
|---|---|---|
| `virtuoso_db` | Virtuoso database files (persists across restarts) | `virtuoso` |
| `pavs_rdf` | Generated RDF/Turtle files (`*.ttl`) | `init` (writes), `loader` + `virtuoso` (read) |

Source data and pre-generated phenopackets are baked into the `Dockerfile.pavs` image — they are **not** stored in volumes.

---

## Running Tests

```bash
# Requires the backend to be running
uv run pytest tests/test_sparql_queries.py -v

# Against a remote deployment
PAVS_API=https://pavs.phenomebrowser.net uv run pytest tests/test_sparql_queries.py -v
```

---

## Troubleshooting

### `init` fails with "No such file or directory"
The `init` service is using a stale image. Rebuild it:
```bash
docker compose -f docker-compose-sparql.yml build init
docker compose -f docker-compose-sparql.yml up -d
```

### Backend starts but returns empty results
The loader may not have completed. Check:
```bash
docker logs pavs-loader-1
docker logs pavs-backend-1
```

### Virtuoso not accepting connections
The `virtuoso` healthcheck retries for up to 2.5 minutes. On slow machines, increase `retries` or `start_period` in `docker-compose-sparql.yml`.

### VITE_API_URL and the API endpoint
- `PAVS_PUBLIC_URL` is passed as `VITE_API_URL` at **build time** to the frontend.
- Leave it blank (default) — the frontend uses relative paths and nginx proxies `/api/` to the backend. This works for any hostname without rebuilding.
- Set it to the full URL only if the API is on a **different origin** than the frontend.
