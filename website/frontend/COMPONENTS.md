# frontend_sparql — Component Reference

React 19 + TypeScript SPA.  Vite build, served by nginx.
Nginx proxies `/api/` → `backend:8000`, so all API calls use relative paths (no hardcoded host).

## App structure

```
src/
  App.tsx              — Tab bar, routing between views, language switcher (EN/AR)
  App.css              — All global styles (BEM-ish classes, CSS variables)
  i18n/
    en.json            — English strings
    ar.json            — Arabic strings (RTL)
  components/
    PhenotypeSearch.tsx
    VariantLookup.tsx
    DiseaseBrowser.tsx
    GeneBrowser.tsx
    CaseDetail.tsx
    AboutPage.tsx
    SourceBadge.tsx
```

## Components

### PhenotypeSearch

Multi-select HPO autocomplete → POST `/api/search/phenotype`.

Options:
- **Method**: Lin+BMA (default) / Resnik+BMA
- **Expand with OMIM/ClinVar disease HPO annotations**: calls backend with `include_disease_phenotypes: true`
- **Include non-Saudi cases (DDD / Literature / ClinVar)**: `include_non_saudi: true`
- **Max results**: 10 / 20 / 50 / 100

Score display: raw decimal (e.g. `0.9286`), not a percentage.  Lin is nominally [0,1];
Resnik is unbounded IC units.  Never multiply by 100.

### CaseDetail

Loaded when user clicks a case ID link (`/case/{id}` → `GET /api/case/{id}`).

Key response fields and how they are rendered:

| Field | Rendering |
|---|---|
| `properties.progressStatus` | "Solved" / "Unsolved" / "Unknown" (mapped by backend) |
| `properties.isSaudi` | Saudi flag badge if `"true"` or `"1"` |
| `properties.sex` | Detected from NCIT URI substring (`C16576`=Female, `C20197`=Male) |
| `phenotypes` | `[{id, label}]` array → amber tag links to HPO browser: "Microcephaly (HP:0001263)" |
| `excluded_phenotypes` | Same format, styled with strikethrough/gray |
| `suggested_diseases` | Shown only when `diseaseLabel` absent; dashed amber box, italic, "⚠ Suggested (not a diagnosis)" note; links to OMIM |
| `variants[].togovar_url` | `/api/togovar-search?chrom=X&pos=Y` proxy path → browser follows redirect to TogoVar |

**Important**: `phenotypes` and `excluded_phenotypes` come from the backend as arrays, not
from `properties`.  The old code extracted from `properties.hasPhenotype` (which was
overwritten to a single value by the backend dict assignment bug — now fixed).

### VariantLookup

Four search modes via tabs: Gene / rsID / HGVS / ACMG class.  All call `GET /api/search/variant`.

### DiseaseBrowser

`GET /api/search/disease?q=` — searches `pavs:diseaseLabel` in all case graphs.

### GeneBrowser

`GET /api/gene/{symbol}` — fetches gene properties from `graph/genes`.

### AboutPage

Fetches `/about.md` (static file served by nginx) and renders via `react-markdown`.
Edit `public/about.md` without rebuilding.

### SourceBadge

Maps `pavs:source` values to coloured badges:
- `ahmed-variants` / other Saudi sources → green ("Saudi cohort")
- `ddd-diagnoses` → gray ("DDD")
- `Literature` → lavender ("Literature")
- `ClinVar` → blue ("ClinVar")

## i18n

All UI strings go through `useTranslation()`.  When adding a new string:
1. Add to `src/i18n/en.json`
2. Add equivalent to `src/i18n/ar.json`

Keys added in this session:
- `search.includeDisease` — "Expand with OMIM/ClinVar disease HPO annotations"
- `search.includeNonSaudi` — "Include non-Saudi cases (DDD / Literature / ClinVar)"
- `case.suggestedDiseaseNote` — "⚠ Suggested (not a diagnosis) — based on gene, not confirmed clinically:"
- `case.basedOnGene` — "via gene"

## Nginx proxy

`/api/` → `http://backend:8000` (production Docker)

Vite dev server: same proxy to `http://localhost:8000` (see `vite.config.ts`).

This means `togovar_url` values like `/api/togovar-search?chrom=3&pos=4417183` work
as `<a href>` links in both environments without any CORS or hardcoded-URL issues.

## CSS conventions

All custom classes in `App.css`.  Key classes for the suggested-disease feature:

```css
.suggested-disease-block   /* dashed amber border, yellow-tint background */
.suggested-label           /* italic, small, dark amber — the disclaimer line */
.suggested-disease-list    /* ul, no bullets */
.suggested-disease-list a  /* italic, bold, dark amber link to OMIM */
.suggested-id              /* gray, normal weight — "(OMIM:618792)" */
.suggested-gene            /* gray, normal weight — "via gene UGDH" */
```
