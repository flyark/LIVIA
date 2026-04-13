# LIVIA Development

## Wiki Reference

The LIVIA development wiki is at:
`/Users/ah-ramkim/Library/CloudStorage/Dropbox/LLM_wiki/livia-dev/wiki/`

**Before starting any implementation task, check these wiki pages first:**

- `lessons/mistakes.md` — 21+ documented mistakes to avoid (PAE averaging, CORS, const re-declaration, etc.)
- `lessons/review-checklist.md` — Pre-push verification checklist
- `architecture/state-objects.md` — uniState, dimerPairResult, flyParsed structures
- `patterns/canvas-charts.md` — Template for adding per-residue charts
- `metrics/lis-formulas.md` — LIS/iLIS calculation (especially PAE bidirectional averaging)
- `workflows/github-deploy.md` — Always `git push vis main`, never bare `git push`

**For feature planning:**
- `ideas/planned-features.md` — Current roadmap
- `ideas/papers.md` — Papers to reference

**For page-specific context:**
- `pages/universal.md`, `pages/flypredictome.md`, `pages/dimer.md`, `pages/ortholog.md`, `pages/monomer.md`

## Git

- **Active remote:** `vis` → `github.com/flyark/LIVIA`
- **Push command:** `git push vis main` (NEVER bare `git push` — origin is archived)
- **Do NOT push** until explicitly approved by user

## Key Rules

- PAE must always be bidirectionally averaged: `(pae[i][j] + pae[j][i]) / 2`
- Never hardcode CIF column indices — use dynamic `_atom_site.` header parsing
- Never hardcode chain colors — use `getColorsForPair()` or color picker values
- Legends go outside canvas as HTML divs, never drawn on canvas
- Check all pages for consistency when adding features
