# BHRad Interactive Website — Research & Planning

This directory collects planning notes, paper references, and repo-to-section mappings for an interactive Next.js website that explains the **BHRad** research program: searching for signatures of **evaporating primordial black holes (PBHs / EBHs)** in gamma-ray burst data.

## Vision

Build a sleek, science-grounded educational site deployable from this monorepo subdirectory on Vercel. Each section should:

1. Teach a core physical idea with minimal jargon (then deepen on scroll/click).
2. Anchor claims in published literature and in code/notebooks from this repo.
3. Offer at least one interactive graphic (sliders, animated light curves, sky maps, detector overlays).

## Research question (one sentence)

Can very short, hard-spectrum GRBs observed by Fermi GBM/LAT be explained as nearby (~parsec-scale) final evaporation of ~10¹⁴ g primordial black holes, rather than conventional astrophysical short GRBs at cosmological distance?

## Repo context

| Path | Role |
|------|------|
| `BHRad/Modelling/` | Theory: Hawking spectra, BlackHawk, detectability, effective areas |
| `BHRad/GBM_Catalog_Searching/` | Catalog mining, hardness/T90 cuts, sky maps |
| `BHRad/Lightcurve_Fitting/` | threeML + ultranest Bayesian fits of candidate GRBs |
| `BHRad/Lightcurve_Download/` | Fermi data acquisition via threeML |
| `BHRad/TransientLATSources/` | LAT transient association, lifetime/distance fits |
| `BHRad/Papers/` | Thesis, JCAP paper (arXiv:2307.06467), Cline 1997, + extended bibliography |
| `BHRad/OrderOfMagEstimations.ipynb` | Original scaling arguments (distance, J-factor) |

A secondary project `multicompDM/` (ultracompact dark matter) is **out of scope** for the first website version unless we add a sidebar “related work” link.

## Planned website sections

See [SECTIONS.md](./SECTIONS.md) for the full breakdown with repo and paper anchors.

| # | Section | Status |
|---|---------|--------|
| 0 | **Prologue** — Why look for exploding black holes? | Planned |
| 1 | Hawking radiation & the black hole explosion | Planned |
| 2 | Small black holes & cosmic expansion | Planned |
| 3 | Dark degrees of freedom — a cosmic collider | Planned |
| 4 | Multispectral detectors & threeML | Planned |
| 5 | Constraints on detection distance | Planned |
| 6 | Fitting sources (parameters & templates) | Planned |
| 7 | Catalogue search & candidate selection | Planned |
| 8 | **Angular distribution & isotropy** *(added)* | Planned |
| 9 | **Radio afterglow & multi-messenger lag** *(added)* | Planned |
| 10 | **LAT transients & GRB association** *(added)* | Planned |

## Documents in this folder

| File | Purpose |
|------|---------|
| [SECTIONS.md](./SECTIONS.md) | Section-by-section narrative, equations, repo refs, paper refs |
| [REPO_MAP.md](./REPO_MAP.md) | File-level index: which notebook/script supports which claim |
| [BIBLIOGRAPHY.md](./BIBLIOGRAPHY.md) | Annotated paper list with links and relevance tags |
| [INTERACTIVE_VISUALIZATIONS.md](./INTERACTIVE_VISUALIZATIONS.md) | Proposed Next.js interactives per section |
| [RESEARCH_NOTES.md](./RESEARCH_NOTES.md) | Condensed notes from papers and LaTeX files |
| [PAPER_CITATIONS.md](./PAPER_CITATIONS.md) | **Section-level citations** — thesis, JCAP paper, Cline 1997 |
| [OPEN_QUESTIONS.md](./OPEN_QUESTIONS.md) | Gaps, tensions, and remaining open items |

## Key numerical anchors (for consistency across sections)

| Quantity | Value | Source |
|----------|-------|--------|
| Target mass | M ~ 10¹⁴ g (10¹¹ kg) | Repo default; Cline fireball range 6–7×10¹³ g |
| Evaporating-now mass M* | ~5×10¹⁴ g, kT* ~ 20 MeV | Page & Hawking 1976; PBH_EBH_Criterion.tex |
| Hawking temperature | kT_BH ≈ 1.06×10¹³ / M [GeV] | Analytical_Modelling.ipynb |
| Fireball duration | ≲ 50–200 ms (BATSE/GBM resolvable) | Cline 1996–1997 |
| MeV burst distance | ≲ 1–10 pc (GBM); Cline: ≲ 10 pc | Cline 1992; Detectability.ipynb |
| LAT GeV stage distance | ≲ 0.02–0.03 pc | Ritz/Atwood/Omodei; Ackermann 2018 |
| Hardness ratio (pure Hawking limit) | H ~ 250 | Cline 1996 |
| PBH local rate (Fermi-LAT limit) | ≲ 6×10⁴ pc⁻³ yr⁻¹ | PBH_EBH_Criterion.tex (Fermi search) |

## Suggested site architecture (future)

```
website/
├── agentic-research/     ← this folder (planning only)
└── app/                  ← future Next.js app (not yet created)
    ├── page.tsx          ← scroll-driven narrative
    ├── sections/
    └── components/       ← D3/Recharts/Three.js interactives
```

Vercel root directory: `website` (once `app/` exists).

## Primary sources (PDF access verified)

| Document | Location |
|----------|----------|
| **Boluna MS Thesis (2023)** | [xboluna.com PDF](https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf) · `BHRad/Papers/Boluna_2023_MS_Thesis_*.pdf` |
| **Boluna et al. JCAP (arXiv:2307.06467)** | `BHRad/Papers/Boluna_et_al_2024_*.pdf` |
| **Cline et al. 1997** | `BHRad/Papers/Evidence for GRBs consistent with PBH evaporation.pdf` |

See [PAPER_CITATIONS.md](./PAPER_CITATIONS.md) for §/page/equation mapping to website sections.

## Next steps

1. Port key plots from notebooks into static JSON/CSV for client-side interactives.
3. Scaffold Next.js app with section routing matching SECTIONS.md order.
4. Add “scientific honesty” panel: modern null searches vs historical Cline association.
