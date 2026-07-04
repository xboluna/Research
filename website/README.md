# Searching for Exploding Black Holes — Interactive Website

An interactive, educational website for the BHRad research program: the search for
evaporating primordial black holes in γ-ray telescope data.

Based on:

- **Boluna, Profumo, Blé & Hennings** — *Searching for Exploding Black Holes* ([arXiv:2307.06467](https://arxiv.org/abs/2307.06467))
- **Boluna (2023)** — *Detection Methods for Discovering Evaporating Primordial Black Holes in Modern Gamma-Ray Telescopes*, M.S. thesis, UCSC
- The `BHRad/` research code in this repository

## Stack

- **Next.js 15** (App Router, TypeScript, fully static output)
- **Tailwind CSS v4** — dark-only design system
- **KaTeX** — server-rendered equations
- Custom SVG charts + a Bayer-dithered canvas hero (no chart libraries)

## Deploying to Vercel

Point Vercel at this repository and set:

| Setting | Value |
|---------|-------|
| **Root Directory** | `website` |
| Framework preset | **Next.js** |
| Build command | `npm run build` (default) |
| Output Directory | **leave empty** (do not set `public`) |

`vercel.json` in this directory pins the Next.js framework so Vercel does not treat
the app as a static site.

No environment variables required. The site is fully static.

### Troubleshooting

If the build succeeds but deploy fails with:

> No Output Directory named "public" found

the Vercel project is misconfigured as a static site. Fix:

1. **Project Settings → Build & Development**: Framework = Next.js, Output Directory blank.
2. Remove any **Production Overrides** (yellow banner) for Output Directory.
3. Redeploy. `website/vercel.json` should enforce the correct framework on future builds.

## Local development

```bash
cd website
npm install
npm run dev        # http://localhost:3000
```

## Data pipeline

All charts are driven by static JSON in `src/data/`, generated from the research
repo by `data-pipeline/generate_all.py`:

```bash
pip install -r data-pipeline/requirements.txt
npm run data       # regenerates src/data/*.json from ../BHRad
```

| Output | Source | Physics |
|--------|--------|---------|
| `spectra.json` | Analytical_Modelling.ipynb parameterization | Ukwatta et al. 2016 Eqs. 31–34 |
| `effective-areas.json` | `BHRad/Modelling/EffectiveAreas/*.dat` | digitized instrument responses |
| `dmax-curves.json` | reimplementation of the repo detectability loop | JCAP Eq. 4.3–4.6, Fig. 4 |
| `spectral-index.json` | pipeline | JCAP Eq. 4.8, Fig. 6 |
| `mass-functions.json` | pipeline | JCAP Eqs. 3.9–3.14 |
| `alpha-dof.json` | SM threshold model | dM/dt = −α(M)/M² |
| `candidates.json` | `SGRB_search/H>1_T90[0.2-5]_RS=0.csv` | catalog search output |
| `transients.json` | `TransientSources_fitted_params.csv` | thesis §3.1 fits |
| `grb-lightcurves.json` | `~100ms_Source_Data/*.csv` | real GBM lightcurves |
| `lightcurve.json` | pipeline | universal PBH lightcurve |

The generated JSON is committed, so Vercel builds don't need Python.

## Site structure

Single scroll page with 11 sections (`src/app/page.tsx`):

1. **Hawking** — temperature/mass dial, spectrum explorer, universal lightcurve
2. **Origins** — formation clock, mass function gallery
3. **Collider** — degrees-of-freedom staircase (dark sector toggle), spectral index
4. **Telescopes** — effective area atlas, energy ruler, threeML pipeline
5. **Distance** — sensitivity frontier (regenerated JCAP Fig. 4), solar neighborhood
6. **Fitting** — fit real GBM bursts by hand with the PBH template
7. **Catalogue** — filter playground on the real 36 candidates
8. **Sky** — Mollweide map, isotropy statistics
9. **Afterglow** — multi-messenger timeline
10. **Transients** — rising-source fitter with τ–d degeneracy
11. **Verdict** — the evidence board

Planning documents and build notes live in `agentic-research/`.
