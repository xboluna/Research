# Build Notes — Implementation log (complete)

Status: **v1 shipped.** All 11 sections implemented, built, and visually verified.

## Architecture decisions (final)

| Decision | Choice | Rationale |
|----------|--------|-----------|
| App location | `website/` is the Next.js root | Vercel root directory = `website` |
| Next.js | 15.5 App Router, TypeScript, static output | Fully prerendered; 179 kB first load |
| Styling | Tailwind CSS v4 (`@theme` tokens in globals.css) | Dark-only |
| Charts | Custom SVG (`src/components/charts/`) | Dithered aesthetic, no chart deps |
| Dithering | Canvas 4×4 Bayer hero (`DitherHero`), SVG/CSS Bayer patterns (`.dither-25/.dither-50`) | User preference, used tastefully |
| Math | KaTeX server-rendered (`ui/Math.tsx`) | zero client JS for equations |
| Data | Python pipeline → committed JSON in `src/data/` | Vercel builds need no Python |
| Fonts | Space Grotesk + IBM Plex Mono (next/font) | |
| Citations | `ui/Cite.tsx` popover with 18 references | Section-level anchors from PAPER_CITATIONS.md |

## Verification performed

- `tsc --noEmit` clean; `next build` clean (static, 4 pages)
- Headless Chrome screenshots of every section, desktop (1440px) + mobile (390px)
- Console error check: **zero errors** (fixed Mollweide pole NaN; added favicon)
- Interaction tests: dark-sector toggle, GRB selector, filter sliders
- Physics sanity: d_max peaks ~0.02–1.3 pc at 10¹⁰–10¹¹ g (matches JCAP Fig. 4 morphology);
  GBM/BATSE dominate ≥10¹⁵ g; spectral index → 1.5 asymptote confirmed in generated data

## Bugs found & fixed during build

1. **dmax tau grid too short** — curves stopped at 1.3e14 g; extended tau to 1e21 s → 1.7e16 g.
2. **Afterglow blew up exponentially** (old shape `u² e^{δ(u−1)/2} e^{−u}` with δ>2 diverges).
   Replaced with gamma-pulse `u^δ e^{−δ(u−1)}` peaking at t_p with value = norm.
3. **Mollweide NaN at poles** — Newton iteration divides by zero at |b|=90°; clamped to ±89.99°
   and guarded the denominator. This was producing SVG path console errors.
4. **TransientFitter flat mock data** — truth τ was 50× the observation window so the rise was
   invisible; retuned to τ=3.5e8 s (≈11 yr) so a decade of data visibly climbs.
5. **404 favicon** — added `src/app/icon.svg` (dithered photon-ring mark).

## Component inventory

| Component | Section | Data |
|-----------|---------|------|
| `DitherHero` | hero | procedural canvas |
| `MassDial` | 01 | closed-form physics.ts |
| `SpectrumExplorer` | 01 | spectra.json |
| `BackwardsBurst` | 01 | lightcurve.json |
| `FormationTimeline` | 02 | closed-form |
| `MassFunctions` | 02 | mass-functions.json |
| `AlphaStaircase` | 03 | alpha-dof.json |
| `SpectralIndex` | 03 | spectral-index.json |
| `DetectorAtlas` | 04 | effective-areas.json + telescope-bands.json + spectra.json |
| `ThreeMLPipeline` | 04 | static content |
| `DistanceFrontier` | 05 | dmax-curves.json + transients.json |
| `Neighborhood` | 05 | closed-form |
| `LightcurveComposer` | 06 | grb-lightcurves.json (real GBM data) |
| `CataloguePlayground` | 07 | candidates.json (real catalog output) |
| `SkyMap` | 08 | candidates.json + transients.json |
| `AfterglowTimeline` | 09 | static content |
| `TransientFitter` | 10 | closed-form Eq. 5.3 |
| `EvidenceBoard` | 11 | curated claims |

## Known scope cuts (deliberate, v2 candidates)

- No CTA effective area (not in repo `EffectiveAreas/`; JCAP Fig. 3 includes it from digitized
  source we don't have) — noted in copy as "future CTA".
- `alpha-dof.json` uses a pedagogical SM threshold model (relative units), not MacGibbon's
  exact greybody-weighted α(M). Correct qualitative staircase; labelled "[rel.]" on axis.
- Hardness–T90 in CataloguePlayground uses the repo's LAT/GBM fluence ratio, not Cline's
  BATSE bands (called out in the caption).
- Mock data in TransientFitter is synthetic (labelled as such); real 1FLT fits appear as
  overlay points in DistanceFrontier instead.
- No OG image yet (metadata is set; image generation would be a nice v2).

## Performance

- First Load JS: 179 kB (all sections); static prerender; no runtime data fetching.
- DitherHero renders at 1/3 resolution with nearest-neighbor upscale — cheap per-frame.
- Largest JSON (grb-lightcurves, 58 kB) is bundled; total data ~176 kB pre-gzip.

## Deploy

Vercel: root directory `website`, no env vars, framework auto-detected.
