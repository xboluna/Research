# Interactive Visualizations — Proposals per Section

Design principles for the Next.js site:
- **Mobile-first** but reward desktop with richer canvases
- **Sliders** for physical parameters (M, d, τ); **presets** for cited papers ("Cline 1997 event 432")
- Export notebook plots as JSON/CSV where possible (see REPO_MAP.md data files)
- Use consistent color palette: deep space background, accent for "PBH signal", muted for backgrounds

Tech stack suggestions: **React Three Fiber** or **D3** for sky maps; **Recharts/Visx** for light curves; **KaTeX** for equations; **Framer Motion** for scroll reveals.

---

## Section 0: Prologue

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Discovery timeline** | Horizontal scrubber: 1973 GRB → 1974 Hawking → 1997 Cline → 2008 Fermi → BlackHawk | Static dates |
| **"Two populations?"** | Toggle overlay: cosmological sGRB skymap (isotropic) vs local PBH (isotropic but only nearby visible) | Synthetic |

---

## Section 1: Hawking radiation & explosion

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Temperature–mass slider** | Drag M on log scale; show T_BH, τ, and "detectable by GBM?" badge | `Tbh_mass()` from Analytical_Modelling |
| **Mass vs time** | Animate M(t) as BH evaporates; spike at QGP | `Mbh(time)` |
| **Fireball schematic** | Cross-section: horizon → QGP shell → photosphere → escaped γ | Illustration + Cline 1997 §2.2 |
| **Hawking spectrum** | Plot d²N/dtdE; compare masses 10¹⁷, 10¹³, 10¹⁰ g | `dir()` + `frag()` or static PNG |
| **Cline Table 2 checklist** | Click each row (duration, hardness, etc.) to see BATSE example | Cline 1997 Table 2, Figure 5 |

---

## Section 2: PBH formation & cosmic expansion

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Early universe timeline** | Click epochs: inflation, QCD, BBN; show M_H(t) | Carr & Kühnel review |
| **Constraint band chart** | f(M) allowed windows (simplified); hover for probe name | PBH_EBH_Criterion.tex constraints |
| **J-factor distance toy** | Slider d; show φ_GC/φ_BH ratio | OrderOfMagEstimations.ipynb |
| **Log-normal mass function** | σ, M* sliders → dn/dM plot | Analytical_Modelling.ipynb |

---

## Section 3: Dark degrees of freedom / collider

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **α(M) staircase** | SM thresholds light up as M decreases; optional "dark dof" boost | Cline 1997 Fig. 1; repo `alpha(M)` |
| **Collider analogy** | Split view: LHC diagram ↔ PBH jet fragmentation | MacGibbon & Webber 1990 |
| **Primary vs secondary photons** | Toggle BlackHawk components on SED | SpectraPlotter.py |
| **GBM/LAT flux ratio vs time** | Animated ratio during evaporation | `plot_telescope_flux_ratio()` |
| **Dark sector contour** | 2D heatmap: τ vs Δα/α (detection significance) | Analytical_Modelling dark_signal_mesh |

---

## Section 4: Detectors & threeML

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **A_eff overlay** | Multi-select instruments; log-log E vs A_eff | `EffectiveAreas/*.dat` |
| **Energy band ruler** | Colored bands on log E axis for GBM/LAT/HAWC/LHAASO | `telescope_energy_bands.csv` |
| **Evaporation instrument timeline** | Horizontal bar: which instrument dominates which τ | Detectability.ipynb logic |
| **threeML pipeline** | Animated flowchart: Download → Model → Likelihood → Posterior | Fitting.ipynb structure |
| **Joint fit mockup** | Fake GBM+LLE data with model overlay toggles | Simplified from Fitting.ipynb |

---

## Section 5: Detection distance

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **d_max(τ) per detector** | Log-log plot; toggle detectors | Analytical_Modelling contours |
| **Solar neighborhood map** | 1 pc, 0.02 pc spheres around Sun; "GBM vs LAT reach" | Static geometry |
| **N_S vs distance slider** | User picks τ, M; see if N_S>10 or S/N>5 | `instantaneous_N_S()` |
| **Proper motion streak** | Animate 4-year LAT track for nearby PBH | propermotion.dat |
| **Data overlay mode** | Plot real transients + GRB candidates on d_max curves | CSV files |

---

## Section 6: Fitting sources

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Light curve composer** | Sliders: power-law index, afterglow t_m, t_p | `afterglow()` from Fitting.ipynb |
| **PRFF vs FRED** | Side-by-side templates | Theoretical_Modelling vs typical sGRB |
| **Bayesian posterior mock** | Corner plot (precomputed) for one GRB | ultranest output from Fitting.ipynb |
| **Spectral index vs T₀** | Compare fitted transients to model band | Analytical_Modelling + TransientSources_fitted_params.csv |
| **~100 ms data viewer** | Zoomable light curve for specific GRBs | `~100ms_Source_Data/*.csv` |

---

## Section 7: Catalogue search

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Filter playground** | User sets T90 range, hardness cut, RS=null → live candidate count | Logic from Exploring_possible_candidates |
| **Hardness vs T90 scatter** | Highlight repo selection; Cline dashed trend | Notebook plots |
| **Candidate table** | Sortable GRB names with RA/Dec links to SIMBAD | CSV exports |
| **Boresight explainer** | Diagram: GBM Earth-blocking angle vs T90 | BoresightSelection.py condition |
| **Before/after cuts** | Animate filter pipeline reducing N sources | Synthetic + real counts |

---

## Section 8: Angular distribution

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Mollweide explorer** | Color by hardness or T90; galactic plane overlay | BoresightSelection / catalog notebooks |
| **Hemisphere balance** | Live N/S, E/W ratios vs isotropic expectation | `generate_mollweise()` stats |
| **Cline anomalous region** | Highlight RHS octant from Cline 2002 | Literature |
| **V/V_max primer** | Simple 2D simulation: homogeneous vs disk population | Educational toy |

---

## Section 9: Radio afterglow

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Multi-messenger clock** | Log time axis: μs γ → days radio | Afterglow.tex |
| **Blandford spectrum** | ν_c slider with γ_f, B | Blandford 1977 equations |
| **Sky position matcher** | Click γ-burst position; search cone for late radio | Conceptual |
| **Afterglow in γ fit** | Toggle afterglow component on fitted GRB | Fitting.ipynb |

---

## Section 10: LAT transients

| Interactive | Description | Data source |
|-------------|-------------|-------------|
| **Fading transient plot** | Power-law decay with T₀, distance sliders | FittingTransients.ipynb |
| **GRB–transient association** | Map transients + GRBs; expanding 1°/yr ROI | GRBs_in_Transients_ROI.py |
| **Joint skymap** | Transients (gold) vs GRB candidates (red) | Analytical_Modelling mollweide code |

---

## Section 11: Epilogue

| Interactive | Description |
|-------------|-------------|
| **Evidence scale** | Spectrum from "ruled out" → "tentative" → "confirmed" for each signature |
| **What would convince you?** | Checklist users can tick — educational reflection |

---

## Global UI components

- **Parameter HUD**: persistent M, d, τ readout synced across sections
- **"Show the math" drawer**: KaTeX for each interactive
- **Repo link footers**: "See `Analytical_Modelling.ipynb` cell X" (deep links when hosted on GitHub)
- **Citation popovers**: hover paper keys → full reference

---

## Asset extraction tasks (pre-implementation)

1. Run notebooks once to export key PNGs → `website/public/figures/`
2. Serialize `EffectiveAreas/*.dat` → JSON
3. Export candidate CSVs → cleaned JSON for client filter
4. Precompute d_max curves for all detectors → single JSON per τ grid
5. Optional: WASM/Python pyodide for live `dir()+frag()` — otherwise precompute on grid
