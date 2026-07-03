# Research Notes — Condensed from Papers & Repo

Synthesized from in-repo PDF, LaTeX notes, notebooks, and literature review. For full section mapping see [SECTIONS.md](./SECTIONS.md).

---

## 1. The Cline hypothesis (spine of the project)

The repo extends a decades-long program (Cline et al. 1992–2009) asking whether **a subset of very short, hard GRBs** are **local PBH evaporation fireballs** rather than distant astrophysical events.

**Distinctive predictions:**
1. Duration ≲ 200 ms (fireball photosphere crossing)
2. Hard spectrum; **anti-correlation** of hardness with duration
3. Single-spike light curves
4. Rise time ≲ 1 μs (unresolved)
5. Spatial distribution consistent with **homogeneous Euclidean** local population (V/V_max ~ 1/2)
6. Low event rate among all GRBs (~2% in BATSE 1B–3B)
7. No cosmological redshift

**Repo team's additions (PBH_EBH_Criterion.tex):**
- Fermi GBM/LAT catalogs are larger than BATSE → extend search
- GBM time resolution ~1 ms (not μs) — same domain as BATSE for duration
- Hardness in **GBM/LAT band ratio** is "somewhat consistent" with Cline
- **Isotropy / anomalous sky region** is the best chance for new progress
- Afterglow at same sky position weeks/years later (radio) as secondary test

---

## 2. Hawking radiation & fireball (Cline 1997 §2)

### Basic relations

- Temperature: T_BH ∝ 1/M; 10¹⁴ g → ~100 MeV
- Mass loss: dM/dt = −α(M)/M²
- SM: α ≈ 4.1×10⁻³ (in units used in paper) at high T
- Critical evaporation mass M* ≈ 7×10¹⁴ g; dn/dt|M* ≈ 2.2×10⁻¹⁰ N pc⁻³ yr⁻¹

### Fireball at QGP

When T_BH → T_QGP ~ 160 MeV:
- Quarks/gluons hadronize → pion-rich fireball → photon trapping until optically thin
- Burst duration ~ ms
- **Lower average photon energy** than naive Hawking → easier to see as GRB than as 100 MeV diffuse background
- Fluence at distance R: L_γc / (4π R²) ~ 10⁻⁷ erg cm⁻² → R ~ pc scale

### Figure 1 (α vs M)

Two uncertainty regions:
- **(I)** QGP phase transition — rapid dM/dt
- **(II)** New particle thresholds — increased dof

Short bursts possible at T corresponding to M ~ 10⁹ g or 10¹⁴ g depending on model.

---

## 3. BATSE analysis reproduced in Cline 1997 §3

**Hardness ratio (Eq. 11):**
\[
H = \frac{F(115\text{–}320\ \mathrm{keV})}{F(55\text{–}115\ \mathrm{keV})}
\]

**Selection cuts (BATSE 3B):**
1. T90 < 250 ms (but TTE Gaussian fit duration used; cut > 100 ms applied)
2. Complete TTE + hardness + Cp data
3. Single spike; peak ≥ 2× background
4. Table 1: 11 events; hardness 1.6–7.5; T90 6–97 ms

**Spatial test:** log N vs log Cp slope ≈ −3/2 after Petrosian bias correction → consistent with HISE.

**Key tension:** Weak burst cut removes distant events → sample is isotropic but not fully homogeneous.

---

## 4. Repo forward modeling (Ukwatta / BlackHawk)

`Analytical_Modelling.ipynb` implements **Eq. 31–34** from arXiv:1510.04372:
- `dir(E, M)` — direct Hawking photons
- `frag(E, M)` — fragmentary component post-QCD
- Temperature: T [GeV] ≈ 1.058×10¹³ / M [g]

`SpectraPlotter.py` ingests **BlackHawk** output:
- Primary + secondary photon spectra
- Integrates flux in GBM (50 μeV – 300 μeV in code units) vs LAT (0.1 – 1000 GeV) bands
- PSD / power vs time for burst energetics

---

## 5. Detectability logic in repo

From `Analytical_Modelling.ipynb` and `OrderOfMagEstimations.ipynb`:

**Detection thresholds:**
- ≥ 10 source photons in band, OR
- N_S/√N_B ≥ 5

**Instrument solid angles (approximate):**
- GBM, LAT, BATSE: Ω ~ 10⁻³
- HAWC, VERITAS: ~ 10⁻⁵
- LHAASO: ~ 10⁻⁴

**Distance scaling:** φ ∝ d⁻²; combine with exposure and background models `N_B`, cosmic-ray rejection `gamma_had_sep`.

**Order-of-magnitude:** For φ_γ ∝ M⁻³, telescope floor 10⁻¹⁵ cm⁻² s⁻¹ gives d_max(M) — plot in OrderOfMagEstimations.

---

## 6. Catalog selection (repo-specific)

**Hardness (repo definition):**
\[
H_{\mathrm{repo}} = \frac{F_{\mathrm{LAT\,fluence}}}{F_{\mathrm{GBM\,fluence}}}
\]
(LAT 20 MeV – 300 MeV over GBM 10 keV – 25 MeV)

**Primary filter export:** `H>1_T90[0.2-5]_RS=0.csv`
- Hardness > 1
- T90 ∈ [0.2, 5] s (GBM or LLE)
- No measured redshift

**Note:** Repo T90 window is **wider** than Cline's < 250 ms — intentional extension to Fermi resolution and statistics, with tradeoff in purity.

**Boresight:** θ > 70° + T90×(10°/min) — selects bursts where LAT was not favorably oriented.

---

## 7. Fitting pipeline

**Model:** Power-law burst + adiabatic afterglow (`Fitting.ipynb`)

**Stack:** threeML + astromodels + ultranest (ReactiveNestedSampler)

**Targets:** GRBs with T90 < 2 s; high-time-resolution ~100 ms CSVs for select sources

**Transient branch:** Power-law flux decay → fit T₀ (MJD), distance (pc), photon index — compare to SEM predictions on index vs τ plots.

---

## 8. Radio afterglow (from Afterglow.tex)

| Reference | Key result |
|-----------|------------|
| Rees 1977 | Radio pulses may be **more conspicuous** than γ at ~10⁴ pc |
| Blandford 1977 | ν_c ~ GHz scale; spectral slopes 0.57 / −4 |
| Cutchin 2015 | ETA: no detection; limit ρ̇ < 2.3×10⁻⁷ pc⁻³ yr⁻¹; γ_f ~ 10⁵ (10¹¹ g / M) |

Repo open question: spherical PBH ejecta vs jet blob — same afterglow shape?

---

## 9. Carr & Kühnel constraints (from notes)

- PBHs < 10¹⁵ g mostly evaporated today
- Repo mass 10¹⁴ g is near **M*** evaporation window
- Strong EGB limit at M* ~ 5×10¹⁴ g: f(M) < 2×10⁻⁸ (M/M*)^(3+ε)
- Voyager: M < 10¹⁶ g → f < 0.001
- Allowed asteroid-mass window: 10¹⁶–10¹⁷ g (different from repo target)

---

## 10. "Dark degrees of freedom" — terminology bridge

The user's phrase maps to:

1. **Greybody dof sum** α(M) in dM/dt (Cline 1997 Fig. 1)
2. **Hagedorn exponential spectrum** (historical; probably excluded)
3. **QCD jet + fragmentation** (modern SEM / BlackHawk)
4. **BSM dark sectors** increasing α(M) — explored in repo via Δα contours; see arXiv:2606.13788 for recent theory

**"Cosmic collider"** = PBH emits quark/gluon jets that hadronize like accelerator events, probing particle content without a collider ring.

---

## 11. Items not resolved in repo notes

- Exact coordinates of Cline "anomalous region" / RHS octant
- Full understanding of V/V_max methodology (noted in PBH_EBH_Criterion.tex)
- Stefano's burst-duration method vs catalog T90
- Jet photon contribution when M < M_q ≈ 0.4 M*
- Whether `Completed_fits/` images exist elsewhere for publication

See [OPEN_QUESTIONS.md](./OPEN_QUESTIONS.md).
