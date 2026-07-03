# Website Sections — Narrative, Science, and Repo Anchors

Each section is designed as a scroll chapter with: **hook → physics → observable → repo evidence → interactive**.

Paper citations use short keys; full entries in [BIBLIOGRAPHY.md](./BIBLIOGRAPHY.md).

---

## Section 0: Prologue — Why look for exploding black holes?

**Hook:** Most GRBs come from cosmological distances. A tiny subclass might be explosions in our cosmic backyard.

**Core ideas:**
- GRBs were discovered (1973) around when Hawking radiation (1974) was proposed — a natural but long-contested association.
- Conventional GRBs have complex multi-pulse light curves; PBH evaporation predicts single-spike, hard, sub-second bursts.
- This repo extends the Cline et al. program with Fermi's larger catalogs and forward modeling (BlackHawk + threeML).

**Repo anchors:**
- `README.md` — project overview
- `BHRad/Papers/Evidence for GRBs consistent with PBH evaporation.pdf` — Cline, Sanders & Hong 1997 (ApJ 486, 169)
- `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` — §Main conclusions

**Key paper quotes (Cline 1997, §1):**
> "It is still possible that there is a sizable density of PBHs in our Galaxy and that some of the GRBs could be due to PBH evaporation."

**Interactive idea:** Timeline slider from 1973 (GRB discovery) → 2026 (Fermi era + BlackHawk).

---

## Section 1: Hawking radiation & the black hole explosion (overview)

**Hook:** A black hole is not perfectly black — it radiates like a thermal body whose temperature rises as it shrinks, ending in a violent final flash.

### Physics

**Hawking temperature:**
\[
k T_{\mathrm{BH}} = \frac{\hbar c^3}{8\pi G M} \approx \left(\frac{10^{16}\,\mathrm{g}}{M}\right)\,\mathrm{MeV}
\]

For M = 10¹⁴ g: T_BH ~ 100 MeV. For M* ~ 5×10¹⁴ g (evaporating today): T ~ 20 MeV.

**Mass loss (greybody):**
\[
\frac{dM}{dt} = -\frac{\alpha(M)}{M^2}
\]

**Lifetime (standard model, order of magnitude):**
\[
\tau \approx 400\left(\frac{M}{10^{10}\,\mathrm{g}}\right)^3\ \mathrm{s}
\]

**Fireball picture (Cline & Hong 1992; Cline 1997 §2.2):**
- Near the QGP phase transition (T ~ 100–200 MeV), quarks/gluons hadronize into an expanding fireball.
- Photons remain optically thick until the photosphere — burst duration ~ light-crossing time → **milliseconds**.
- Total γ energy at fireball: L ~ 10³³–10³⁴ erg for M ~ 10¹⁴ g.
- Fluence constraint for BATSE detection: L_γ / (4π R²) ~ 10⁻⁷ erg cm⁻² → **R ≲ few pc**.

### Observable signatures (Table 2, Cline 1997)

| Property | PBH expectation |
|----------|-----------------|
| Duration | ≪ 1 s; fireball picture *q ~ 200 ms |
| Spectrum | Hard γ; anti-correlation: shorter → harder |
| Light curve | Single spike |
| Rise time | ≲ 1 μs (below GBM resolution) |
| Afterglow (γ) | None from vacuum evaporation |
| Spatial distribution | Homogeneous & isotropic if local (V/V_max ~ 1/2) |

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/Modelling/Theoretical_Modelling.ipynb` | Simplified light curve (power-law index ~ −0.52 + afterglow); hardness; T90 |
| `BHRad/Modelling/Analytical_Modelling.ipynb` | Analytical photon flux `dir()` + `frag()` from Ukwatta et al. 2016 (arXiv:1510.04372); T(M) |
| `BHRad/Modelling/BHawk_Modelling.ipynb` | BlackHawk primary/secondary spectra |
| `BHRad/Modelling/Blackhawk_Modelling/SpectraPlotter.py` | SED, PSD, mass evolution, GBM/LAT flux ratio |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | §Burst duration, §Hardness, §Rise time |

### Papers
- Hawking 1974, 1975 — quantum emission
- Cline & Hong 1992 (ApJ 401, L57) — unique detection criteria
- Cline & Hong 1996 (ApJ 464, L79) — hardness vs duration
- **Cline, Sanders & Hong 1997** — main repo PDF
- Ukwatta et al. 2016 — modern SEM spectra (arXiv:1510.04372)

---

## Section 2: Small black holes & how cosmic expansion could create them

**Hook:** The universe may have forged black holes in its first fraction of a second — some could be evaporating near us right now.

### Physics

**Horizon mass at formation time t:**
\[
M_H(t) \approx 10^{15}\left(\frac{t}{10^{-23}\,\mathrm{s}}\right)\,\mathrm{g}
\]

Cosmic expansion sets the Hubble horizon; overdensities above the Jeans threshold can collapse into PBHs.

**Formation mechanisms** (Carr & Kühnel 2020; notes in `PBH_EBH_Criterion.tex` §Carr and Kuhnel):
1. Inflationary density fluctuations (δ > δ_c ~ 0.45 in radiation era)
2. Critical collapse — broad low-mass tail
3. QCD / electroweak phase transitions — enhanced abundance at specific masses
4. Cosmic string loops, bubble collisions (speculative)

**Mass function & constraints:**
- PBHs with M < 10¹⁵ g have largely evaporated by today.
- Repo target M ~ 10¹⁴ g is **hotter than M*** but still has τ ~ 10⁴–10⁵ yr remaining — not necessarily in the final millisecond.
- Allowed windows for PBH dark matter are narrow; evaporation constraints from Voyager positrons, extragalactic γ background (EGB), etc.

**Local density argument (repo):**
- `OrderOfMagEstimations.ipynb`: if one PBH sits at d ~ 0.015 pc, compare flux to galactic-center signal via J-factor scaling.
- Number density n ~ (4/3 π d³)⁻¹ for "one object in a sphere of radius d."

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/OrderOfMagEstimations.ipynb` | M_BH, ρ_BH, f_DM, J_GC, distance for equal flux |
| `BHRad/Modelling/Analytical_Modelling.ipynb` | Log-normal PBH mass function; evolved mass function with α(M) |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | §Carr and Kuhnel — formation, constraints figure refs |
| `BHRad/TransientLATSources/predictions_number_PBH.nb` | Mathematica number-density predictions |

### Papers
- Carr & Kühnel 2020 (arXiv:2006.02838) — review
- Jedamzik 1997 — QCD epoch PBHs
- Meszaros 1975 — Poisson fluctuations / structure seeding (mentioned in notes)

### Website note
Clarify for readers: **"cosmic expansion" here means PBH formation tied to the expanding early universe**, not expansion *creating* black holes today. Today's evaporation is mass-loss physics, not Hubble expansion.

---

## Section 3: Dark degrees of freedom — an exploding PBH as a cosmic-scale collider

**Hook:** As a black hole heats up, it emits every particle species it is hot enough to produce — like a particle collider whose beam energy is set by gravity alone.

### Physics

**Running greybody factor α(M)** counts effective particle degrees of freedom. In the Standard Model, α steps up at each mass threshold. Near the QGP transition, models disagree:

| Model | Behavior of α(M) near ~150 MeV |
|-------|-------------------------------|
| Standard Model + QCD jets | Step-like increase; fragmentation photons |
| Hagedorn (historical) | Exponential growth ρ(m) ∝ m^−β exp(m/Λ) — microsecond bursts |
| Beyond-SM "dark dof" | Extra contributions to α → faster evaporation, altered photon yield |

**Figure 1 (Cline 1997):** Shows α(M) with uncertainty regions (I) QGP transition, (II) new particle thresholds. Intense short bursts possible when dM/dt spikes.

**"Cosmic collider" framing:**
- Quark/gluon jets hadronize like accelerator jets (MacGibbon & Webber 1990).
- Probes energy scales up to T_BH without a beam pipe — sensitive to **total dof count**, not portal couplings.
- Repo explores **dark-sector perturbations to α** in detectability contours (`Analytical_Modelling.ipynb` dark_signal_mesh).

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/Modelling/Analytical_Modelling.ipynb` | `alpha(M)`; dark-sector Δα/α contour plots; spectral index vs mass |
| `BHRad/Modelling/BHawk_Modelling.ipynb` | Full fragmentation spectra from BlackHawk |
| `BHRad/Modelling/Blackhawk_Modelling/SpectraPlotter.py` | Primary vs secondary photons; `plot_telescope_flux_ratio()` |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | Hagedorn notes; jet photon question at M_q |

### Papers
- MacGibbon & Webber 1990 (PRD 41, 3052) — QCD jets from PBHs
- Arbey et al. 2019 — BlackHawk code (arXiv:1905.04268)
- Ukwatta et al. 2016 — SEM observational characteristics
- Coogan, Morrison & Profumo 2026 (arXiv:2606.13788) — dark dof & mass-function tail (external; recommended reading)

---

## Section 4: High-energy detectors & threeML as a multispectral framework

**Hook:** No single telescope sees the whole evaporation story — MeV, GeV, TeV, and radio each probe different stages.

### Instrument landscape (repo effective areas)

| Detector | Energy band (repo) | Peak A_eff (order) | Role for ~10¹⁴ g PBH |
|----------|-------------------|-------------------|----------------------|
| **Fermi GBM** | ~8 keV – 40 MeV | ~10² cm² | Fireball / MeV burst — **primary** |
| **Fermi LAT** | 0.1 – 1000 GeV | ~0.5–1 m² | Months–years evaporation; proper motion |
| **BATSE** (historical) | keV – MeV | — | Cline's original catalog |
| **HAWC** | ~0.3 – 100 TeV | 10⁵–10⁹ cm² | Final sub-second TeV tail |
| **LHAASO** | ~20 GeV – PeV | 10⁶–10¹⁰ cm² | Highest-energy photons |
| **VERITAS** | TeV | — | Final-stage ground arrays |

Data files: `BHRad/Modelling/EffectiveAreas/*.dat`, `telescope_energy_bands.csv`.

### threeML workflow

**threeML** = multi-mission maximum likelihood framework used throughout BHRad:

1. **Download:** `FermiGBMBurstCatalog`, `FermiLLEBurstCatalog`, `TimeSeriesBuilder` (`Lightcurve_Download.ipynb`, `Download.py`)
2. **Model:** astromodels spectral + temporal components (power-law burst + afterglow)
3. **Fit:** ultranest nested sampling (`Fitting.ipynb`, `Generate_Images.py`)
4. **Compare:** joint GBM + LAT (LLE) likelihood vs astrophysical GRB templates

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/Modelling/Detectability.ipynb` | LAT, HAWC, LHAASO detectability |
| `BHRad/Modelling/Analytical_Modelling.ipynb` | Overlay all A_eff curves; N_S vs lifetime |
| `BHRad/Lightcurve_Download/Lightcurve_Download.ipynb` | GBM TTE/cspec + LAT LLE download |
| `BHRad/Lightcurve_Download/Jan9Workaround/Download.py` | Standalone download helpers |
| `BHRad/Lightcurve_Fitting/Fitting.ipynb` | Main Bayesian fitting pipeline |
| `BHRad/Lightcurve_Fitting/AfterglowExamples/Generate_Images.py` | Example fit for GRB150902733 |
| `BHRad/Modelling/Blackhawk_Modelling/SpectraPlotter.py` | `apply_telescope_bands()`, GBM/LAT ratio |

### Papers
- Atwood et al. 2009 — LAT instrument
- Ackermann et al. 2012 — LAT performance
- Vianello et al. 2015 — HAWC
- Vasiliev et al. 2015 — threeML ([docs](https://threeml.readthedocs.io))
- Ritz, Atwood, Omodei — Fermi EBH search (cited in PBH_EBH_Criterion.tex)

---

## Section 5: Constraints on detection distance

**Hook:** A PBH flash is bright but brief — only objects within parsecs (or milliparsecs at GeV) are detectable.

### Physics

Detected photon count scales as:
\[
N_S \propto \frac{1}{4\pi d^2} \int \frac{dN_\gamma}{dE} A_{\mathrm{eff}}(E)\, dE
\]

**Repo detection criteria** (`Analytical_Modelling.ipynb`):
- N_S ≥ 10 photons, **or**
- N_S / √N_B ≥ 5 (signal-to-background)
- Solid angle factor dΩ per instrument

**Representative limits:**

| Stage | Instrument | r_max (order) |
|-------|------------|---------------|
| MeV fireball | GBM | ~1–10 pc |
| GeV evaporation | LAT | ~0.02 pc |
| TeV final burst | HAWC | ~0.015 pc |
| Proper-motion streak | LAT (4 yr) | ~0.02 pc at 16 GeV |

**Distance vs mass** (`OrderOfMagEstimations.ipynb`):
\[
d_\mathrm{max} \sim \sqrt{\frac{4\pi \times 10^{-15}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}}{|\phi_\gamma|}}
\]
with φ_γ ∝ M⁻³ scaling from luminosity.

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/OrderOfMagEstimations.ipynb` | d(M) plot; flux ratio φ_GC/φ_BH |
| `BHRad/Modelling/Analytical_Modelling.ipynb` | Distance contours per detector; proper-motion limits |
| `BHRad/Modelling/Detectability.ipynb` | Telescope-specific detectability |
| `BHRad/Modelling/PBH_detectability.nb` | Mathematica counterpart |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | Fermi-LAT ≤ 0.02 pc for 16 GeV; FLDetectability.png |

### Overlay data products (for interactives)
- Fitted LAT transients plotted on detection curves (`TransientSources_fitted_params.csv`)
- Selected GRBs from catalog (`selectedGBMCatalogGRBs.csv`)

---

## Section 6: Fitting sources — parameters & templates

**Hook:** Is a candidate burst shaped like an evaporating black hole, or like a conventional short GRB?

### Model components (repo)

**Primary burst (PRFF-like):**
- Power-law rise/fall: dN/dt ∝ τ^α with α ~ −0.52 (Theoretical_Modelling.ipynb, Fitting.ipynb)

**Afterglow (adiabatic blob):**
```python
# Fitting.ipynb — afterglow component
def afterglow(tau, delta, t_m, t_p, normalization):
    ...
```
Inspired by Cutchin et al. 2021 blob expansion; repo asks whether spherical PBH ejecta produce similar late emission.

**Bayesian parameters (ultranest):**
- Burst normalization, power-law index, afterglow turn-on (t_m, t_p), spectral indices
- Separate priors for GBM vs LAT (LLE) channels (`prior_14_GBM`, `prior_14_LLE`, etc.)

**Transient fits** (`FittingTransients.ipynb`):
- Power-law decay in flux vs time → extract T₀ (evaporation epoch) and distance
- Fit quality via reduced χ² with degrees of freedom

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/Lightcurve_Fitting/Fitting.ipynb` | Batch fits for T90 < 2 s; ultranest |
| `BHRad/Lightcurve_Fitting/AfterglowExamples/Generate_Images.py` | GRB150902733 example |
| `BHRad/Lightcurve_Fitting/~100ms_Source_Data/*.csv` | High-time-resolution light curves |
| `BHRad/TransientLATSources/FittingTransients.ipynb` | LAT transient power-law fits |
| `BHRad/TransientLATSources/TransientSources_fitted_params.csv` | T₀, distance, photon index |
| `BHRad/Modelling/Theoretical_Modelling.ipynb` | Analytical template + 10% peak crossing |
| `BHRad/Papers/Notes/Afterglow.tex` | Blandford '77, Rees '77, Cutchin '15 |

### Distinguishing PBH from astrophysical sGRB

| Feature | PBH template | Astrophysical sGRB |
|---------|--------------|-------------------|
| Pulses | Single | Multi-peaked FRED |
| Spectral evolution | Soft → hard (T rises) | Hard → soft |
| Duration | ms–100 ms (fireball) | 0.1–2 s common |
| Redshift | None (local) | Often extragalactic |
| Afterglow | Radio (weeks), not X-ray sGRB afterglow | X-ray/optical/radio |

---

## Section 7: Catalogue search & candidate selection

**Hook:** Start with thousands of GRBs — filter to the few that look like nearby explosions.

### Selection pipeline (repo)

**Data sources:**
- Fermi 2nd LAT GRB Catalog (`FullCatalog.fits`, `RS=null.fits`)
- GBM burst/trigger catalogs (`GBMCatalog.fits`, `GBMTriggerCatalog.fits`, `LLECatalog.fits`)
- Merged table: `GBM_BurstTrig_merge.csv` (`BoresightSelection.py`)

**Cuts (Cline-inspired, extended in repo):**

| Criterion | Cline 1997 (BATSE 3B) | Repo (`Exploring_possible_candidates.ipynb`) |
|-----------|----------------------|-----------------------------------------------|
| Duration | T90 < 250 ms (TTE fit) | T90 ∈ [0.2, 5] s (GBM and/or LLE) |
| Hardness | F(115–320 keV)/F(55–115 keV) | LAT fluence / GBM fluence > 1 |
| Redshift | — | No redshift (local) |
| Light curve | Single spike | High GBM TS; GBM-LAT timing |
| Sky position | Isotropic | Mollweide maps; boresight θ cut |
| Weak bursts | Peak ≥ 2× background | Explored in BoresightSelection |

**Hardness definition in repo:**
```python
# Exploring_possible_candidates.ipynb
targets['hardness'] = targets['like_lat_fluence'] / targets['like_gbm_fluence']
# LAT: 20 MeV – 300 MeV; GBM: 10 keV – 25 MeV
```

**Output catalogs:**
- `H>1_T90[0.2-5]_RS=0.csv`
- `BOTH_T90_in_[0.2-5].csv`
- `Single_bin_<2s_sources.csv`
- `<2s_t90_dataset.csv`

**Boresight cut** (`BoresightSelection.py`):
```python
GBM['condition'] = (GBM['THETA'] > (70 + GBM['T90']*(1/6)))  # degrees
```
Selects bursts where the LAT would not have been boresighted — explains GBM-only or delayed LAT coverage.

### Repo anchors

| Asset | Content |
|-------|---------|
| `BHRad/GBM_Catalog_Searching/SGRB_search/Exploring_possible_candidates.ipynb` | Main catalog analysis |
| `BHRad/GBM_Catalog_Searching/BoresightSelection.py` | Merge catalogs, Mollweide maps, T90/T50 |
| `BHRad/GBM_Catalog_Searching/Before_Sensitivity_Roadblock/Widened_GRB_search.ipynb` | Earlier broader search |
| `BHRad/GBM_Catalog_Searching/SGRB_search/Stefano_Theory/*.nb` | Hardness, T90, anisotropy theory |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | §Isotropy, §V/Vmax, Cline paper notes |

### Papers
- Cline 1997 §3 — BATSE analysis methodology
- Cline 2002 (ApJ) — galactic anomalous region (20 excess sources)
- Cline 2009 — final evaporation evidence

---

## Section 8: Angular distribution & isotropy *(motivated addition)*

**Why add this section:** The Cline program's most distinctive claim is not hardness alone — it is that very short bursts appear **isotropic and locally homogeneous** (V/V_max ~ 1/2), unlike Galactic neutron-star populations.

**Physics:**
- PBHs in the Galactic halo would be spatially correlated with the halo, but **only nearby** events are detectable → apparent isotropy on the sky.
- Cline 2002: ~20 excess sources in one octant ("anomalous region") while avoiding Galactic plane/center.
- Test: log N vs log C_p slope ≈ −3/2 for homogeneous Euclidean distribution.

**Repo anchors:**
- `BoresightSelection.py` — `generate_mollweise()`, center/anticenter ratios
- `Exploring_possible_candidates.ipynb` — Mollweide hardness/T90 maps
- `Analytical_Modelling.ipynb` — `GRB_and_transients_source_dist.png` production
- `PBH_EBH_Criterion.tex` — §Isotropy/angular distribution

**Interactive:** Rotating Mollweide map with user-adjustable sky density; compare to Galactic plane overlay.

---

## Section 9: Radio afterglow & multi-messenger lag *(motivated addition)*

**Why add this section:** Repo afterglow work spans γ-ray fitting *and* radio transient physics — a rich multi-messenger story for the website.

**Physics:**
- **Rees 1977:** e⁺e⁻ pairs from explosion → relativistic conducting shell → synchrotron radio pulse; detectable to ~10⁴ pc.
- **Blandford 1977:** Spectrum ν_c ~ γ_f^(8/3) B^(2/3) E^(1/3) GHz.
- **Cutchin 2015:** ETA null → ρ̇_PBH < 2.3×10⁻⁷ pc⁻³ yr⁻¹; arrives weeks–years after γ burst.
- Repo question (`Afterglow.tex`): Does spherical PBH ejecta produce the same blob expansion as jet afterglows?

**Repo anchors:**
- `BHRad/Papers/Notes/Afterglow.tex`
- `BHRad/Lightcurve_Fitting/Fitting.ipynb` — afterglow component in γ fit
- `BHRad/Lightcurve_Fitting/AfterglowExamples/Generate_Images.py`

**Interactive:** Dual timeline — ms γ burst vs weeks–years radio afterglow at same sky position.

---

## Section 10: LAT transients & GRB association *(motivated addition)*

**Why add this section:** Independent path to local PBHs — long-lived GeV transients and spatial/temporal association with short GRBs.

**Physics:**
- Unassociated LAT sources with power-law fading may be months–years evaporation stage.
- Fit T₀ and distance from flux evolution (`FittingTransients.ipynb`).
- Search for GRBs within expanding ROI: 1°/yr from last transient detection (`GRBs_in_Transients_ROI.py`).

**Repo anchors:**
- `BHRad/TransientLATSources/FittingTransients.ipynb`
- `BHRad/TransientLATSources/Modelling.ipynb` — lifetime τ modeling
- `BHRad/TransientLATSources/GRBs_in_Transients_ROI.py`
- `BHRad/TransientLATSources/GRBs_Transients_proximity.csv`

**Paper:** [Transient LAT Sources](https://iopscience.iop.org/article/10.3847/1538-4365/ac072a) (referenced in Modelling.ipynb).

---

## Section 11: Epilogue — Open questions & scientific tension *(recommended)*

**Why add this section:** Educational integrity — the field is contested.

- Historical Cline association (BATSE short bursts) vs modern template searches with full SEM (some null results).
- Repo extends with forward modeling rather than phenomenology alone — position this clearly.
- Missing from repo: thesis PDF, `Completed_fits/`, full BlackHawk output folders.

See [OPEN_QUESTIONS.md](./OPEN_QUESTIONS.md).

---

## Suggested narrative order on the website

```
Prologue → Hawking & explosion → PBH formation → Dark dof / collider
    → Detectors & threeML → Distance limits → Catalogue search
    → Fitting → Isotropy → Radio afterglow → LAT transients → Epilogue
```

Alternative: interleave **Catalogue search** after **Detectors** so readers see *why* multispectral coverage matters before filtering.
