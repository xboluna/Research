# Repository → Website Section Map

Quick index from code assets to website sections. Paths relative to repo root.

## Legend

| Tag | Meaning |
|-----|---------|
| **S1** | Hawking & explosion |
| **S2** | Small BH & formation |
| **S3** | Dark dof / collider |
| **S4** | Detectors & threeML |
| **S5** | Detection distance |
| **S6** | Fitting sources |
| **S7** | Catalogue search |
| **S8** | Isotropy |
| **S9** | Radio afterglow |
| **S10** | LAT transients |

---

## Python modules

| File | Sections | What it does |
|------|----------|--------------|
| `BHRad/Modelling/Blackhawk_Modelling/SpectraPlotter.py` | S1, S3, S4 | Load BlackHawk outputs; SED/PSD; GBM/LAT flux ratio; telescope bands |
| `BHRad/GBM_Catalog_Searching/BoresightSelection.py` | S7, S8 | Merge GBM catalogs; boresight θ cut; Mollweide maps; T90/T50 vs fluence |
| `BHRad/Lightcurve_Download/Jan9Workaround/Download.py` | S4, S6 | threeML GBM/LLE lightcurve download |
| `BHRad/Lightcurve_Fitting/AfterglowExamples/Generate_Images.py` | S4, S6, S9 | threeML + ultranest fit for GRB150902733 |
| `BHRad/TransientLATSources/GRBs_in_Transients_ROI.py` | S7, S10 | Cross-match GRBs with LAT transients; proximity plots |

---

## Jupyter notebooks

| Notebook | Sections | Key outputs / concepts |
|----------|----------|------------------------|
| `BHRad/OrderOfMagEstimations.ipynb` | S2, S5 | d(M); φ_GC/φ_BH; Stefano scaling; J-factor |
| `BHRad/Modelling/Analytical_Modelling.ipynb` | S1, S3, S4, S5, S8 | `dir()`, `frag()` spectra; α(M); A_eff overlays; distance contours; dark-sector Δα; hardness index |
| `BHRad/Modelling/Theoretical_Modelling.ipynb` | S1, S6 | Analytical light curve; hardness; T90; anisotropy |
| `BHRad/Modelling/BHawk_Modelling.ipynb` | S1, S3 | BlackHawk spectrum processing |
| `BHRad/Modelling/Detectability.ipynb` | S4, S5 | LAT/HAWC/LHAASO detectability |
| `BHRad/GBM_Catalog_Searching/SGRB_search/Exploring_possible_candidates.ipynb` | S7, S8 | Hardness vs T90; RS=null; candidate CSVs; Mollweide |
| `BHRad/GBM_Catalog_Searching/Before_Sensitivity_Roadblock/Widened_GRB_search.ipynb` | S7 | Earlier broad GRB search |
| `BHRad/Lightcurve_Download/Lightcurve_Download.ipynb` | S4, S6 | Fermi data download via threeML |
| `BHRad/Lightcurve_Fitting/Fitting.ipynb` | S4, S6, S9 | ultranest Bayesian fits; afterglow model; batch T90<2s |
| `BHRad/TransientLATSources/FittingTransients.ipynb` | S5, S6, S10 | Power-law transient fits; T₀, distance |
| `BHRad/TransientLATSources/Modelling.ipynb` | S10 | Lifetime τ; transient paper link |

---

## Data files (for future static import)

| File | Sections | Description |
|------|----------|-------------|
| `BHRad/Modelling/EffectiveAreas/*.dat` | S4, S5 | GBM, LAT, BATSE, HAWC, LHAASO, VERITAS A_eff |
| `BHRad/Modelling/EffectiveAreas/propermotion.dat` | S5 | LAT proper-motion sensitivity |
| `BHRad/Modelling/Blackhawk_Modelling/telescope_energy_bands.csv` | S4 | Instrument energy bands |
| `BHRad/GBM_Catalog_Searching/SGRB_search/H>1_T90[0.2-5]_RS=0.csv` | S7 | Filtered candidates |
| `BHRad/GBM_Catalog_Searching/SGRB_search/BOTH_T90_in_[0.2-5].csv` | S7 | Both GBM & LLE T90 in range |
| `BHRad/GBM_Catalog_Searching/GBM_BurstTrig_merge.csv` | S7 | Merged GBM catalog |
| `BHRad/TransientLATSources/TransientSources_fitted_params.csv` | S5, S6, S10 | Fitted transient parameters |
| `BHRad/TransientLATSources/GRBs_Transients_proximity.csv` | S10 | Angular proximity table |
| `BHRad/Lightcurve_Fitting/~100ms_Source_Data/*.csv` | S6 | High-cadence GRB light curves |

---

## LaTeX & PDF notes

| File | Sections | Content |
|------|----------|---------|
| `BHRad/Papers/Evidence for GRBs consistent with PBH evaporation.pdf` | S1, S7, S8 | Cline, Sanders & Hong 1997 — primary paper |
| `BHRad/Papers/Notes/PBH_EBH_Criterion.tex` | S1–S8 | Literature review; selection criteria; Carr & Kühnel notes |
| `BHRad/Papers/Notes/Afterglow.tex` | S6, S9 | Blandford, Rees, Cutchin radio notes |

---

## Mathematica notebooks

| File | Sections | Content |
|------|----------|---------|
| `BHRad/Modelling/PBH_detectability.nb` | S5 | Detectability calculations |
| `BHRad/GBM_Catalog_Searching/SGRB_search/Stefano_Theory/BH_hardness_ratios.nb` | S1, S7 | Hardness theory |
| `BHRad/GBM_Catalog_Searching/SGRB_search/Stefano_Theory/BH_t90.nb` | S1, S7 | T90 predictions |
| `BHRad/GBM_Catalog_Searching/SGRB_search/Stefano_Theory/GCanisotropy.nb` | S8 | Galactic anisotropy |
| `BHRad/TransientLATSources/predictions_number_PBH.nb` | S2 | PBH number density |

---

## Missing / referenced but not in repo

| Referenced path | Impact on website |
|-----------------|-------------------|
| `Notes_on_BH_Criterion.pdf` | Likely duplicate of `.tex` notes |
| `Completed_fits/` | Fitted-source figures for S6 |
| `Spectrum_Gen&Plot/` | BlackHawk plot pipeline for S3 |
| `Modelling.ipynb` (BHRad root) | Split across `Modelling/` subfolder |
| Thesis PDF (Arne Christian Johnson) | threeML methodology detail for S4 |
| BlackHawk output folders | Raw spectra for S3 interactives |
| `correct_unassociated/` CSVs | LAT transient light curves for S10 |

---

## multicompDM (out of scope v1)

| File | Note |
|------|------|
| `multicompDM/MulticomponentDMNotebook.ipynb` | Tremaine–Gunn u-cDM; Colossus/CAMB — separate physics program |
