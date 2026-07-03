# Paper & Thesis Citation Reference

How to cite primary sources for the BHRad website. Each website section should link to **printed/thesis page numbers** where possible, plus equation and figure numbers.

## PDF access status (verified 2026-07-03)

| Document | Access | Quality | Notes |
|----------|--------|---------|-------|
| **Boluna MS Thesis (2023)** | ✅ [Author PDF](https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf) · [ProQuest](https://www.proquest.com/dissertations-theses/detection-methods-discovering-evaporating/docview/2812309821/se-2) · 📁 `BHRad/Papers/Boluna_2023_MS_Thesis_*.pdf` | **Excellent** — 63/65 pages, ~74k chars, LaTeX PDF | **Primary source** for website narrative |
| **Boluna et al. JCAP draft (arXiv:2307.06467v2)** | ✅ 📁 `BHRad/Papers/Boluna_et_al_2024_Searching_for_Exploding_Black_Holes_arXiv_2307.06467v2.pdf` | **Excellent** — 28/29 pages, ~67k chars | Journal expansion of thesis; same science, more polished |
| **Cline, Sanders & Hong 1997 (ApJ 486, 169)** | ✅ 📁 `BHRad/Papers/Evidence for GRBs consistent with PBH evaporation.pdf` | **Good** — 10/10 pages; ~118 encoding artifacts (`Ð`, `Ï`) from 1997 scan | Historical BATSE evidence; use for context & comparison |

**Extraction method:** `pypdf` text extraction. Equations, figures, and section numbers are reliably recoverable. Cline 1997 needs light cleanup when quoting (`Ð` → `—`, etc.).

**threeML:** Named explicitly in **thesis** (Fig. 2.7 caption, §3.2 analysis). The JCAP paper uses "non-linear least squares" for transient fits and does not name threeML in extracted text — website should cite thesis for the fitting pipeline.

---

## Document hierarchy for the website

```
1. Boluna MS Thesis (2023)     ← authoritative for methods + repo mapping
2. Boluna et al. arXiv:2307.06467  ← authoritative for published equations/figures
3. Cline et al. 1997           ← historical motivation + BATSE selection criteria
4. Repo notebooks/code         ← implementation details & updated exploratory work
```

---

## Boluna MS Thesis (2023)

**Full title:** *Detection Methods for Discovering Evaporating Primordial Black Holes in Modern Gamma-Ray Telescopes*

**Author:** Xavier D. Boluna · **Degree:** M.S. Physics, UCSC · **Date:** March 2023

### Structure (printed page → topic)

| Printed p. | Section | Website section |
|----------|---------|-----------------|
| ix | Abstract | S0 Prologue |
| **Part I — Properties of evaporating PBHs** | | |
| 2–11 | **Ch. 1 Introduction** (§1.1–1.6) | S1, S2 |
| 2 | §1.1 Dark matter | S2 |
| 4 | §1.2 Early universe perturbations | S2 |
| 5 | §1.3 PBHs as dark matter | S2 |
| 7 | §1.4 Black hole evaporation | S1 |
| 8 | §1.5 Photon emission evolution | S1, S3 |
| 11 | §1.6 This work | S0 |
| 12–30 | **Ch. 2 Viability for direct-detection** | S3–S6, S9 |
| 12 | §2.1 Photon lightcurve parameterization | S1, S6 |
| 15 | §2.2 Maximum possible distance | S5 |
| 19 | §2.3 Expected population density and isotropy | S2, S8 |
| 20 | §2.4 Effect of proper motion | S5 |
| 21 | §2.5 Effect of motion on luminosity | S5 |
| 22 | §2.6 Spectral energy index | S3, S6 |
| 24 | §2.7 Advantages in multi-mission detection | S4 |
| 24 | §2.8 Effect of additional **dark degrees of freedom** | **S3** |
| 26 | §2.9 Signal duration vs parameterized index | S6 |
| 28 | §2.10 Possibility for **afterglow** | S9 |
| **Part II — Searching in γ-ray catalogs** | | |
| 32–45 | **Ch. 3 Introduction** (§3.1–3.4) | S7, S10 |
| 33 | §3.1 Analysis of **long-term γ-ray sources** (LAT transients) | S10 |
| 39 | §3.2 Analysis of **short-duration GRB sources** | S7 |
| 40 | §3.3 **Associating transient and GRB sources** | S10 |
| 41 | §3.4 Discussion of results | S11 |
| 46 | **Ch. 4 Conclusions** | S11 |

### PDF page index (for Ctrl+F in bundled PDF)

| Section | PDF page (approx.) |
|---------|-------------------|
| Abstract | 9 |
| §1.4 Black hole evaporation | 18 |
| §2.1 Lightcurve parameterization | 23 |
| BlackHawk discussion | 24–25 |
| §2.2 Max distance | 26 |
| §2.8 Dark degrees of freedom | 35 |
| §2.10 Afterglow | 39 |
| threeML figure/caption | 41 |
| Part II / §3.1 LAT transients | 44 |
| §3.2 GBM GRB catalog | 50 |
| §3.3 Transient–GRB association | 51 |

### Key thesis quotes (for website copy)

**Abstract (p. ix):**
> "We investigate also the possibility and effect of dark particle radiation arising from dark degrees of freedom in black hole evaporation, and the possibility of a multi-spectral afterglow from PBH evaporation products. Lastly, we apply these novel constraints to the Fermi mission catalogs and produce several candidates..."

**§2.8 title:** "Effect of additional dark degrees of freedom" — maps directly to website Section 3 ("cosmic collider").

**threeML (Fig. 2.7, ~p. 30 printed):**
> "Lightcurve retrieval and multi-parameter fitting was performed using ThreeML [52]."

---

## Boluna et al. — *Searching for Exploding Black Holes* (arXiv:2307.06467v2)

**Authors:** Xavier Boluna, Stefano Profumo, Juliette Blé, Dana Hennings · **Prepared for:** JCAP

### Section map (printed/JCAP page from TOC)

| § | Title | PDF p. | Website section |
|---|-------|--------|-----------------|
| 1 | Introduction | 2 | S0 |
| 2 | Photon emission from exploding black holes | 3 | S1, S3 |
| 3 | Density and explosion rates | 5 | S2 |
| 4 | Sensitivity estimates and search strategies | 9 | S4, S5 |
| 4.1 | Spectral energy index | 13 | S3, S6 |
| 4.2 | Effect of proper motion | 14 | S5 |
| 5 | Searches with Fermi LAT and GBM | 14 | S7, S10 |
| 5.1 | Long-term gamma-ray sources | 14 | S10 |
| 5.2 | Short-duration GRB-like sources | 21 | S7 |
| 5.3 | Association of transients and GRBs | 21 | S10 |
| 6 | Discussion and conclusions | 23 | S11 |

### Key equations to cite on website

| Eq. | Content | Website use |
|-----|---------|-------------|
| (2.1)–(2.2) | Hawking greybody + secondary flux | S1 |
| (2.3)–(2.5) | Photon rate & flux vs T_BH, M | S1, S5 |
| (2.6) | dM/dt = −α(M)/M² | S1, S3 |
| (3.8) | ṅ_PBH = ρ_DM ψ_i(M_U) / (3 t_U) | S2 |
| (3.9)–(3.14) | Log-normal, power-law, critical-collapse mass functions | S2 |
| (4.1)–(4.7) | Background, N_S, N_B, detection criteria, rate limits | S4, S5 |
| (4.8) | Spectral index γ(E₀) | S3, S6 |
| (4.9) | Proper motion angle θ | S5 |
| (5.1)–(5.3) | Transient lightcurve fit model F_γ(t) ∝ (τ−t)^−0.533 | S6, S10 |

### Key figures to reproduce interactively

| Fig. | Content | Repo asset |
|------|---------|------------|
| Fig. 1 | f_PBH & explosion rate vs M_* (log-normal) | `Analytical_Modelling.ipynb` |
| Fig. 3 | Effective areas (GBM, LAT, HAWC, …) | `EffectiveAreas/*.dat` |
| Fig. 4 | Max distance vs BH mass per detector | `Analytical_Modelling.ipynb`, `Detectability.ipynb` |
| Fig. 6 | Spectral index γ vs mass | `Analytical_Modelling.ipynb` |
| Table 1 | Fitted LAT transient parameters | `TransientSources_fitted_params.csv` |

### Introduction bullet list (§1, PDF p. 2) — use for S1/S4

Four PBH vs astrophysical GRB differences:
1. Luminosity **never decreases** (universal lightcurve)
2. Universal spectrum (barring extra dof)
3. Must be **local** (≲ few pc); may show **proper motion**
4. **Same intrinsic luminosity** for every event

---

## Cline, Sanders & Hong 1997 (ApJ 486, 169–178)

### Section → PDF page → ApJ page

| Section | PDF p. | ApJ p. |
|---------|--------|--------|
| §1 Introduction | 1 | 169 |
| §2 PBH Evaporation | 1–3 | 169–171 |
| §2.1 QGP phase transition model | 1 | 169 |
| §2.2 PBH fireball | 2–3 | 170–171 |
| §3 BATSE analysis | 4–5 | 172–174 |
| §4 GRB results | 6–8 | 174–177 |
| §5 Conclusion | 9 | 177–178 |

### Artifacts to cite in website

| Ref | Content | Website |
|-----|---------|---------|
| Eq. (1) | Hawking emission rate d²N/dtdE | S1 |
| Eq. (2) | dM/dt = −α(M)/M² | S1, S3 |
| **Fig. 1** | α(M) with QGP regions I & II | **S3** |
| Table 1 | Hardness vs duration (11 BATSE events) | S7 |
| **Table 2** | PBH vs GRB characteristics checklist | S7, S1 |
| Fig. 2 | Hardness anticorrelation with duration | S7 interactive |
| Fig. 4 | Sky distribution (isotropy) | S8 |
| Table 3 | Detection efficiencies vs PBH density model | S5 |

---

## Website section → primary citation (quick lookup)

| Website § | Thesis | JCAP paper | Cline 1997 | Repo |
|-----------|--------|------------|------------|------|
| S0 Prologue | Abstract, §1.6 | §1 | §1 | README |
| S1 Hawking & explosion | §1.4–1.5, §2.1 | §2 | §2, Table 2 | `Theoretical_Modelling.ipynb` |
| S2 Formation & expansion | §1.2–1.3 | §3 | — | `OrderOfMagEstimations.ipynb` |
| S3 Dark dof / collider | **§2.8** | §2, §4.1, Fig. 6 | **Fig. 1** | `Analytical_Modelling.ipynb` |
| S4 Detectors & threeML | **§2.7** | §4, Fig. 3 | — | `EffectiveAreas/`, **thesis Fig. 2.7 threeML** |
| S5 Distance limits | §2.2–2.5 | §4, Fig. 4–5 | §2.2 (≲10 pc) | `Detectability.ipynb` |
| S6 Fitting | §2.1, §2.9, Fig. 2.7 | §5.1, Eq. 5.1–5.3 | — | `Fitting.ipynb` |
| S7 Catalogue search | §3.2 | §5.2 | §3, Table 1 | `Exploring_possible_candidates.ipynb` |
| S8 Isotropy | §2.3 | — | §3, Fig. 4, V/V | `BoresightSelection.py` |
| S9 Afterglow | **§2.10** | — | — | `Afterglow.tex`, `Fitting.ipynb` |
| S10 LAT transients | §3.1, §3.3 | §5.1, §5.3 | — | `TransientLATSources/` |
| S11 Epilogue | §3.4, Ch. 4 | §6 | §4–5 | `OPEN_QUESTIONS.md` |

---

## Citation format for website UI

Use a consistent popover pattern:

```text
Boluna (2023), §2.8, p. 24  — dark degrees of freedom
Boluna et al. (2024), Eq. (4.3), §4  — detection criterion N_S
Cline et al. (1997), Table 2, ApJ 486, p. 176  — PBH vs GRB checklist
```

Link DOIs/URLs:
- Thesis: https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf
- arXiv: https://arxiv.org/abs/2307.06467
- Cline 1997: ApJ 486, 169 (PDF in repo)
