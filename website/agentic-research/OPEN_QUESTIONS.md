# Open Questions & Gaps

Items to resolve before or during website implementation. Track resolution here.

---

## Missing documents

| Item | Referenced in | Impact | Action |
|------|---------------|--------|--------|
| ~~**Thesis PDF**~~ | ✅ [Author PDF](https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf) · `BHRad/Papers/Boluna_2023_MS_Thesis_*.pdf` | **Resolved** — see `PAPER_CITATIONS.md` |
| **Johnson PhD thesis** (threeML methodology) | PBH_EBH_Criterion.tex | Secondary threeML reference | Optional |
| `Notes_on_BH_Criterion.pdf` | README | Likely redundant with `.tex` | Confirm or add to repo |
| `Completed_fits/` | README | S6 figure assets | Locate or re-run Fitting.ipynb |
| `Spectrum_Gen&Plot/` | README | S3 BlackHawk figures | Locate or regenerate |
| BlackHawk output folders | BHawk_Modelling, SpectraPlotter | Live spectra interactives | Run BlackHawk or ship precomputed JSON |
| `correct_unassociated/*.csv` | FittingTransients, Analytical_Modelling | S10 transient light curves | Restore from author |

---

## Scientific clarifications for accurate website copy

### 1. Mass target vs "evaporating today"

Repo uses M ~ 10¹⁴ g, but M* (evaporating at cosmic age) ~ 5×10¹⁴ g. **Website must explain:** 10¹⁴ g objects are hotter and shorter-lived remaining time (~10⁴ yr), not necessarily in the final millisecond today.

### 2. Cline cuts vs repo cuts

| Parameter | Cline 1997 | Repo catalog |
|-----------|------------|--------------|
| T90 | < 250 ms (often < 100 ms) | 0.2 – 5 s |
| Hardness bands | BATSE keV | LAT/GBM fluence ratio |

**Decision needed:** Present repo cuts as "extended search" with explicit tradeoff, or add a "strict Cline mode" toggle in S7 interactive.

### 3. Scientific tension / null results

Modern Fermi searches with full SEM templates may disfavor PBH interpretation of Swift sGRBs (literature post-2018). Website should include **Epilogue (S11)** acknowledging:
- Historical Cline evidence is suggestive, not definitive
- Forward modeling (this repo) is the right modern approach
- Cite Ackermann 2018 and recent template searches

### 4. "Dark degrees of freedom" naming

User terminology is pedagogically good but non-standard. Glossary should map to:
- α(M) greybody factor
- QCD / Hagedorn / BSM enhancements

### 5. Cosmic expansion wording

Clarify: expansion enables **formation** in early universe, not ongoing creation of 10¹⁴ g BHs today.

---

## Technical debt for website

| Issue | Detail |
|-------|--------|
| No `requirements.txt` | Pin threeML, ultranest, astropy for reproducibility |
| Notebook paths are machine-specific | `/Users/xboluna/...` in outputs |
| Large data gitignored | FITS catalogs present in some folders only |
| Mathematica `.nb` | Need Python reimplementation or static exports for web |
| PDF figures in LaTeX | `Contraints.png`, `FLDetectability.png` missing from repo |

---

## Sections user listed vs final plan

| User suggestion | Status |
|-----------------|--------|
| Hawking & explosion | ✅ S1 |
| Small BHs & cosmic expansion | ✅ S2 (wording refined) |
| Dark dof / cosmic collider | ✅ S3 |
| Detectors & threeML | ✅ S4 |
| Detection distance | ✅ S5 |
| Fitting sources | ✅ S6 |
| Catalogue search | ✅ S7 |
| — Angular isotropy | ➕ S8 (motivated) |
| — Radio afterglow | ➕ S9 (motivated; repo has Afterglow.tex) |
| — LAT transients | ➕ S10 (motivated; major code folder) |
| — Epilogue / tension | ➕ S11 (recommended) |

---

## Questions for paper/thesis author

1. What is the official thesis title and which chapters map to each repo folder?
2. Is there a preferred candidate GRB list for the website hero example?
3. Should the site foreground **confirmation** or **exploration** framing?
4. Are there published figures we must use vs regenerate?
5. Access to BlackHawk runs at M = 10¹⁴ g for SpectraPlotter demos?

---

## Agent / literature follow-ups completed

- ✅ Cline 1997 PDF extracted and mapped
- ✅ PBH_EBH_Criterion.tex and Afterglow.tex read
- ✅ Key notebooks scanned for physics and cuts
- ✅ External literature review: Carr & Kühnel, Ukwatta, Coogan 2026 (dark dof), Ackermann 2018

## Recommended future agent tasks

1. **Thesis ingestion** — when PDF available, map chapters → sections
2. **BlackHawk export** — generate JSON spectra for M = 10¹³, 10¹⁴, 10¹⁵ g
3. **Catalog re-run** — reproduce candidate counts with documented environment
