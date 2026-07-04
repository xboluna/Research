/** GitHub URLs for the xboluna/Research monorepo (BHRad + website). */
export const REPO_ROOT = "https://github.com/xboluna/Research";
export const REPO_BRANCH = "master";

export function repoUrl(path: string) {
  return `${REPO_ROOT}/tree/${REPO_BRANCH}/${path}`;
}

export function repoFileUrl(path: string) {
  return `${REPO_ROOT}/blob/${REPO_BRANCH}/${path}`;
}

/** Well-known paths cited throughout the site. */
export const REPO_PATHS = {
  bhrad: "BHRad",
  analyticalModelling: "BHRad/Modelling/Analytical_Modelling.ipynb",
  effectiveAreas: "BHRad/Modelling/EffectiveAreas",
  lightcurveData: "BHRad/Lightcurve_Fitting/~100ms_Source_Data",
  catalogSearch: "BHRad/GBM_Catalog_Searching/SGRB_search/Exploring_possible_candidates.ipynb",
  candidatesCsv: "BHRad/GBM_Catalog_Searching/SGRB_search/H>1_T90[0.2-5]_RS=0.csv",
  boresightSelection: "BHRad/GBM_Catalog_Searching/BoresightSelection.py",
  fittingPipeline: "BHRad/Lightcurve_Fitting/Fitting.ipynb",
  transientFits: "BHRad/TransientLATSources/FittingTransients.ipynb",
  afterglowNotes: "BHRad/Papers/Notes/Afterglow.tex",
  dataPipeline: "website/data-pipeline/generate_all.py",
} as const;
