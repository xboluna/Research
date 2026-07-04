#!/usr/bin/env python3
"""
Data pipeline for the BHRad interactive website.

Reads source data from the Research repo (../../BHRad) and emits static JSON
into ../src/data/ for consumption by the Next.js app.

Physics mirrors:
  - BHRad/Modelling/Analytical_Modelling.ipynb  (dir/frag spectra, detectability)
  - Boluna et al., arXiv:2307.06467 (JCAP)      (Eqs. 2.3-2.6, 4.1-4.9, 5.1-5.3)
  - Ukwatta et al. 2016, arXiv:1510.04372       (Eqs. 31-34 parameterization)

Run:  python3 generate_all.py
"""

import json
import math
import os
from glob import glob

import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
BHRAD = os.path.join(REPO, "BHRad")
OUT = os.path.abspath(os.path.join(HERE, "..", "src", "data"))
os.makedirs(OUT, exist_ok=True)

PC_TO_CM = 3.0856775814913673e18
YEAR_S = 3.156e7


def write(name: str, obj) -> None:
    path = os.path.join(OUT, name)
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"))
    print(f"  wrote {name} ({os.path.getsize(path)/1024:.1f} kB)")


# ----------------------------------------------------------------------------
# Hawking photon spectrum: direct + fragmentation
# (Analytical_Modelling.ipynb; Ukwatta et al. Eqs. 31-34; all GeV, seconds)
# ----------------------------------------------------------------------------
A_FRAG = 6.339e23
B_FRAG = 1.1367e24


def x_gamma(E, M):
    return E / (1058.0 * 1e10 / M)


def theta_s(u):
    return 0.5 * (1.0 + np.tanh(10.0 * u))


def Ffunc(y):
    y = np.asarray(y, dtype=float)
    val = np.exp(-0.0962 - 1.982 * (np.log(y) - 1.908) * (1 + np.tanh(20 * (np.log(y) - 1.908))))
    val = np.where(y <= 2.0, 1.0, val)
    return val


def spec_dir(E, M):
    x = x_gamma(E, M)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        out = 1.13e19 * x**6 / ((np.exp(np.clip(x, None, 700)) - 1.0) * Ffunc(x))
    return np.nan_to_num(out, nan=0.0, posinf=0.0)


def spec_frag(E, M):
    x = x_gamma(E, M)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        out = (
            A_FRAG * x ** (-1.5) * (1.0 - theta_s(x - 0.3))
            + B_FRAG * np.exp(-np.clip(x, None, 700)) * theta_s(x - 0.3) / (x * (x + 1.0))
        )
    return np.nan_to_num(out, nan=0.0, posinf=0.0)


def spec_total(E, M):
    return spec_dir(E, M) + spec_frag(E, M)


def T_of_M_GeV(M_g):
    return 1.058e13 / M_g


def M_of_tau(tau_s):
    return 1e10 * (tau_s / 407.0) ** (1.0 / 3.0)


def tau_of_M(M_g):
    return 407.0 * (M_g / 1e10) ** 3


# ----------------------------------------------------------------------------
# 1. Photon spectra grid (S1 spectrum explorer)
# ----------------------------------------------------------------------------
def gen_spectra():
    masses = [1e17, 1e15, 1e14, 1e13, 1e12, 1e11, 1e10, 1e9]
    E = np.logspace(-8, 6, 240)
    out = {"energy_GeV": [float(f"{v:.6g}") for v in E], "masses": []}
    for M in masses:
        d = spec_dir(E, M)
        f = spec_frag(E, M)
        out["masses"].append(
            {
                "mass_g": M,
                "T_GeV": float(f"{T_of_M_GeV(M):.4g}"),
                "tau_s": float(f"{tau_of_M(M):.4g}"),
                "direct": [float(f"{v:.4g}") for v in d],
                "frag": [float(f"{v:.4g}") for v in f],
            }
        )
    write("spectra.json", out)


# ----------------------------------------------------------------------------
# 2. Effective areas (S4)
# ----------------------------------------------------------------------------
def load_effective_areas():
    areas = {}
    for path in glob(os.path.join(BHRAD, "Modelling", "EffectiveAreas", "*.dat")):
        name = os.path.basename(path)[:-4]
        if name == "propermotion":
            continue
        df = pd.read_csv(path, sep=" ", header=None, names=["E", "A"])
        areas[name] = df
    return areas


def gen_effective_areas(areas):
    out = {}
    for name, df in areas.items():
        out[name] = {
            "E_GeV": [float(f"{v:.5g}") for v in df["E"]],
            "Aeff_cm2": [float(f"{v:.5g}") for v in df["A"]],
        }
    write("effective-areas.json", out)


# ----------------------------------------------------------------------------
# 3. Detectability curves: d_max(M) per detector (S5; JCAP Fig. 4 analog)
#    N_S >= 10  OR  N_S / sqrt(N_B) >= 5, Tobs = 1 yr
# ----------------------------------------------------------------------------
ERANGES = {  # log10 GeV, from Analytical_Modelling.ipynb
    "LAT": (-2, 3),
    "GBM": (-4, -2),
    "BATSE": (-5, -2),
    "HAWC": (1, 2),
    "VERITAS": (1, 4),
    "LHAASO": (1, 6),
}
DOMEGA = {
    "LAT": 1e-3,
    "GBM": 1e-3,
    "BATSE": 1e-3,
    "HAWC": 1e-5,
    "VERITAS": 1e-5,
    "LHAASO": 1e-4,
}
CR_REJECT = {  # gamma/hadron separation rejection factor (notebook gamma_had_sep)
    "LAT": 0.0,
    "GBM": 0.0,
    "BATSE": 0.0,
    "HAWC": 1e-2,
    "VERITAS": 1e-1,
    "LHAASO": 1e-5,
}


def phi_bkg_egb(E):
    """Isotropic extragalactic background, JCAP Eq. 4.1 [cm-2 GeV-1 s-1 sr-1]."""
    return 1.4e-6 * E ** (-2.1)


def phi_bkg_total(E, frej):
    """EGB + misidentified cosmic rays with rejection frej (notebook CRbkg)."""
    egb = phi_bkg_egb(E)
    if frej <= 0:
        return egb
    cr = frej * 1.2 * 7900e-4 * E ** (-2.65)
    return np.maximum(cr, egb)


def signal_integral(Erange, tau, Aeff, Tobs=YEAR_S, n_t=80):
    """S = ∫dt ∫dE  d²N/dEdt(M(u)) Aeff(E)  over remaining-lifetime window.

    Remaining lifetime runs from tau down to max(tau - Tobs, eps).
    Returns S such that N_S = S / (4 pi d^2).
    """
    u_lo = max(tau - Tobs, 1e-4)
    if u_lo >= tau:
        return 0.0
    us = np.logspace(math.log10(u_lo), math.log10(tau), n_t)
    A = Aeff(Erange)
    per_t = np.array(
        [simpson(spec_total(Erange, M_of_tau(u)) * A, x=Erange) for u in us]
    )
    return float(simpson(per_t, x=us))


def gen_dmax_curves(areas):
    # tau up to 1e21 s covers masses up to ~1.7e16 g (JCAP Fig. 4 range)
    taus = np.logspace(-3, 21, 110)
    out = {"detectors": [], "notes": "d in pc; mass in g; Tobs=1yr; N_S>=10 or 5sigma"}
    for det, df in areas.items():
        if det not in ERANGES:
            continue
        lo, hi = ERANGES[det]
        Erange = np.logspace(lo, hi, 120)
        Aeff = interp1d(df["E"], df["A"], bounds_error=False, fill_value=0.0)
        dOm = DOMEGA[det]
        bkg_rate = float(
            simpson(phi_bkg_total(Erange, CR_REJECT[det]) * Aeff(Erange) * dOm, x=Erange)
        )  # photons / s
        d_pc, m_g = [], []
        for tau in taus:
            S = signal_integral(Erange, tau, Aeff)
            if S <= 0:
                d_pc.append(0.0)
                m_g.append(float(M_of_tau(tau)))
                continue
            NB = bkg_rate * min(YEAR_S, tau)
            d1 = math.sqrt(S / (4 * math.pi * 10.0))  # cm
            d2 = math.sqrt(S / (4 * math.pi * 5.0 * math.sqrt(max(NB, 1e-30)))) if NB > 0 else d1
            d = min(d1, d2) / PC_TO_CM
            d_pc.append(float(f"{d:.4g}"))
            m_g.append(float(f"{M_of_tau(tau):.4g}"))
        out["detectors"].append(
            {"name": det, "mass_g": m_g, "tau_s": [float(f"{t:.4g}") for t in taus], "dmax_pc": d_pc}
        )
    # Proper motion contours: theta = (180/pi) v min(tau,Tobs) / d  (JCAP Eq. 4.9)
    v_cm_s = 220e5  # 220 km/s
    pm = {"theta_deg": [0.1, 1.0, 10.0], "curves": []}
    for theta in pm["theta_deg"]:
        d_pm = [
            float(
                f"{(180.0/math.pi) * v_cm_s * min(t, YEAR_S) / (theta * PC_TO_CM):.4g}"
            )
            for t in taus
        ]
        pm["curves"].append(d_pm)
    pm["mass_g"] = [float(f"{M_of_tau(t):.4g}") for t in taus]
    out["proper_motion"] = pm
    write("dmax-curves.json", out)


# ----------------------------------------------------------------------------
# 4. Spectral index gamma(E0) vs mass (S3/S6; JCAP Eq. 4.8, Fig. 6)
# ----------------------------------------------------------------------------
def gen_spectral_index():
    masses = np.logspace(8, 17, 200)

    def gamma_at(E0, M):
        e1, e2 = E0 * 0.99, E0 * 1.01
        f1 = spec_total(np.array([e1]), M)[0]
        f2 = spec_total(np.array([e2]), M)[0]
        if f1 <= 0 or f2 <= 0:
            return None
        return (math.log(f1) - math.log(f2)) / (math.log(e2) - math.log(e1))

    out = {"mass_g": [float(f"{m:.4g}") for m in masses], "pivots": []}
    for E0, label in [(1.0, "1 GeV"), (1e-3, "1 MeV")]:
        g = [gamma_at(E0, m) for m in masses]
        g = [None if (v is None or not math.isfinite(v) or abs(v) > 50) else float(f"{v:.4g}") for v in g]
        out["pivots"].append({"E0_GeV": E0, "label": label, "gamma": g})
    write("spectral-index.json", out)


# ----------------------------------------------------------------------------
# 5. Mass functions (S2; JCAP Eqs. 3.9-3.14)
# ----------------------------------------------------------------------------
def gen_mass_functions():
    M = np.logspace(9, 19, 220)
    MU = 5.1e14

    def lognormal(Mstar, sigma):
        return np.exp(-np.log(M / Mstar) ** 2 / (2 * sigma**2)) / (
            math.sqrt(2 * math.pi) * sigma * M
        )

    def powerlaw(Mstar, gamma):
        psi = np.zeros_like(M)
        if gamma > 0:
            mask = M < Mstar
            psi[mask] = gamma / M[mask] * (M[mask] / Mstar) ** gamma
        else:
            mask = M > Mstar
            psi[mask] = -gamma / M[mask] * (M[mask] / Mstar) ** gamma
        return psi

    def critical(Mc):
        return M**2.85 * np.exp(-((M / Mc) ** 2.85))

    def norm(psi):
        s = simpson(psi, x=M)
        return psi / s if s > 0 else psi

    out = {
        "mass_g": [float(f"{m:.4g}") for m in M],
        "MU_g": MU,
        "families": {
            "lognormal": [
                {
                    "label": f"σ={s}",
                    "sigma": s,
                    "Mstar": 5.1e14,
                    "psi": [float(f"{v:.4g}") for v in norm(lognormal(5.1e14, s))],
                }
                for s in [0.1, 0.5, 1.0]
            ],
            "powerlaw": [
                {
                    "label": f"γ={g}",
                    "gamma": g,
                    "Mstar": 1e15,
                    "psi": [float(f"{v:.4g}") for v in norm(powerlaw(1e15, g))],
                }
                for g in [2.0, -2.0]
            ],
            "critical": [
                {
                    "label": "Mc=5.7e15 g",
                    "Mc": 5.7e15,
                    "psi": [float(f"{v:.4g}") for v in norm(critical(5.7e15))],
                }
            ],
        },
    }
    write("mass-functions.json", out)


# ----------------------------------------------------------------------------
# 6. alpha(M) degrees-of-freedom staircase (S3)
#    Pedagogical: Page coefficient grows as species with m < T_BH(M) unlock.
#    Weights from MacGibbon 1991 per-dof contributions (s=1/2: 4.197e-4/dof,
#    s=1: 1.023e-4/dof in units of 1e-3 * (1e10g)^3/s ... we use relative).
# ----------------------------------------------------------------------------
SM_SPECIES = [
    # (label, mass_GeV, dof-weighted alpha contribution, kind)
    ("photon", 0.0, 0.060, "boson"),
    ("neutrinos ×3", 0.0, 0.147, "fermion"),
    ("graviton", 0.0, 0.007, "boson"),
    ("e⁺e⁻", 5.11e-4, 0.142, "fermion"),
    ("μ⁺μ⁻", 0.1057, 0.142, "fermion"),
    ("u,d quarks + gluons (QCD)", 0.2, 0.815, "qcd"),
    ("s quark", 0.2, 0.10, "qcd"),
    ("τ, c", 1.4, 0.24, "fermion"),
    ("b quark", 4.2, 0.10, "qcd"),
    ("W,Z", 85.0, 0.12, "boson"),
    ("t, Higgs", 125.0, 0.12, "fermion"),
]


def gen_alpha_dof():
    masses = np.logspace(6, 18, 400)  # g
    T = 1.058e13 / masses  # GeV
    base = []
    for Tv in T:
        a = sum(w for (_, m, w, _) in SM_SPECIES if Tv > m)
        base.append(a)
    base = np.array(base)
    thresholds = [
        {"label": lbl, "mass_GeV": m, "M_bh_g": float(f"{1.058e13/max(m,1e-12):.3g}") if m > 0 else None}
        for (lbl, m, w, _) in SM_SPECIES
        if m > 0
    ]
    out = {
        "mass_g": [float(f"{m:.4g}") for m in masses],
        "T_GeV": [float(f"{t:.4g}") for t in T],
        "alpha_rel_SM": [float(f"{a:.4g}") for a in base],
        "alpha_dark_x2": [float(f"{a*2:.4g}") for a in base],
        "thresholds": thresholds,
        "note": "alpha in relative units normalized to full SM ≈ 2.0; dark curve doubles dof",
    }
    write("alpha-dof.json", out)


# ----------------------------------------------------------------------------
# 7. GRB candidates (S7/S8) from catalog search CSVs
# ----------------------------------------------------------------------------
def gen_candidates():
    path = os.path.join(BHRAD, "GBM_Catalog_Searching", "SGRB_search", "H>1_T90[0.2-5]_RS=0.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        h = r["hardness"]
        rows.append(
            {
                "name": str(r["name"]),
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "hardness": None if (isinstance(h, float) and (math.isinf(h) or math.isnan(h))) else float(f"{h:.4g}"),
                "lle_t90": float(f"{r['lle_t90']:.4g}"),
                "gbm_t90": float(f"{r['gbm_cat_t90']:.4g}"),
            }
        )
    write("candidates.json", {"selection": "H>1, T90∈[0.2,5]s, no redshift", "sources": rows})


# ----------------------------------------------------------------------------
# 8. LAT transients (S10) with fitted tau/distance from thesis pipeline
# ----------------------------------------------------------------------------
def gen_transients():
    path = os.path.join(BHRAD, "TransientLATSources", "TransientSources_fitted_params.csv")
    df = pd.read_csv(path)
    rows = []
    for _, r in df.iterrows():
        name = str(r["NAME"]).replace("b'", "").replace("'", "").strip()
        row = {
            "name": name,
            "ra": float(r["RA"]),
            "dec": float(r["DEC"]),
            "lii": float(r["LII"]),
            "bii": float(r["BII"]),
            "index": float(f"{r['PLAW_PHOTON_INDEX']:.4g}"),
            "index_err": float(f"{r['PLAW_PHOTON_INDEX_ERROR']:.4g}"),
        }
        d = r.get("distance (pc)")
        t0 = r.get("T0 (MJD)")
        if pd.notna(d) and pd.notna(t0):
            row["distance_pc"] = float(f"{d:.4g}")
            row["dist_min"] = float(f"{r['dist_min']:.4g}") if pd.notna(r["dist_min"]) else None
            row["dist_max"] = float(f"{r['dist_max']:.4g}") if pd.notna(r["dist_max"]) else None
            row["T0_mjd"] = float(f"{t0:.6g}")
            # remaining lifetime in seconds at detection: tau ~ T0(MJD)*86400 (per notebook convention)
            row["tau_s"] = float(f"{t0*86400.0:.4g}")
            row["mass_g"] = float(f"{M_of_tau(row['tau_s']):.4g}")
        rows.append(row)
    write("transients.json", {"catalog": "1FLT unassociated, thesis §3.1", "sources": rows})


# ----------------------------------------------------------------------------
# 9. Real GBM lightcurves (S6), downsampled + background-subtracted
# ----------------------------------------------------------------------------
def gen_grb_lightcurves():
    src_dir = os.path.join(BHRAD, "Lightcurve_Fitting", "~100ms_Source_Data")
    out = {"grbs": []}
    for path in sorted(glob(os.path.join(src_dir, "*resolution_0.1.csv"))):
        base = os.path.basename(path)
        name = base.split("_")[0]
        df = pd.read_csv(path, index_col=0)
        tcol = "time"
        det_cols = [c for c in df.columns if not c.endswith("_bkg") and c != tcol]
        sig = np.zeros(len(df))
        for c in det_cols:
            b = f"{c}_bkg"
            if b in df.columns:
                sig += df[c].to_numpy() - df[b].to_numpy()
            else:
                sig += df[c].to_numpy()
        t = df[tcol].to_numpy()
        # trim to interesting window and downsample to <=500 points
        step = max(1, len(t) // 500)
        t_d = t[::step]
        s_d = sig[::step]
        out["grbs"].append(
            {
                "name": name,
                "detectors": det_cols,
                "t_s": [float(f"{v:.4g}") for v in t_d],
                "counts": [float(f"{v:.5g}") for v in s_d],
            }
        )
    write("grb-lightcurves.json", out)


# ----------------------------------------------------------------------------
# 10. Telescope energy bands (S4 ruler)
# ----------------------------------------------------------------------------
def gen_bands():
    bands = [
        {"name": "Fermi GBM", "lo_GeV": 8e-6, "hi_GeV": 0.04, "kind": "space"},
        {"name": "BATSE", "lo_GeV": 2e-5, "hi_GeV": 0.002, "kind": "space", "historical": True},
        {"name": "Fermi LAT", "lo_GeV": 0.02, "hi_GeV": 300.0, "kind": "space"},
        {"name": "VERITAS", "lo_GeV": 85.0, "hi_GeV": 30000.0, "kind": "ground"},
        {"name": "HAWC", "lo_GeV": 300.0, "hi_GeV": 100000.0, "kind": "ground"},
        {"name": "LHAASO", "lo_GeV": 20.0, "hi_GeV": 1e6, "kind": "ground"},
    ]
    write("telescope-bands.json", {"bands": bands})


# ----------------------------------------------------------------------------
# 11. Precomputed hero lightcurve: total flux vs remaining lifetime (S1)
# ----------------------------------------------------------------------------
def gen_lightcurve():
    taus = np.logspace(-3, 9, 160)
    E = np.logspace(-6, 6, 200)
    flux = []
    for tau in taus:
        M = M_of_tau(tau)
        flux.append(float(f"{simpson(spec_total(E, M), x=E):.4g}"))
    # per-band lightcurves for GBM & LAT bands (flux ratio interactive, S3)
    E_gbm = np.logspace(-4, -2, 100)
    E_lat = np.logspace(-1, 3, 100)
    f_gbm = [float(f"{simpson(spec_total(E_gbm, M_of_tau(t)), x=E_gbm):.4g}") for t in taus]
    f_lat = [float(f"{simpson(spec_total(E_lat, M_of_tau(t)), x=E_lat):.4g}") for t in taus]
    write(
        "lightcurve.json",
        {
            "tau_s": [float(f"{t:.4g}") for t in taus],
            "mass_g": [float(f"{M_of_tau(t):.4g}") for t in taus],
            "total_ph_per_s": flux,
            "gbm_band_ph_per_s": f_gbm,
            "lat_band_ph_per_s": f_lat,
        },
    )


if __name__ == "__main__":
    print("Generating website data from BHRad repo …")
    areas = load_effective_areas()
    gen_spectra()
    gen_effective_areas(areas)
    gen_dmax_curves(areas)
    gen_spectral_index()
    gen_mass_functions()
    gen_alpha_dof()
    gen_candidates()
    gen_transients()
    gen_grb_lightcurves()
    gen_bands()
    gen_lightcurve()
    print("Done.")
