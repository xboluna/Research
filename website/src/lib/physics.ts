/**
 * Physics helpers ported from BHRad/Modelling/Analytical_Modelling.ipynb
 * and Boluna et al. arXiv:2307.06467 (JCAP).
 *
 * Units: mass in grams, energy in GeV, time in seconds unless noted.
 */

/** Hawking temperature in GeV for a BH of mass M grams (notebook Tbh_mass). */
export function temperatureGeV(massG: number): number {
  return 1.058e13 / massG;
}

/** Remaining lifetime in seconds for mass M grams (JCAP Eq. 4.4). */
export function lifetimeS(massG: number): number {
  return 407 * Math.pow(massG / 1e10, 3);
}

/** Mass in grams for a remaining lifetime tau seconds. */
export function massOfTau(tauS: number): number {
  return 1e10 * Math.pow(tauS / 407, 1 / 3);
}

const A_FRAG = 6.339e23;
const B_FRAG = 1.1367e24;

function xGamma(E: number, M: number): number {
  return E / ((1058 * 1e10) / M);
}

function thetaS(u: number): number {
  return 0.5 * (1 + Math.tanh(10 * u));
}

function Ffunc(y: number): number {
  if (y <= 2) return 1;
  const ly = Math.log(y);
  return Math.exp(-0.0962 - 1.982 * (ly - 1.908) * (1 + Math.tanh(20 * (ly - 1.908))));
}

/** Direct (primary) Hawking photon spectrum d²N/dEdt [GeV⁻¹ s⁻¹]. */
export function specDirect(E: number, M: number): number {
  const x = xGamma(E, M);
  if (x > 500) return 0;
  const denom = (Math.exp(x) - 1) * Ffunc(x);
  if (!isFinite(denom) || denom <= 0) return 0;
  return (1.13e19 * Math.pow(x, 6)) / denom;
}

/** Fragmentation (secondary) photon spectrum d²N/dEdt [GeV⁻¹ s⁻¹]. */
export function specFrag(E: number, M: number): number {
  const x = xGamma(E, M);
  if (x <= 0) return 0;
  const soft = A_FRAG * Math.pow(x, -1.5) * (1 - thetaS(x - 0.3));
  const hard = x < 500 ? (B_FRAG * Math.exp(-x) * thetaS(x - 0.3)) / (x * (x + 1)) : 0;
  return soft + hard;
}

export function specTotal(E: number, M: number): number {
  return specDirect(E, M) + specFrag(E, M);
}

/**
 * Long-term transient flux model, JCAP Eq. (5.3):
 * F(t) = 2.7e-8 (pc/d)^2 ((tau - t)/s)^(-0.533) cm^-2 s^-1  (E > 0.1 GeV)
 */
export function transientFlux(tS: number, tauS: number, dPc: number): number {
  const rem = tauS - tS;
  if (rem <= 0) return NaN;
  return 2.7e-8 * Math.pow(1 / dPc, 2) * Math.pow(rem, -0.533);
}

/** Proper motion in degrees over min(tau, Tobs), JCAP Eq. (4.9). ψ=90°. */
export function properMotionDeg(tauS: number, dPc: number, TobsS = 3.156e7): number {
  const v = 220e5; // cm/s
  const pcCm = 3.0857e18;
  return ((180 / Math.PI) * (v * Math.min(tauS, TobsS))) / (dPc * pcCm);
}

/**
 * PBH lightcurve power-law used in fitting (thesis §2.1, Fitting.ipynb):
 * counts ∝ (tau - t)^(-0.52) before explosion.
 */
export function burstLightcurve(tS: number, tauS: number, norm: number, index = -0.52): number {
  const rem = tauS - tS;
  if (rem <= 0) return 0;
  return norm * Math.pow(rem, index);
}

/**
 * Afterglow component from adiabatic blob expansion (Fitting.ipynb).
 * Rises after t_m, peaks at t_p (value = norm), then decays; delta sets
 * how sharply the pulse rises and falls (gamma-distribution-like shape).
 */
export function afterglow(
  tS: number,
  delta: number,
  tM: number,
  tP: number,
  norm: number
): number {
  if (tS <= tM) return 0;
  const u = (tS - tM) / Math.max(tP - tM, 1e-9);
  // peaks at u = 1 with value norm
  return norm * Math.pow(u, delta) * Math.exp(-delta * (u - 1));
}

/** Horizon mass at cosmic time t seconds: M_H ≈ 1e15 (t / 1e-23 s) g. */
export function horizonMassG(tS: number): number {
  return 1e15 * (tS / 1e-23);
}

export const M_EVAPORATING_NOW = 5.1e14; // grams (JCAP §3, M_U)
