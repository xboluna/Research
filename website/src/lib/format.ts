/** Formatting helpers for scientific display. */

const SUP: Record<string, string> = {
  "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
  "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻",
};

export function sup(n: number | string): string {
  return String(n)
    .split("")
    .map((c) => SUP[c] ?? c)
    .join("");
}

/** 3.2×10⁻⁵ style formatting. */
export function sci(v: number, digits = 2): string {
  if (v === 0) return "0";
  if (!isFinite(v)) return "—";
  const exp = Math.floor(Math.log10(Math.abs(v)));
  if (exp >= -2 && exp <= 3) {
    return v.toPrecision(digits + 1).replace(/\.?0+$/, "");
  }
  const mant = v / Math.pow(10, exp);
  const m = mant.toFixed(digits - 1).replace(/\.?0+$/, "");
  return `${m === "1" ? "" : m + "×"}10${sup(exp)}`;
}

/** Human-friendly duration from seconds. */
export function duration(s: number): string {
  if (!isFinite(s)) return "—";
  if (s < 1e-3) return `${sci(s * 1e6)} µs`;
  if (s < 1) return `${sci(s * 1e3)} ms`;
  if (s < 120) return `${sci(s)} s`;
  if (s < 7200) return `${sci(s / 60)} min`;
  if (s < 172800) return `${sci(s / 3600)} hr`;
  if (s < 6.31e7) return `${sci(s / 86400)} days`;
  return `${sci(s / 3.156e7)} yr`;
}

/** Energy label from GeV. */
export function energy(gev: number): string {
  if (gev >= 1e6) return `${sci(gev / 1e6)} PeV`;
  if (gev >= 1e3) return `${sci(gev / 1e3)} TeV`;
  if (gev >= 1) return `${sci(gev)} GeV`;
  if (gev >= 1e-3) return `${sci(gev * 1e3)} MeV`;
  if (gev >= 1e-6) return `${sci(gev * 1e6)} keV`;
  return `${sci(gev * 1e9)} eV`;
}

export function mass(g: number): string {
  return `${sci(g)} g`;
}

export function distance(pc: number): string {
  if (pc >= 1e-2) return `${sci(pc)} pc`;
  if (pc >= 1e-5) return `${sci(pc * 1e3)} mpc`;
  const au = pc * 206265;
  return `${sci(au)} AU`;
}
