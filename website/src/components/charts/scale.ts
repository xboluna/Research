export type ScaleKind = "log" | "linear";

export interface Scale {
  (v: number): number;
  invert: (px: number) => number;
  kind: ScaleKind;
  domain: [number, number];
  range: [number, number];
}

export function makeScale(
  kind: ScaleKind,
  domain: [number, number],
  range: [number, number]
): Scale {
  const [d0, d1] = kind === "log" ? [Math.log10(domain[0]), Math.log10(domain[1])] : domain;
  const [r0, r1] = range;
  const f = ((v: number) => {
    const t = kind === "log" ? Math.log10(v) : v;
    return r0 + ((t - d0) / (d1 - d0)) * (r1 - r0);
  }) as Scale;
  f.invert = (px: number) => {
    const t = d0 + ((px - r0) / (r1 - r0)) * (d1 - d0);
    return kind === "log" ? Math.pow(10, t) : t;
  };
  f.kind = kind;
  f.domain = domain;
  f.range = range;
  return f;
}

/** Nice log ticks (powers of 10, decimated to <= maxTicks). */
export function logTicks(domain: [number, number], maxTicks = 8): number[] {
  const lo = Math.ceil(Math.log10(domain[0]));
  const hi = Math.floor(Math.log10(domain[1]));
  const all: number[] = [];
  for (let e = lo; e <= hi; e++) all.push(Math.pow(10, e));
  const step = Math.max(1, Math.ceil(all.length / maxTicks));
  return all.filter((_, i) => i % step === 0);
}

export function linTicks(domain: [number, number], count = 6): number[] {
  const [lo, hi] = domain;
  const span = hi - lo;
  const raw = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 5, 10].map((m) => m * mag).find((s) => span / s <= count) ?? mag * 10;
  const start = Math.ceil(lo / step) * step;
  const out: number[] = [];
  for (let v = start; v <= hi + step * 0.001; v += step) out.push(Number(v.toPrecision(10)));
  return out;
}

/** SVG path from x/y arrays through a pair of scales, skipping non-finite. */
export function linePath(
  xs: ArrayLike<number>,
  ys: ArrayLike<number | null>,
  sx: Scale,
  sy: Scale
): string {
  let d = "";
  let pen = false;
  for (let i = 0; i < xs.length; i++) {
    const yv = ys[i];
    const xv = xs[i];
    const ok =
      yv != null &&
      isFinite(yv) &&
      isFinite(xv) &&
      (sy.kind !== "log" || yv > 0) &&
      (sx.kind !== "log" || xv > 0);
    if (!ok) {
      pen = false;
      continue;
    }
    const X = sx(xv).toFixed(2);
    const Y = sy(yv).toFixed(2);
    d += pen ? `L${X},${Y}` : `M${X},${Y}`;
    pen = true;
  }
  return d;
}

const SUP: Record<string, string> = {
  "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴",
  "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻",
};

export function tickLabel(v: number, kind: ScaleKind): string {
  if (kind === "log") {
    const e = Math.round(Math.log10(v));
    if (e === 0) return "1";
    if (e === 1) return "10";
    return `10${String(e).split("").map((c) => SUP[c] ?? c).join("")}`;
  }
  if (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0)) return v.toExponential(0);
  return String(Number(v.toPrecision(4)));
}
