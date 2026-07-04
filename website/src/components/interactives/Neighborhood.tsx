"use client";

import { useState } from "react";
import { C } from "@/lib/colors";
import { sci } from "@/lib/format";

/**
 * How close is "a few parsecs"? Log-scale map from the Sun outward with the
 * detection horizons drawn as dithered shells.
 */
const LANDMARKS = [
  { d: 4.85e-6, label: "Earth (1 AU)", dim: true },
  { d: 1.45e-4, label: "Voyager 1 (~165 AU)", dim: true },
  { d: 0.0194, label: "Oort cloud edge (~4000 AU)", dim: true },
  { d: 1.3, label: "Proxima Centauri", dim: false },
  { d: 3.6, label: "Barnard's Star", dim: true },
  { d: 12, label: "~30 nearest stars", dim: true },
];

const HORIZONS = [
  { d: 0.08, label: "LAT horizon (peak)", color: C.signal },
  { d: 0.1, label: "HAWC horizon (peak)", color: C.theory },
  { d: 1.3, label: "LHAASO horizon (peak)", color: "#6AAEE8" },
];

const D_MIN = 1e-6;
const D_MAX = 30;

function toX(d: number): number {
  return (Math.log10(d) - Math.log10(D_MIN)) / (Math.log10(D_MAX) - Math.log10(D_MIN));
}

export function Neighborhood() {
  const [probe, setProbe] = useState(0.1);

  // expected number of PBHs within probe distance if they saturate f_PBH at M_U:
  // n ~ f · rho_DM / M. rho_DM = 0.4 GeV/cm3 -> in g/pc3: 0.4*1.78e-24 g/cm3 * (3.086e18)^3 cm3/pc3
  const rhoDM_g_pc3 = 0.4 * 1.78e-24 * Math.pow(3.086e18, 3);
  const fmax = 1e-8; // approximate EGB bound at M_U
  const nPBH = (fmax * rhoDM_g_pc3) / 5.1e14; // per pc^3
  const expected = nPBH * (4 / 3) * Math.PI * Math.pow(probe, 3);

  return (
    <div>
      <div className="relative h-40 overflow-hidden rounded border border-line bg-void">
        {/* dithered distance shells */}
        <div className="dither-25 absolute inset-0" />
        {/* sun */}
        <div className="absolute left-[2%] top-1/2 -translate-y-1/2">
          <div className="h-3 w-3 rounded-full bg-burst shadow-[0_0_12px_2px_rgba(232,180,90,0.5)]" />
          <div className="mt-1 font-mono text-[9px] text-burst">Sun</div>
        </div>
        {/* landmarks */}
        {LANDMARKS.map((l) => (
          <div
            key={l.label}
            className="absolute top-0 h-full"
            style={{ left: `${2 + toX(l.d) * 95}%` }}
          >
            <div className={`h-full w-px ${l.dim ? "bg-line" : "bg-line-bright"}`} />
            <div
              className={`absolute top-2 -translate-x-1/2 whitespace-nowrap font-mono text-[8.5px] ${
                l.dim ? "text-ink-faint" : "text-ink-dim"
              }`}
              style={{ writingMode: "vertical-rl" }}
            >
              {l.label}
            </div>
          </div>
        ))}
        {/* horizons */}
        {HORIZONS.map((h) => (
          <div
            key={h.label}
            className="absolute top-0 h-full"
            style={{ left: `${2 + toX(h.d) * 95}%` }}
          >
            <div className="h-full w-0.5" style={{ background: h.color, opacity: 0.75 }} />
            <div
              className="absolute bottom-1 -translate-x-1/2 whitespace-nowrap rounded px-1 font-mono text-[8.5px]"
              style={{ color: h.color }}
            >
              ◆
            </div>
          </div>
        ))}
        {/* probe */}
        <div
          className="absolute top-0 h-full w-px bg-danger"
          style={{ left: `${2 + toX(probe) * 95}%` }}
        />
      </div>

      <div className="mt-2 flex justify-between font-mono text-[9px] text-ink-faint">
        <span>1 AU</span>
        <span>10⁻⁴ pc</span>
        <span>10⁻² pc</span>
        <span>1 pc</span>
        <span>30 pc</span>
      </div>

      <div className="mt-4 flex items-center gap-3">
        <span className="font-mono text-[11px] uppercase tracking-widest text-ink-faint">
          search radius
        </span>
        <input
          type="range"
          min={Math.log10(1e-4)}
          max={Math.log10(20)}
          step={0.02}
          value={Math.log10(probe)}
          onChange={(e) => setProbe(Math.pow(10, Number(e.target.value)))}
          className="flex-1"
          aria-label="Search radius in parsecs"
        />
        <span className="w-24 text-right font-mono text-sm text-danger">
          {probe < 0.01 ? `${sci(probe * 206265)} AU` : `${sci(probe)} pc`}
        </span>
      </div>

      <div className="mt-3 grid gap-3 font-mono text-xs sm:grid-cols-2">
        <div className="rounded border border-line bg-void px-3 py-2.5">
          <span className="text-ink-faint">volume · </span>
          <span className="text-ink">{sci((4 / 3) * Math.PI * Math.pow(probe, 3))} pc³</span>
        </div>
        <div className="rounded border border-line bg-void px-3 py-2.5">
          <span className="text-ink-faint">expected PBHs near explosion (EGB-limited) · </span>
          <span className={expected > 1 ? "text-signal" : "text-ink"}>{sci(expected)}</span>
        </div>
      </div>
      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Detection horizons (◆) are the peak reach from the sensitivity frontier above. The tension
        in one picture: γ-ray telescopes only see explosions within a bubble much smaller than the
        distance to the nearest star, while abundance limits make such nearby events rare.
      </p>
    </div>
  );
}
