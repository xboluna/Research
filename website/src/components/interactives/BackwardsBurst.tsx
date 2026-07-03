"use client";

import { useMemo, useState } from "react";
import lc from "@/data/lightcurve.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { C } from "@/lib/colors";
import { duration, mass, sci } from "@/lib/format";

/**
 * The universal PBH lightcurve: photon output vs remaining lifetime.
 * Time axis is reversed (remaining lifetime shrinking → explosion at right).
 */
export function BackwardsBurst() {
  const [cursor, setCursor] = useState(0.35); // fraction along tau axis (0 = long lifetime)

  const taus = lc.tau_s as number[];
  const flux = lc.total_ph_per_s as number[];
  const masses = lc.mass_g as number[];

  // reversed x: plot vs tau but with domain flipped via negative mapping
  const idx = Math.min(
    taus.length - 1,
    Math.max(0, Math.round((1 - cursor) * (taus.length - 1)))
  );

  const revX = useMemo(() => taus.map((t) => 1 / t), [taus]); // 1/tau grows toward explosion

  return (
    <div>
      <Chart
        xDomain={[1e-9, 1e3]}
        yDomain={[1e27, 1e35]}
        xLabel="1 / remaining lifetime  [s⁻¹]  →  explosion"
        yLabel="photon output [s⁻¹]"
      >
        {({ sx, sy }) => (
          <>
            <path d={linePath(revX, flux, sx, sy)} fill="none" stroke={C.signal} strokeWidth="2" />
            {/* cursor */}
            <line
              x1={sx(revX[idx])}
              x2={sx(revX[idx])}
              y1={sy.range[0]}
              y2={sy.range[1]}
              stroke={C.burst}
              strokeWidth="1.5"
              strokeDasharray="4 3"
            />
            <circle cx={sx(revX[idx])} cy={sy(flux[idx])} r="4" fill={C.burst} />
          </>
        )}
      </Chart>

      <div className="mt-3 flex items-center gap-3">
        <span className="font-mono text-[11px] uppercase tracking-widest text-ink-faint">
          scrub
        </span>
        <input
          type="range"
          min={0}
          max={1}
          step={0.002}
          value={cursor}
          onChange={(e) => setCursor(Number(e.target.value))}
          className="flex-1"
          aria-label="Scrub along the lightcurve"
        />
      </div>

      <div className="mt-3 grid gap-3 font-mono text-xs sm:grid-cols-3">
        <div className="rounded border border-line bg-void px-3 py-2">
          <span className="text-ink-faint">remaining lifetime · </span>
          <span className="text-ink">{duration(taus[idx])}</span>
        </div>
        <div className="rounded border border-line bg-void px-3 py-2">
          <span className="text-ink-faint">mass left · </span>
          <span className="text-ink">{mass(masses[idx])}</span>
        </div>
        <div className="rounded border border-line bg-void px-3 py-2">
          <span className="text-ink-faint">output · </span>
          <span className="text-signal">{sci(flux[idx])} γ/s</span>
        </div>
      </div>

      <Legend items={[{ label: "universal lightcurve — it only ever brightens", color: C.signal }]} />
    </div>
  );
}
