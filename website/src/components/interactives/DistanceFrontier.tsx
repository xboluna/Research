"use client";

import { useState } from "react";
import dmax from "@/data/dmax-curves.json";
import transients from "@/data/transients.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { RepoLink } from "@/components/ui/RepoLink";
import { REPO_PATHS } from "@/lib/repo";
import { DETECTOR_COLORS, C } from "@/lib/colors";

/**
 * The sensitivity frontier: max distance vs BH mass per detector
 * (JCAP Fig. 4 analog, regenerated from the repo's detectability pipeline),
 * with proper-motion contours and fitted LAT transients overlaid.
 */
export function DistanceFrontier() {
  const [showPM, setShowPM] = useState(true);
  const [showTransients, setShowTransients] = useState(false);

  const dets = dmax.detectors as {
    name: string;
    mass_g: number[];
    dmax_pc: number[];
  }[];
  const pm = dmax.proper_motion as {
    theta_deg: number[];
    mass_g: number[];
    curves: number[][];
  };
  const tr = (transients.sources as {
    name: string;
    mass_g?: number;
    distance_pc?: number;
  }[]).filter((s) => s.mass_g && s.distance_pc);

  return (
    <div>
      <div className="mb-4 flex flex-wrap gap-x-6 gap-y-2">
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={showPM}
            onChange={(e) => setShowPM(e.target.checked)}
            className="accent-[#8B91A7]"
          />
          proper-motion contours (0.1° · 1° · 10° per obs.)
        </label>
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={showTransients}
            onChange={(e) => setShowTransients(e.target.checked)}
            className="accent-[#E86A6A]"
          />
          overlay fitted 1FLT transients (§10)
        </label>
      </div>

      <Chart
        xDomain={[1e9, 1e16]}
        yDomain={[1e-7, 10]}
        xLabel="black hole mass [g]"
        yLabel="max visible distance [pc]"
        height={420}
      >
        {({ sx, sy }) => (
          <>
            {showPM &&
              pm.theta_deg.map((th, i) => (
                <g key={th}>
                  <path
                    d={linePath(pm.mass_g, pm.curves[i], sx, sy)}
                    fill="none"
                    stroke={C.inkFaint}
                    strokeWidth="1"
                    strokeDasharray="3 4"
                  />
                  <text
                    x={sx.range[1] - 34}
                    y={sy(pm.curves[i][pm.curves[i].length - 1]) - 4}
                    fontSize="9"
                    fill={C.inkFaint}
                    fontFamily="var(--font-plex-mono)"
                  >
                    {th}°
                  </text>
                </g>
              ))}
            {dets.map((d) => (
              <path
                key={d.name}
                d={linePath(d.mass_g, d.dmax_pc, sx, sy)}
                fill="none"
                stroke={DETECTOR_COLORS[d.name] ?? C.ink}
                strokeWidth="2"
              />
            ))}
            {showTransients &&
              tr.map((s) => (
                <g key={s.name}>
                  <circle
                    cx={sx(s.mass_g!)}
                    cy={sy(s.distance_pc!)}
                    r="3.5"
                    fill="none"
                    stroke={C.danger}
                    strokeWidth="1.5"
                  />
                  <line
                    x1={sx(s.mass_g!) - 3.5}
                    y1={sy(s.distance_pc!)}
                    x2={sx(s.mass_g!) + 3.5}
                    y2={sy(s.distance_pc!)}
                    stroke={C.danger}
                    strokeWidth="1"
                  />
                </g>
              ))}
          </>
        )}
      </Chart>

      <Legend
        items={[
          ...dets.map((d) => ({
            label: d.name,
            color: DETECTOR_COLORS[d.name] ?? C.ink,
          })),
          ...(showPM
            ? [{ label: "proper motion threshold", color: C.inkFaint, dash: "3 4" }]
            : []),
          ...(showTransients
            ? [{ label: "fitted 1FLT transients", color: C.danger }]
            : []),
        ]}
      />

      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Anything below a curve is detectable by that instrument (N_S ≥ 10 photons or 5σ over
        background, 1-yr observation). Regenerated from the{" "}
        <RepoLink path={REPO_PATHS.analyticalModelling} file>
          repo detectability pipeline
        </RepoLink>{" "}
        (Analytical_Modelling.ipynb; cf. Boluna et al. Fig. 4). Above the dashed contours the
        source would visibly streak across the sky — a smoking-gun signature, but also a
        challenge for catalog association.
      </p>
    </div>
  );
}
