"use client";

import { useState } from "react";
import dof from "@/data/alpha-dof.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { C } from "@/lib/colors";

/**
 * The Page coefficient α(M): every particle species lighter than T_BH adds a
 * channel. Toggle a hypothetical dark sector to see evaporation accelerate.
 */
export function AlphaStaircase() {
  const [dark, setDark] = useState(false);

  const M = dof.mass_g as number[];
  const T = dof.T_GeV as number[];
  const aSM = dof.alpha_rel_SM as number[];
  const aDark = dof.alpha_dark_x2 as number[];

  return (
    <div>
      <div className="mb-4 flex items-center justify-between gap-4">
        <label className="flex cursor-pointer items-center gap-2.5 font-mono text-xs text-ink-dim">
          <input
            type="checkbox"
            checked={dark}
            onChange={(e) => setDark(e.target.checked)}
            className="accent-[#9D8CFF]"
          />
          add a dark sector (×2 degrees of freedom)
        </label>
        {dark && (
          <span className="font-mono text-[11px] text-theory">
            lifetime shortens: τ ∝ 1/α — the hole dies faster, dimmer in photons
          </span>
        )}
      </div>

      <Chart
        xDomain={[1e6, 1e18]}
        yDomain={[0, 4.6]}
        yKind="linear"
        xLabel="black hole mass M [g]   (← hotter · lighter)"
        yLabel="α(M) — evaporation channels [rel.]"
      >
        {({ sx, sy }) => (
          <>
            {/* threshold markers */}
            {(dof.thresholds as { label: string; M_bh_g: number | null }[]).map(
              (t) =>
                t.M_bh_g &&
                t.M_bh_g > 1e6 &&
                t.M_bh_g < 1e18 && (
                  <g key={t.label}>
                    <line
                      x1={sx(t.M_bh_g)}
                      x2={sx(t.M_bh_g)}
                      y1={sy.range[0]}
                      y2={sy.range[1]}
                      stroke={C.line}
                    />
                    <text
                      x={sx(t.M_bh_g)}
                      y={sy.range[1] + 10}
                      fontSize="8.5"
                      fill={C.inkFaint}
                      fontFamily="var(--font-plex-mono)"
                      transform={`rotate(-38 ${sx(t.M_bh_g)} ${sy.range[1] + 10})`}
                    >
                      {t.label}
                    </text>
                  </g>
                )
            )}
            <path d={linePath(M, aSM, sx, sy)} fill="none" stroke={C.signal} strokeWidth="2" />
            {dark && (
              <path
                d={linePath(M, aDark, sx, sy)}
                fill="none"
                stroke={C.theory}
                strokeWidth="2"
                strokeDasharray="6 4"
              />
            )}
          </>
        )}
      </Chart>

      <Legend
        items={[
          { label: "Standard Model", color: C.signal },
          ...(dark ? [{ label: "SM + dark sector", color: C.theory, dash: "6 4" }] : []),
        ]}
      />

      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Each step: the shrinking hole gets hot enough to radiate a new species (T = 1.058×10¹³
        GeV·g / M). dM/dt = −α(M)/M², so more channels → faster evaporation. Dark degrees of
        freedom leave gravity no place to hide — they must be radiated too, whether or not they
        couple to light. Thesis §2.8.
      </p>
    </div>
  );
}
