"use client";

import { useState } from "react";
import mf from "@/data/mass-functions.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { C } from "@/lib/colors";

const FAMILY_META: Record<
  string,
  { label: string; color: string; blurb: string }
> = {
  lognormal: {
    label: "Log-normal",
    color: C.signal,
    blurb:
      "Expected from a smooth, symmetric peak in the inflationary power spectrum. Width σ controls how much of the distribution leaks into constrained regions.",
  },
  powerlaw: {
    label: "Power-law",
    color: C.burst,
    blurb:
      "Arises from scale-invariant fluctuations or collapsing cosmic strings; the index γ is set by the equation of state w when the holes formed: γ = −2w/(1+w).",
  },
  critical: {
    label: "Critical collapse",
    color: C.theory,
    blurb:
      "ψ ∝ M².⁸⁵ exp(−(M/Mc)².⁸⁵) — the universal outcome of near-threshold collapse, extending to arbitrarily low masses below the horizon-mass cutoff.",
  },
};

export function MassFunctions() {
  const [family, setFamily] = useState<string>("lognormal");
  const M = mf.mass_g as number[];
  const fam = (mf.families as Record<string, { label: string; psi: number[] }[]>)[family];
  const meta = FAMILY_META[family];

  return (
    <div>
      <div className="mb-4 flex flex-wrap gap-2">
        {Object.keys(FAMILY_META).map((k) => (
          <button
            key={k}
            onClick={() => setFamily(k)}
            className={`rounded border px-3 py-1.5 font-mono text-xs transition-colors ${
              family === k
                ? "border-signal/60 bg-signal/10 text-signal"
                : "border-line text-ink-dim hover:border-line-bright"
            }`}
          >
            {FAMILY_META[k].label}
          </button>
        ))}
      </div>

      <Chart
        xDomain={[1e9, 1e19]}
        yDomain={[1e-22, 1e-12]}
        xLabel="PBH mass M [g]"
        yLabel="ψ(M) [normalized]"
      >
        {({ sx, sy }) => (
          <>
            {/* MU line: holes dying today */}
            <line
              x1={sx(mf.MU_g)}
              x2={sx(mf.MU_g)}
              y1={sy.range[0]}
              y2={sy.range[1]}
              stroke={C.danger}
              strokeDasharray="4 4"
              strokeWidth="1.2"
            />
            <text
              x={sx(mf.MU_g) + 5}
              y={sy.range[1] + 14}
              fontSize="10"
              fill={C.danger}
              fontFamily="var(--font-plex-mono)"
            >
              M_U — exploding now
            </text>
            {fam.map((curve, i) => (
              <path
                key={curve.label}
                d={linePath(M, curve.psi, sx, sy)}
                fill="none"
                stroke={meta.color}
                strokeWidth={1.8}
                strokeDasharray={i === 0 ? undefined : i === 1 ? "6 4" : "2 3"}
                opacity={1 - i * 0.15}
              />
            ))}
          </>
        )}
      </Chart>

      <Legend
        items={fam.map((c, i) => ({
          label: c.label,
          color: meta.color,
          dash: i === 0 ? undefined : i === 1 ? "6 4" : "2 3",
        }))}
      />
      <p className="mt-3 max-w-2xl text-[12px] leading-relaxed text-ink-faint">{meta.blurb}</p>
      <p className="mt-2 font-mono text-[11px] text-ink-faint">
        The explosion rate today is set by how much of ψ sits at M_U:{" "}
        <span className="text-ink-dim">ṅ = ρ_DM ψ(M_U) / 3t_U</span> (Boluna et al. Eq. 3.8). For
        every allowed mass function it stays far below the HAWC bound of 3400 pc⁻³ yr⁻¹ — unless
        the distribution is spiked almost exactly at M_U.
      </p>
    </div>
  );
}
