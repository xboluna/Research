"use client";

import si from "@/data/spectral-index.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { C } from "@/lib/colors";

/** JCAP Fig. 6 — spectral index vs mass, converging to universal γ→1.5. */
export function SpectralIndex() {
  const M = si.mass_g as number[];
  const pivots = si.pivots as { label: string; gamma: (number | null)[] }[];

  return (
    <div>
      <Chart
        xDomain={[1e8, 1e17]}
        yDomain={[-1, 6]}
        yKind="linear"
        xLabel="black hole mass [g]   (evaporation runs right → left)"
        yLabel="spectral index γ  (dN/dE ∝ E^−γ)"
      >
        {({ sx, sy }) => (
          <>
            {/* universal asymptote */}
            <line
              x1={sx.range[0]}
              x2={sx.range[1]}
              y1={sy(1.5)}
              y2={sy(1.5)}
              stroke={C.danger}
              strokeDasharray="5 4"
              strokeWidth="1.2"
            />
            <text
              x={sx.range[0] + 8}
              y={sy(1.5) - 6}
              fontSize="10"
              fill={C.danger}
              fontFamily="var(--font-plex-mono)"
            >
              γ → 1.5 universal endpoint
            </text>
            <path
              d={linePath(M, pivots[0].gamma, sx, sy)}
              fill="none"
              stroke={C.signal}
              strokeWidth="2"
            />
            <path
              d={linePath(M, pivots[1].gamma, sx, sy)}
              fill="none"
              stroke={C.burst}
              strokeWidth="2"
              strokeDasharray="6 4"
            />
          </>
        )}
      </Chart>
      <Legend
        items={[
          { label: "pivot 1 GeV", color: C.signal },
          { label: "pivot 1 MeV", color: C.burst, dash: "6 4" },
        ]}
      />
      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        As the hole approaches its final moments, the measured power-law slope at any energy
        converges to γ = 1.5 — a fingerprint no astrophysical source shares. This is the primary
        template used to screen the Fermi transient catalog. (Boluna et al. Eq. 4.8, Fig. 6.)
      </p>
    </div>
  );
}
