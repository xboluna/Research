"use client";

import { useMemo, useState } from "react";
import spectra from "@/data/spectra.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { C } from "@/lib/colors";
import { duration, energy, mass } from "@/lib/format";

export function SpectrumExplorer() {
  const [mi, setMi] = useState(2); // index into spectra.masses (1e14 g default)
  const [showComponents, setShowComponents] = useState(true);

  const E = spectra.energy_GeV;
  const entry = spectra.masses[mi];

  const total = useMemo(
    () => entry.direct.map((d: number, i: number) => d + entry.frag[i]),
    [entry]
  );

  return (
    <div>
      <div className="mb-4 flex flex-wrap items-center gap-x-6 gap-y-3">
        <div className="flex items-center gap-3">
          <span className="font-mono text-[11px] uppercase tracking-widest text-ink-faint">
            mass
          </span>
          <input
            type="range"
            min={0}
            max={spectra.masses.length - 1}
            step={1}
            value={mi}
            onChange={(e) => setMi(Number(e.target.value))}
            className="w-44"
            aria-label="Select black hole mass"
          />
          <span className="font-mono text-sm text-signal">{mass(entry.mass_g)}</span>
        </div>
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={showComponents}
            onChange={(e) => setShowComponents(e.target.checked)}
            className="accent-[#7FD4C1]"
          />
          show direct / fragmentation split
        </label>
      </div>

      <Chart
        xDomain={[1e-7, 1e5]}
        yDomain={[1e-2, 1e28]}
        xLabel="photon energy E [GeV]"
        yLabel="d²N/dE dt  [GeV⁻¹ s⁻¹]"
      >
        {({ sx, sy }) => (
          <>
            {showComponents && (
              <>
                <path
                  d={linePath(E, entry.direct, sx, sy)}
                  fill="none"
                  stroke={C.theory}
                  strokeWidth="1.5"
                  strokeDasharray="5 4"
                />
                <path
                  d={linePath(E, entry.frag, sx, sy)}
                  fill="none"
                  stroke={C.burst}
                  strokeWidth="1.5"
                  strokeDasharray="2 3"
                />
              </>
            )}
            <path d={linePath(E, total, sx, sy)} fill="none" stroke={C.signal} strokeWidth="2" />
            {/* temperature marker */}
            <line
              x1={sx(entry.T_GeV)}
              x2={sx(entry.T_GeV)}
              y1={sy.range[0]}
              y2={sy.range[1]}
              stroke={C.inkFaint}
              strokeDasharray="3 4"
            />
            <text
              x={sx(entry.T_GeV) + 5}
              y={sy.range[1] + 14}
              fontSize="10"
              fill={C.inkDim}
              fontFamily="var(--font-plex-mono)"
            >
              k_BT = {energy(entry.T_GeV)}
            </text>
          </>
        )}
      </Chart>

      <Legend
        items={[
          { label: "total", color: C.signal },
          ...(showComponents
            ? [
                { label: "direct (primary Hawking)", color: C.theory, dash: "5 4" },
                { label: "fragmentation (quark/gluon jets → π⁰ → γγ)", color: C.burst, dash: "2 3" },
              ]
            : []),
        ]}
      />

      <p className="mt-3 font-mono text-[11px] leading-relaxed text-ink-faint">
        T = {energy(entry.T_GeV)} · remaining lifetime ≈ {duration(entry.tau_s)}. Parameterization
        from Ukwatta et al. (2016) Eqs. 31–34, as implemented in{" "}
        <span className="text-ink-dim">Analytical_Modelling.ipynb</span>.
      </p>
    </div>
  );
}
