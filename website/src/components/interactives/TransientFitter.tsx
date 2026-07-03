"use client";

import { useMemo, useState } from "react";
import transients from "@/data/transients.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { transientFlux } from "@/lib/physics";
import { C } from "@/lib/colors";
import { duration, sci } from "@/lib/format";

/**
 * The long-game search: a steadily brightening source. Fit F(t) ∝ (τ−t)^−0.533
 * by hand and see how distance/τ trade off (JCAP Eq. 5.3, thesis §3.1).
 */
export function TransientFitter() {
  const [logTau, setLogTau] = useState(8.6);
  const [logD, setLogD] = useState(-2.55);

  const tau = Math.pow(10, logTau);
  const d = Math.pow(10, logD);

  // synthetic decade of observations of a rising source (mock 1FLT cadence)
  const tObs = useMemo(() => Array.from({ length: 60 }, (_, i) => 1e6 + (i * 3.15e8) / 59), []);
  const model = tObs.map((t) => transientFlux(t, tau, d) ?? null);

  // "truth": tau ≈ 11 yr so a decade of data climbs visibly (1FLT group-1-like)
  const truth = useMemo(() => tObs.map((t) => transientFlux(t, 3.5e8, 2.74e-3)), [tObs]);
  const noisy = useMemo(
    () =>
      truth.map((v, i) => {
        const r = Math.sin(i * 12.9898) * 43758.5453;
        const noise = 1 + 0.35 * ((r - Math.floor(r)) - 0.5);
        return v * noise;
      }),
    [truth]
  );

  const fitted = (transients.sources as { distance_pc?: number }[]).filter(
    (s) => s.distance_pc
  ).length;

  return (
    <div>
      <Chart
        xDomain={[1e6, 3.5e8]}
        yDomain={[3e-8, 3e-6]}
        xLabel="time since first detection [s]  (~decade of LAT data)"
        yLabel="flux > 0.1 GeV  [cm⁻² s⁻¹]"
      >
        {({ sx, sy }) => (
          <>
            {/* mock data points */}
            {tObs.map((t, i) =>
              t > 1e6 && noisy[i] ? (
                <circle
                  key={i}
                  cx={sx(t)}
                  cy={sy(noisy[i])}
                  r="2.5"
                  fill={C.inkDim}
                  opacity="0.7"
                />
              ) : null
            )}
            <path d={linePath(tObs, model, sx, sy)} fill="none" stroke={C.signal} strokeWidth="2" />
          </>
        )}
      </Chart>

      <Legend
        items={[
          { label: "mock unassociated transient (1FLT-like)", color: C.inkDim },
          { label: "your model: F ∝ d⁻² (τ−t)^−0.533", color: C.signal },
        ]}
      />

      <div className="mt-4 grid gap-x-8 gap-y-3 sm:grid-cols-2">
        <div className="flex items-center gap-3">
          <span className="w-40 font-mono text-[11px] text-ink-faint">remaining lifetime τ</span>
          <input
            type="range"
            min={8.52}
            max={11.5}
            step={0.02}
            value={logTau}
            onChange={(e) => setLogTau(Number(e.target.value))}
            className="flex-1"
            aria-label="Remaining lifetime"
          />
          <span className="w-20 text-right font-mono text-xs text-signal">{duration(tau)}</span>
        </div>
        <div className="flex items-center gap-3">
          <span className="w-40 font-mono text-[11px] text-ink-faint">distance d</span>
          <input
            type="range"
            min={-4.5}
            max={-1}
            step={0.02}
            value={logD}
            onChange={(e) => setLogD(Number(e.target.value))}
            className="flex-1"
            aria-label="Distance"
          />
          <span className="w-20 text-right font-mono text-xs text-signal">{sci(d * 1e3)} mpc</span>
        </div>
      </div>

      <p className="mt-4 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Notice the degeneracy: raising τ and shrinking d can produce nearly identical decade-long
        light curves — exactly why the {fitted} fitted 1FLT sources carry wide error bars
        (TransientSources_fitted_params.csv). The tell that finally breaks the degeneracy is
        proper motion: at milliparsec distances the source should drift degrees per year, which
        the LAT catalog does not observe for these sources.
      </p>
    </div>
  );
}
