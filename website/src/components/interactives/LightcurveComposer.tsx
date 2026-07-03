"use client";

import { useMemo, useState } from "react";
import grbs from "@/data/grb-lightcurves.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { afterglow, burstLightcurve } from "@/lib/physics";
import { C } from "@/lib/colors";

/**
 * Fit-by-hand: real GBM counts for candidate GRBs with the PBH power-law +
 * afterglow template overlaid; user drives the parameters the sampler fits.
 */
export function LightcurveComposer() {
  const list = grbs.grbs as { name: string; t_s: number[]; counts: number[] }[];
  const [gi, setGi] = useState(5); // GRB141222298 default if present
  const [tau, setTau] = useState(1.5);
  const [norm, setNorm] = useState(3.2);
  const [agOn, setAgOn] = useState(true);
  const [agT, setAgT] = useState(12);
  const [agN, setAgN] = useState(2.0);

  const g = list[Math.min(gi, list.length - 1)];

  const { model, peak } = useMemo(() => {
    const normAbs = Math.pow(10, norm);
    const agAbs = agOn ? Math.pow(10, agN) * 100 : 0;
    const m = g.t_s.map((t) => {
      const b = burstLightcurve(t, tau, normAbs);
      const a = agOn ? afterglow(t, 2.5, tau + 1.5, tau + agT, agAbs) : 0;
      return b + a;
    });
    return { model: m, peak: Math.max(...g.counts) };
  }, [g, tau, norm, agOn, agT, agN]);

  const yMax = Math.max(peak * 1.6, 100);
  const yMin = Math.min(...g.counts) < 0 ? Math.min(...g.counts) * 1.3 : 0;

  return (
    <div>
      <div className="mb-4 flex flex-wrap items-center gap-2">
        {list.map((x, i) => (
          <button
            key={x.name}
            onClick={() => setGi(i)}
            className={`rounded border px-2.5 py-1 font-mono text-[11px] transition-colors ${
              i === gi
                ? "border-burst/60 bg-burst/10 text-burst"
                : "border-line text-ink-faint hover:border-line-bright"
            }`}
          >
            {x.name}
          </button>
        ))}
      </div>

      <Chart
        xDomain={[Math.min(...g.t_s), Math.max(...g.t_s)]}
        yDomain={[yMin, yMax]}
        xKind="linear"
        yKind="linear"
        xLabel="time since trigger [s]"
        yLabel="background-subtracted counts / 100 ms"
        height={360}
      >
        {({ sx, sy }) => (
          <>
            {/* data as dithered histogram strokes */}
            <path
              d={linePath(g.t_s, g.counts, sx, sy)}
              fill="none"
              stroke={C.inkDim}
              strokeWidth="1"
              opacity="0.8"
            />
            {/* model */}
            <path
              d={linePath(g.t_s, model, sx, sy)}
              fill="none"
              stroke={C.signal}
              strokeWidth="2"
            />
            {/* explosion epoch marker */}
            <line
              x1={sx(tau)}
              x2={sx(tau)}
              y1={sy.range[0]}
              y2={sy.range[1]}
              stroke={C.burst}
              strokeDasharray="4 3"
            />
            <text
              x={sx(tau) + 4}
              y={sy.range[1] + 12}
              fontSize="9"
              fill={C.burst}
              fontFamily="var(--font-plex-mono)"
            >
              τ (burst end)
            </text>
          </>
        )}
      </Chart>

      <Legend
        items={[
          { label: `${g.name} — real GBM data (100 ms bins)`, color: C.inkDim },
          { label: "PBH template: (τ−t)^−0.52 + afterglow", color: C.signal },
        ]}
      />

      <div className="mt-4 grid gap-x-8 gap-y-3 sm:grid-cols-2">
        <Param label="burst epoch τ" min={-5} max={30} step={0.1} v={tau} set={setTau} unit="s" />
        <Param label="log₁₀ norm" min={1} max={5} step={0.05} v={norm} set={setNorm} unit="" />
        <div className="flex items-center gap-2">
          <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
            <input
              type="checkbox"
              checked={agOn}
              onChange={(e) => setAgOn(e.target.checked)}
              className="accent-[#E8B45A]"
            />
            afterglow component
          </label>
        </div>
        {agOn && (
          <>
            <Param label="afterglow spread" min={2} max={40} step={0.5} v={agT} set={setAgT} unit="s" />
            <Param label="log₁₀ afterglow norm" min={0} max={4} step={0.05} v={agN} set={setAgN} unit="" />
          </>
        )}
      </div>

      <p className="mt-4 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        These are the actual light curves fitted in the thesis pipeline
        (Lightcurve_Fitting/Fitting.ipynb). In the real analysis, ultranest explores this exact
        parameter space against every GBM detector simultaneously via threeML — you are doing by
        hand what nested sampling does with 400 live points.
      </p>
    </div>
  );
}

function Param({
  label,
  min,
  max,
  step,
  v,
  set,
  unit,
}: {
  label: string;
  min: number;
  max: number;
  step: number;
  v: number;
  set: (n: number) => void;
  unit: string;
}) {
  return (
    <div className="flex items-center gap-3">
      <span className="w-40 font-mono text-[11px] text-ink-faint">{label}</span>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={v}
        onChange={(e) => set(Number(e.target.value))}
        className="flex-1"
        aria-label={label}
      />
      <span className="w-16 text-right font-mono text-xs text-signal">
        {v.toFixed(step < 0.1 ? 2 : 1)}
        {unit}
      </span>
    </div>
  );
}
