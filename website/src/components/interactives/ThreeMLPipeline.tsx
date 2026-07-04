"use client";

import { useState } from "react";

const STAGES = [
  {
    id: "data",
    label: "Mission data",
    mono: "FermiGBMBurstCatalog · TimeSeriesBuilder",
    detail:
      "Raw time-tagged photon events from each GBM scintillator and LAT tracker are pulled straight from the Fermi archives. Every detector has its own response, background, and quirks.",
  },
  {
    id: "plugins",
    label: "Instrument plugins",
    mono: "DispersionSpectrumLike · FermiLATLike",
    detail:
      "threeML wraps each instrument in a plugin that knows its response matrix. This is what makes joint fitting honest: each dataset keeps its own likelihood, no lossy conversions.",
  },
  {
    id: "model",
    label: "One shared model",
    mono: "powerlaw(τ−t)^−0.52 + afterglow(t)",
    detail:
      "A single astromodels source — the PBH lightcurve plus an optional afterglow — is fit simultaneously against every instrument. Same parameters, all wavelengths.",
  },
  {
    id: "sample",
    label: "Bayesian sampling",
    mono: "ultranest.ReactiveNestedSampler",
    detail:
      "Nested sampling maps the full posterior over normalization, decay index, and afterglow timing — with honest uncertainties and Bayesian evidence for model comparison.",
  },
  {
    id: "verdict",
    label: "Verdict per source",
    mono: "logZ · posterior corners",
    detail:
      "Does the PBH template beat an ordinary GRB template? For candidates in this work, evidence ratios and residuals decide which sources stay on the list.",
  },
];

export function ThreeMLPipeline() {
  const [active, setActive] = useState(0);

  return (
    <div>
      <div className="flex flex-col gap-0 sm:flex-row">
        {STAGES.map((s, i) => (
          <div key={s.id} className="flex flex-1 items-stretch">
            <button
              onClick={() => setActive(i)}
              className={`flex-1 border px-3 py-3 text-left transition-colors ${
                active === i
                  ? "border-signal/60 bg-signal/10"
                  : "border-line bg-void hover:border-line-bright"
              } ${i > 0 ? "sm:-ml-px" : ""} ${i === 0 ? "rounded-l" : ""} ${
                i === STAGES.length - 1 ? "rounded-r" : ""
              }`}
            >
              <div className="font-mono text-[9px] uppercase tracking-widest text-ink-faint">
                {String(i + 1).padStart(2, "0")}
              </div>
              <div
                className={`mt-0.5 text-xs font-medium ${
                  active === i ? "text-signal" : "text-ink-dim"
                }`}
              >
                {s.label}
              </div>
            </button>
            {i < STAGES.length - 1 && (
              <div className="hidden items-center px-1 text-ink-faint sm:flex">→</div>
            )}
          </div>
        ))}
      </div>
      <div className="mt-4 rounded border border-line bg-void px-4 py-3.5">
        <div className="font-mono text-[11px] text-theory">{STAGES[active].mono}</div>
        <p className="mt-1.5 text-[13px] leading-relaxed text-ink-dim">{STAGES[active].detail}</p>
      </div>
    </div>
  );
}
