"use client";

import { useState } from "react";
import { horizonMassG } from "@/lib/physics";
import { duration, mass, sci } from "@/lib/format";
import { C } from "@/lib/colors";

const EPOCHS = [
  {
    logT: -36,
    name: "Inflation ends",
    detail:
      "Quantum fluctuations stretched to cosmic scales seed all future structure — including, possibly, PBHs. Holes formed here would be lighter than 1 g and long gone.",
  },
  {
    logT: -23,
    name: "PBHs of ~10¹⁵ g",
    detail:
      "An overdense patch collapsing when the universe is 10⁻²³ s old traps roughly a horizon mass: ~10¹⁵ g. These are the holes whose explosions we could witness today.",
  },
  {
    logT: -12,
    name: "Electroweak transition",
    detail:
      "The Higgs field switches on. Horizon-mass holes formed now would weigh ~10²⁶ g — sublunar mass, and still viable dark-matter candidates.",
  },
  {
    logT: -5,
    name: "QCD transition",
    detail:
      "Quarks bind into protons. The equation of state softens, making collapse easier — a natural epoch for enhanced PBH formation near one solar mass.",
  },
  {
    logT: 1,
    name: "Nucleosynthesis begins",
    detail:
      "Any holes light enough to evaporate during BBN would disrupt the light-element abundances — this is why 10⁹–10¹³ g PBHs are tightly constrained.",
  },
];

export function FormationTimeline() {
  const [logT, setLogT] = useState(-23);
  const t = Math.pow(10, logT);
  const MH = horizonMassG(t);
  const active = EPOCHS.reduce((best, e) =>
    Math.abs(e.logT - logT) < Math.abs(best.logT - logT) ? e : best
  );
  const nearActive = Math.abs(active.logT - logT) < 2.5;

  return (
    <div>
      {/* timeline track */}
      <div className="relative mb-2 h-14">
        <div className="absolute left-0 right-0 top-1/2 h-px bg-line-bright" />
        {EPOCHS.map((e) => {
          const x = ((e.logT + 38) / 41) * 100;
          return (
            <button
              key={e.name}
              onClick={() => setLogT(e.logT)}
              className="group absolute top-1/2 -translate-x-1/2 -translate-y-1/2"
              style={{ left: `${x}%` }}
              aria-label={e.name}
            >
              <span
                className={`block h-2.5 w-2.5 rotate-45 border transition-colors ${
                  active.name === e.name && nearActive
                    ? "border-signal bg-signal"
                    : "border-line-bright bg-surface group-hover:border-signal"
                }`}
              />
            </button>
          );
        })}
        {/* cursor */}
        <div
          className="absolute top-0 h-full w-px bg-burst/70"
          style={{ left: `${((logT + 38) / 41) * 100}%` }}
        />
      </div>

      <input
        type="range"
        min={-38}
        max={3}
        step={0.1}
        value={logT}
        onChange={(e) => setLogT(Number(e.target.value))}
        className="w-full"
        aria-label="Cosmic time"
      />

      <div className="mt-4 grid gap-3 sm:grid-cols-2">
        <div className="rounded border border-line bg-void px-4 py-3">
          <div className="font-mono text-[10px] uppercase tracking-widest text-ink-faint">
            cosmic time
          </div>
          <div className="mt-1 font-mono text-lg text-ink">{duration(t)}</div>
          <div className="mt-2 font-mono text-[10px] uppercase tracking-widest text-ink-faint">
            horizon mass → PBH mass
          </div>
          <div className="mt-1 font-mono text-lg" style={{ color: C.signal }}>
            {mass(MH)}
          </div>
          <div className="mt-1 text-[11px] text-ink-faint">
            M_H ≈ 10¹⁵ g × (t / 10⁻²³ s) — collapse when a fluctuation enters the horizon
          </div>
        </div>
        <div
          className={`rounded border px-4 py-3 transition-colors ${
            nearActive ? "border-signal/40 bg-signal/5" : "border-line bg-void"
          }`}
        >
          <div className="font-mono text-[10px] uppercase tracking-widest text-ink-faint">
            {nearActive ? "epoch" : "nearest epoch"}
          </div>
          <div className="mt-1 text-sm font-medium text-ink">{active.name}</div>
          <p className="mt-1.5 text-[12px] leading-relaxed text-ink-dim">{active.detail}</p>
        </div>
      </div>

      {Math.abs(MH / 5.1e14 - 1) < 0.7 && (
        <p className="mt-3 rounded border border-burst/40 bg-burst/5 px-3 py-2 font-mono text-[11px] text-burst">
          ★ Holes born right now (≈{sci(5.1e14)} g) have a lifetime equal to the age of the
          universe — they are exploding <em>today</em>.
        </p>
      )}
    </div>
  );
}
