"use client";

import { useState } from "react";
import { lifetimeS, temperatureGeV } from "@/lib/physics";
import { duration, energy, mass, sci } from "@/lib/format";
import { C } from "@/lib/colors";

const COMPARISONS: { mass: number; label: string }[] = [
  { mass: 1e9, label: "≈ a small asteroid" },
  { mass: 1e11, label: "≈ Halley's comet nucleus" },
  { mass: 5.1e14, label: "≈ a large mountain — evaporating today" },
  { mass: 1e17, label: "≈ 1/100 the Moon's mass" },
  { mass: 2e33, label: "≈ the Sun" },
];

function bestDetector(T_GeV: number): { name: string; color: string } {
  if (T_GeV > 100) return { name: "ground arrays (HAWC · LHAASO · VERITAS)", color: C.theory };
  if (T_GeV > 0.05) return { name: "Fermi LAT", color: C.signal };
  if (T_GeV > 1e-5) return { name: "Fermi GBM", color: C.burst };
  return { name: "X-ray / longer-wavelength missions", color: C.inkDim };
}

export function MassDial() {
  const [logM, setLogM] = useState(14.7); // ≈ M_U

  const M = Math.pow(10, logM);
  const T = temperatureGeV(M);
  const tau = lifetimeS(M);
  const det = bestDetector(T);
  const cmp = COMPARISONS.reduce((best, c) =>
    Math.abs(Math.log10(c.mass) - logM) < Math.abs(Math.log10(best.mass) - logM) ? c : best
  );
  const ageU = 4.35e17; // s

  return (
    <div>
      <div className="mb-5 flex items-center gap-4">
        <span className="font-mono text-[11px] uppercase tracking-widest text-ink-faint">
          Black hole mass
        </span>
        <input
          type="range"
          min={8}
          max={17}
          step={0.05}
          value={logM}
          onChange={(e) => setLogM(Number(e.target.value))}
          className="flex-1"
          aria-label="Black hole mass (log grams)"
        />
        <span className="w-24 text-right font-mono text-sm text-signal">{mass(M)}</span>
      </div>

      <div className="grid gap-3 sm:grid-cols-3">
        <Stat label="Hawking temperature" value={energy(T)} sub={`k_B T = ħc³ / 8πGM`} />
        <Stat
          label="Remaining lifetime"
          value={duration(tau)}
          sub={
            tau > ageU
              ? `${sci(tau / ageU)}× the age of the universe`
              : tau > 3.156e10
                ? "longer than human history"
                : tau > 3.156e7
                  ? "watchable within a career"
                  : "explodes on human timescales"
          }
          highlight={tau < 3.156e7 * 100}
        />
        <Stat
          label="Best current detector"
          value={det.name}
          color={det.color}
          sub={`peak emission near ${energy(T * 4.3)}`}
        />
      </div>

      <p className="mt-4 font-mono text-[11px] text-ink-faint">
        scale: {cmp.label}
        {Math.abs(logM - 14.71) < 0.15 && (
          <span className="ml-2 text-signal">← M_U = 5.1×10¹⁴ g: lifetime = age of the universe</span>
        )}
      </p>
    </div>
  );
}

function Stat({
  label,
  value,
  sub,
  color,
  highlight,
}: {
  label: string;
  value: string;
  sub?: string;
  color?: string;
  highlight?: boolean;
}) {
  return (
    <div
      className={`rounded border px-3.5 py-3 ${
        highlight ? "border-burst/50 bg-burst/5" : "border-line bg-void"
      }`}
    >
      <div className="font-mono text-[10px] uppercase tracking-widest text-ink-faint">{label}</div>
      <div className="mt-1 font-mono text-base" style={{ color: color ?? C.ink }}>
        {value}
      </div>
      {sub && <div className="mt-0.5 text-[11px] text-ink-faint">{sub}</div>}
    </div>
  );
}
