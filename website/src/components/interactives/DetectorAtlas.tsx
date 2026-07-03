"use client";

import { useState } from "react";
import areas from "@/data/effective-areas.json";
import bands from "@/data/telescope-bands.json";
import spectra from "@/data/spectra.json";
import { Chart, Legend } from "@/components/charts/Chart";
import { linePath } from "@/components/charts/scale";
import { DETECTOR_COLORS, C } from "@/lib/colors";
import { mass } from "@/lib/format";

type AreaMap = Record<string, { E_GeV: number[]; Aeff_cm2: number[] }>;

export function DetectorAtlas() {
  const [enabled, setEnabled] = useState<Record<string, boolean>>({
    GBM: true,
    BATSE: false,
    LAT: true,
    HAWC: true,
    VERITAS: true,
    LHAASO: true,
  });
  const [spectrumMass, setSpectrumMass] = useState<number | null>(2); // index; null = off

  const A = areas as AreaMap;
  const names = Object.keys(A).sort();
  const entry = spectrumMass != null ? spectra.masses[spectrumMass] : null;
  // scale spectrum into Aeff axis for visual comparison (arbitrary norm)
  const specScaled = entry
    ? spectra.energy_GeV.map((E: number, i: number) => {
        const v = (entry.direct[i] + entry.frag[i]) * E; // E·dN/dE shape
        return v > 0 ? v * 1e-16 : null;
      })
    : null;

  return (
    <div>
      <div className="mb-4 flex flex-wrap items-center gap-2">
        {names.map((n) => (
          <button
            key={n}
            onClick={() => setEnabled((s) => ({ ...s, [n]: !s[n] }))}
            className={`rounded border px-2.5 py-1 font-mono text-[11px] transition-colors ${
              enabled[n]
                ? "border-transparent text-void"
                : "border-line text-ink-faint hover:border-line-bright"
            }`}
            style={enabled[n] ? { background: DETECTOR_COLORS[n] } : undefined}
          >
            {n}
          </button>
        ))}
        <span className="mx-2 h-4 w-px bg-line-bright" />
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={spectrumMass != null}
            onChange={(e) => setSpectrumMass(e.target.checked ? 2 : null)}
            className="accent-[#7FD4C1]"
          />
          overlay PBH spectrum shape
        </label>
        {spectrumMass != null && (
          <>
            <input
              type="range"
              min={0}
              max={spectra.masses.length - 1}
              step={1}
              value={spectrumMass}
              onChange={(e) => setSpectrumMass(Number(e.target.value))}
              className="w-28"
              aria-label="Spectrum mass"
            />
            <span className="font-mono text-[11px] text-signal">
              {mass(spectra.masses[spectrumMass].mass_g)}
            </span>
          </>
        )}
      </div>

      <Chart
        xDomain={[1e-6, 1e6]}
        yDomain={[1, 1e11]}
        xLabel="photon energy [GeV]"
        yLabel="effective area [cm²]"
      >
        {({ sx, sy }) => (
          <>
            {names.map(
              (n) =>
                enabled[n] && (
                  <path
                    key={n}
                    d={linePath(A[n].E_GeV, A[n].Aeff_cm2, sx, sy)}
                    fill="none"
                    stroke={DETECTOR_COLORS[n]}
                    strokeWidth="1.8"
                  />
                )
            )}
            {specScaled && (
              <path
                d={linePath(spectra.energy_GeV, specScaled, sx, sy)}
                fill="none"
                stroke={C.ink}
                strokeWidth="1.2"
                strokeDasharray="2 3"
                opacity="0.65"
              />
            )}
          </>
        )}
      </Chart>

      <Legend
        items={[
          ...names.filter((n) => enabled[n]).map((n) => ({ label: n, color: DETECTOR_COLORS[n] })),
          ...(specScaled
            ? [{ label: "E·dN/dE spectrum shape (arb. norm)", color: C.ink, dash: "2 3" }]
            : []),
        ]}
      />

      {/* energy band ruler */}
      <div className="mt-6">
        <div className="mb-1 font-mono text-[10px] uppercase tracking-widest text-ink-faint">
          energy coverage
        </div>
        <div className="space-y-1">
          {(bands.bands as { name: string; lo_GeV: number; hi_GeV: number; kind: string }[]).map(
            (b) => {
              const lo = (Math.log10(b.lo_GeV) + 6) / 12;
              const hi = (Math.log10(b.hi_GeV) + 6) / 12;
              const key = b.name.replace("Fermi ", "");
              return (
                <div key={b.name} className="flex items-center gap-2">
                  <span className="w-20 text-right font-mono text-[10px] text-ink-faint">
                    {b.name}
                  </span>
                  <div className="relative h-3 flex-1 rounded-sm bg-void dither-25">
                    <div
                      className="absolute h-full rounded-sm opacity-80"
                      style={{
                        left: `${lo * 100}%`,
                        width: `${(hi - lo) * 100}%`,
                        background: DETECTOR_COLORS[key] ?? C.inkDim,
                      }}
                    />
                  </div>
                </div>
              );
            }
          )}
          <div className="ml-[5.5rem] flex justify-between font-mono text-[9px] text-ink-faint">
            <span>keV</span>
            <span>MeV</span>
            <span>GeV</span>
            <span>TeV</span>
            <span>PeV</span>
          </div>
        </div>
      </div>
    </div>
  );
}
