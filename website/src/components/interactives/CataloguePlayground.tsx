"use client";

import { useMemo, useState } from "react";
import candidates from "@/data/candidates.json";
import { Chart } from "@/components/charts/Chart";
import { RepoLink } from "@/components/ui/RepoLink";
import { C } from "@/lib/colors";
import { REPO_PATHS } from "@/lib/repo";

/**
 * Filter playground on the real candidate list produced by the repo's
 * catalog search (H>1, T90∈[0.2,5]s, no redshift → 36 sources).
 */
type Src = {
  name: string;
  ra: number;
  dec: number;
  hardness: number | null;
  lle_t90: number;
  gbm_t90: number;
};

export function CataloguePlayground() {
  const all = candidates.sources as Src[];
  const [t90Max, setT90Max] = useState(5);
  const [hMin, setHMin] = useState(1);
  const [useLLE, setUseLLE] = useState(true);
  const [selected, setSelected] = useState<string | null>(null);

  const passing = useMemo(
    () =>
      all.filter((s) => {
        const t90 = useLLE && s.lle_t90 > 0 ? s.lle_t90 : s.gbm_t90;
        const hOk = s.hardness == null ? hMin <= 1 : s.hardness >= hMin;
        return t90 <= t90Max && t90 > 0 && hOk;
      }),
    [all, t90Max, hMin, useLLE]
  );

  const sel = all.find((s) => s.name === selected);

  return (
    <div>
      <div className="mb-4 grid gap-x-8 gap-y-3 sm:grid-cols-2">
        <div className="flex items-center gap-3">
          <span className="w-32 font-mono text-[11px] text-ink-faint">max T90</span>
          <input
            type="range"
            min={0.2}
            max={5}
            step={0.1}
            value={t90Max}
            onChange={(e) => setT90Max(Number(e.target.value))}
            className="flex-1"
            aria-label="Maximum T90"
          />
          <span className="w-14 text-right font-mono text-xs text-signal">{t90Max.toFixed(1)}s</span>
        </div>
        <div className="flex items-center gap-3">
          <span className="w-32 font-mono text-[11px] text-ink-faint">min hardness</span>
          <input
            type="range"
            min={0}
            max={10}
            step={0.25}
            value={hMin}
            onChange={(e) => setHMin(Number(e.target.value))}
            className="flex-1"
            aria-label="Minimum hardness"
          />
          <span className="w-14 text-right font-mono text-xs text-signal">{hMin.toFixed(2)}</span>
        </div>
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={useLLE}
            onChange={(e) => setUseLLE(e.target.checked)}
            className="accent-[#7FD4C1]"
          />
          prefer LAT-LLE T90 when available
        </label>
        <div className="text-right font-mono text-xs">
          <span className="text-ink-faint">passing cuts · </span>
          <span className="text-2xl text-signal">{passing.length}</span>
          <span className="text-ink-faint"> / {all.length}</span>
        </div>
      </div>

      <Chart
        xDomain={[0.15, 60]}
        yDomain={[0.3, 400]}
        xLabel="T90 [s]"
        yLabel="hardness (LAT fluence / GBM fluence)"
        height={360}
      >
        {({ sx, sy }) => (
          <>
            {/* cut region shading */}
            <rect
              x={sx.range[0]}
              y={sy.range[1]}
              width={sx(t90Max) - sx.range[0]}
              height={sy(Math.max(hMin, 0.3)) - sy.range[1]}
              fill={C.signal}
              opacity="0.05"
            />
            <line
              x1={sx(t90Max)}
              x2={sx(t90Max)}
              y1={sy.range[0]}
              y2={sy.range[1]}
              stroke={C.signal}
              strokeDasharray="4 4"
              opacity="0.5"
            />
            {hMin > 0.3 && (
              <line
                x1={sx.range[0]}
                x2={sx.range[1]}
                y1={sy(hMin)}
                y2={sy(hMin)}
                stroke={C.signal}
                strokeDasharray="4 4"
                opacity="0.5"
              />
            )}
            {all.map((s) => {
              const t90 = useLLE && s.lle_t90 > 0 ? s.lle_t90 : s.gbm_t90;
              const h = s.hardness ?? 300; // inf plotted at top
              const pass = passing.includes(s);
              if (t90 <= 0) return null;
              return (
                <g
                  key={s.name}
                  onClick={() => setSelected(s.name === selected ? null : s.name)}
                  className="cursor-pointer"
                >
                  <circle
                    cx={sx(t90)}
                    cy={sy(Math.min(h, 350))}
                    r={selected === s.name ? 6 : 4}
                    fill={pass ? C.signal : C.line}
                    stroke={selected === s.name ? C.burst : pass ? C.signalDim : C.lineBright}
                    strokeWidth="1.5"
                    opacity={pass ? 0.95 : 0.55}
                  />
                  {s.hardness == null && (
                    <text
                      x={sx(t90)}
                      y={sy(350) - 8}
                      fontSize="8"
                      fill={C.inkFaint}
                      textAnchor="middle"
                      fontFamily="var(--font-plex-mono)"
                    >
                      ∞
                    </text>
                  )}
                </g>
              );
            })}
          </>
        )}
      </Chart>

      {sel && (
        <div className="mt-3 rounded border border-burst/40 bg-burst/5 px-4 py-3 font-mono text-xs">
          <span className="text-burst">{sel.name}</span>
          <span className="text-ink-faint">
            {" "}
            · RA {sel.ra.toFixed(2)}° · Dec {sel.dec.toFixed(2)}° · LLE T90 {sel.lle_t90}s · GBM T90{" "}
            {sel.gbm_t90}s · hardness {sel.hardness ?? "∞ (no GBM fluence)"}
          </span>
        </div>
      )}

      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Real sources from the{" "}
        <RepoLink path={REPO_PATHS.catalogSearch} file>
          repo&apos;s Fermi catalog search
        </RepoLink>{" "}
        (
        <RepoLink path={REPO_PATHS.candidatesCsv} file>
          H&gt;1_T90[0.2-5]_RS=0.csv
        </RepoLink>
        ): each passed hardness &gt; 1, T90 within 0.2–5 s in GBM or LAT-LLE, and no measured
        redshift. Sources marked ∞ had zero GBM-band fluence — the hardest events in the sample.
        Click a point to inspect it.
      </p>
    </div>
  );
}
