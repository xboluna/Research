"use client";

import { useState, type ReactNode } from "react";
import { RepoLink } from "@/components/ui/RepoLink";
import { REPO_PATHS } from "@/lib/repo";
import { C } from "@/lib/colors";

/**
 * Multi-messenger timeline: the γ flash is over in milliseconds; the radio
 * afterglow of the same event may arrive weeks to years later.
 */
const EVENTS: {
  t: number;
  band: string;
  color: string;
  label: string;
  detail: ReactNode;
}[] = [
  {
    t: -6,
    band: "γ",
    color: C.burst,
    label: "Final γ-ray flash",
    detail:
      "The last ~10⁹ g evaporate in about a microsecond–millisecond of TeV-scale emission. This is the 'burst' a GBM-like detector triggers on.",
  },
  {
    t: -1,
    band: "γ",
    color: C.burst,
    label: "Burst tail",
    detail:
      "Cline's fireball picture stretches the photon escape over the photosphere light-crossing time — up to a few hundred milliseconds.",
  },
  {
    t: 5.9,
    band: "radio",
    color: C.signal,
    label: "Relativistic shell sweeps ISM (~weeks)",
    detail:
      "Rees (1977): the e⁺e⁻ fireball acts as an expanding conductor, expelling the ambient magnetic field and radiating a coherent low-frequency pulse. Blandford (1977) put its spectrum near 1 GHz.",
  },
  {
    t: 7.5,
    band: "radio",
    color: C.signal,
    label: "Blazar-like afterglow peak (~months–years)",
    detail: (
      <>
        In adiabatic blob expansion, the GeV flare leads the radio by 40–140 days (observed for Mrk
        421). The repo asks: does a spherical PBH shell do the same? (
        <RepoLink path={REPO_PATHS.afterglowNotes} file>
          Afterglow.tex
        </RepoLink>
        )
      </>
    ),
  },
];

export function AfterglowTimeline() {
  const [sel, setSel] = useState(0);
  const e = EVENTS[sel];

  const toX = (logT: number) => ((logT + 7) / 16) * 100;

  return (
    <div>
      <div className="relative h-24">
        {/* dithered gradient track */}
        <div className="dither-50 absolute left-0 right-0 top-1/2 h-6 -translate-y-1/2 rounded-sm bg-void" />
        <div className="absolute left-0 right-0 top-1/2 h-px bg-line-bright" />
        {EVENTS.map((ev, i) => (
          <button
            key={i}
            onClick={() => setSel(i)}
            className="absolute top-1/2 -translate-x-1/2 -translate-y-1/2"
            style={{ left: `${toX(ev.t)}%` }}
            aria-label={ev.label}
          >
            <span
              className={`block rounded-full border-2 transition-all ${
                sel === i ? "h-5 w-5" : "h-3 w-3 hover:h-4 hover:w-4"
              }`}
              style={{
                borderColor: ev.color,
                background: sel === i ? ev.color : "transparent",
              }}
            />
          </button>
        ))}
        {/* band labels */}
        <span className="absolute left-[8%] top-1 font-mono text-[10px] text-burst">
          γ-ray (ms)
        </span>
        <span className="absolute right-[6%] top-1 font-mono text-[10px] text-signal">
          radio (weeks–years)
        </span>
      </div>

      <div className="flex justify-between font-mono text-[9px] text-ink-faint">
        <span>1 µs</span>
        <span>1 ms</span>
        <span>1 s</span>
        <span>1 hr</span>
        <span>1 mo</span>
        <span>30 yr</span>
      </div>

      <div className="mt-4 rounded border border-line bg-void px-4 py-3.5">
        <div className="flex items-baseline gap-3">
          <span
            className="rounded px-1.5 py-0.5 font-mono text-[10px] uppercase text-void"
            style={{ background: e.color }}
          >
            {e.band}
          </span>
          <span className="text-sm font-medium text-ink">{e.label}</span>
        </div>
        <p className="mt-2 text-[13px] leading-relaxed text-ink-dim">{e.detail}</p>
      </div>

      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        The observational strategy: catalog the position of every candidate burst, then watch
        those coordinates for a late-arriving radio transient. A match would be extraordinary —
        no astrophysical short GRB should produce this particular γ→radio sequence from a
        stationary point at parsec distance. Radio non-detections (ETA: &lt; 2.3×10⁻⁷ pc⁻³ yr⁻¹,
        Cutchin 2015) already constrain the Rees–Blandford channel.
      </p>
    </div>
  );
}
