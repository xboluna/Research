"use client";

import { useState } from "react";
import { C } from "@/lib/colors";

const ROWS: {
  claim: string;
  status: "supports" | "neutral" | "against";
  weight: string;
  detail: string;
}[] = [
  {
    claim: "Very short GRBs get harder as they get shorter",
    status: "supports",
    weight: "suggestive",
    detail:
      "Seen in BATSE (Cline 1996–97) and echoed in the repo's Fermi hardness–T90 analysis. Expected for PBH fireballs, but selection effects can mimic it.",
  },
  {
    claim: "Candidate sky distribution is isotropic, off-plane",
    status: "supports",
    weight: "suggestive",
    detail:
      "Both the GRB candidates and unassociated transients avoid the galactic plane — inconsistent with neutron-star populations, consistent with a local isotropic population.",
  },
  {
    claim: "36 Fermi sources pass all PBH burst cuts",
    status: "neutral",
    weight: "candidates only",
    detail:
      "Passing cuts is survival, not confirmation. Each needs the full template fit and, ideally, an afterglow or proper-motion cross-check.",
  },
  {
    claim: "1FLT transients fit the rising PBH lightcurve",
    status: "neutral",
    weight: "degenerate",
    detail:
      "35 unassociated LAT transients admit (τ, d) fits — but with huge error bars, and the two solution clusters imply proper motion the catalog doesn't see.",
  },
  {
    claim: "Expected proper motion is absent",
    status: "against",
    weight: "strong",
    detail:
      "At the fitted milliparsec distances the sources should streak degrees per year across the sky. The LAT transient catalog shows no such motion for these sources.",
  },
  {
    claim: "Theoretical explosion rate ≪ detector reach",
    status: "against",
    weight: "strong",
    detail:
      "For every allowed mass function, ṅ_PBH stays orders of magnitude below the HAWC bound of 3400 pc⁻³ yr⁻¹ — unless ψ(M) is finely spiked at M_U (Boluna et al. §3).",
  },
  {
    claim: "No radio afterglow counterparts found",
    status: "against",
    weight: "model-dependent",
    detail:
      "Low-frequency searches (ETA, Arecibo) constrain the Rees–Blandford radio pulse to < 2.3×10⁻⁷ pc⁻³ yr⁻¹ — though the modern SEM predicts a weaker pulse than Rees assumed.",
  },
];

const STATUS_META = {
  supports: { label: "hints", color: C.signal },
  neutral: { label: "open", color: C.burst },
  against: { label: "tension", color: C.danger },
};

export function EvidenceBoard() {
  const [open, setOpen] = useState<number | null>(null);

  return (
    <div className="space-y-1.5">
      {ROWS.map((r, i) => {
        const meta = STATUS_META[r.status];
        return (
          <button
            key={i}
            onClick={() => setOpen(open === i ? null : i)}
            className={`block w-full rounded border px-4 py-3 text-left transition-colors ${
              open === i ? "border-line-bright bg-raise" : "border-line bg-void hover:border-line-bright"
            }`}
          >
            <div className="flex items-center justify-between gap-3">
              <span className="text-[13px] text-ink">{r.claim}</span>
              <span className="flex shrink-0 items-center gap-2">
                <span className="font-mono text-[10px] text-ink-faint">{r.weight}</span>
                <span
                  className="rounded px-2 py-0.5 font-mono text-[10px] uppercase tracking-wider"
                  style={{ color: meta.color, border: `1px solid ${meta.color}44` }}
                >
                  {meta.label}
                </span>
              </span>
            </div>
            {open === i && (
              <p className="mt-2 text-[12.5px] leading-relaxed text-ink-dim">{r.detail}</p>
            )}
          </button>
        );
      })}
    </div>
  );
}
