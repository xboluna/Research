"use client";

import { useState } from "react";

const REFS: Record<string, { label: string; detail: string; href?: string }> = {
  thesis: {
    label: "Boluna 2023",
    detail: "Boluna, X. D. (2023). Detection Methods for Discovering Evaporating Primordial Black Holes in Modern Gamma-Ray Telescopes. M.S. thesis, UC Santa Cruz.",
    href: "https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf",
  },
  jcap: {
    label: "Boluna et al. 2024",
    detail: "Boluna, X., Profumo, S., Blé, J. & Hennings, D. Searching for Exploding Black Holes. arXiv:2307.06467.",
    href: "https://arxiv.org/abs/2307.06467",
  },
  cline97: {
    label: "Cline et al. 1997",
    detail: "Cline, D. B., Sanders, D. A. & Hong, W. (1997). Further Evidence for Some Gamma-Ray Bursts Consistent with Primordial Black Hole Evaporation. ApJ 486, 169.",
    href: "https://ui.adsabs.harvard.edu/abs/1997ApJ...486..169C",
  },
  hawking74: {
    label: "Hawking 1974",
    detail: "Hawking, S. W. (1974). Black hole explosions? Nature 248, 30.",
    href: "https://www.nature.com/articles/248030a0",
  },
  carr20: {
    label: "Carr & Kühnel 2020",
    detail: "Carr, B. & Kühnel, F. (2020). Primordial Black Holes as Dark Matter: Recent Developments. Annu. Rev. Nucl. Part. Sci. 70, 355. arXiv:2006.02838.",
    href: "https://arxiv.org/abs/2006.02838",
  },
  ukwatta16: {
    label: "Ukwatta et al. 2016",
    detail: "Ukwatta, T. N. et al. (2016). Primordial Black Holes: Observational characteristics of the final evaporation. Astropart. Phys. 80, 90. arXiv:1510.04372.",
    href: "https://arxiv.org/abs/1510.04372",
  },
  macgibbon90: {
    label: "MacGibbon & Webber 1990",
    detail: "MacGibbon, J. H. & Webber, B. R. (1990). Quark- and gluon-jet emission from primordial black holes. Phys. Rev. D 41, 3052.",
    href: "https://journals.aps.org/prd/abstract/10.1103/PhysRevD.41.3052",
  },
  blackhawk: {
    label: "BlackHawk",
    detail: "Arbey, A. & Auffinger, J. (2019). BlackHawk: a public code for calculating the Hawking evaporation spectra. Eur. Phys. J. C 79, 693. arXiv:1905.04268.",
    href: "https://arxiv.org/abs/1905.04268",
  },
  threeml: {
    label: "threeML",
    detail: "Vianello, G. et al. (2015). The Multi-Mission Maximum Likelihood framework (3ML). arXiv:1507.08343.",
    href: "https://threeml.readthedocs.io",
  },
  rees77: {
    label: "Rees 1977",
    detail: "Rees, M. J. (1977). A better way of searching for black-hole explosions? Nature 266, 333.",
    href: "https://www.nature.com/articles/266333a0",
  },
  blandford77: {
    label: "Blandford 1977",
    detail: "Blandford, R. D. (1977). Spectrum of a radio pulse from an exploding black hole. MNRAS 181, 489.",
  },
  cutchin15: {
    label: "Cutchin et al. 2015",
    detail: "Cutchin, S. E. et al. (2015). Constraining the Rate of Primordial Black Hole Explosions Using Low-Frequency Radio Observations. PASP 127, 1269.",
  },
  hawc: {
    label: "HAWC 2020",
    detail: "Albert, A. et al. (HAWC) (2020). Constraining the local burst rate density of primordial black holes with HAWC. JCAP 04, 026.",
    href: "https://arxiv.org/abs/1911.04356",
  },
  page76: {
    label: "Page & Hawking 1976",
    detail: "Page, D. N. & Hawking, S. W. (1976). Gamma rays from primordial black holes. ApJ 206, 1.",
  },
  cline02: {
    label: "Cline et al. 2002",
    detail: "Cline, D. B., Matthey, C. & Otwinowski, S. (2003). Evidence for a Galactic Origin of Very Short Gamma-Ray Bursts. Astropart. Phys. 18, 531.",
  },
  fermilat18: {
    label: "Fermi-LAT 2018",
    detail: "Ackermann, M. et al. (2018). Search for Gamma-Ray Emission from Local Primordial Black Holes with the Fermi Large Area Telescope. ApJ 857, 49.",
    href: "https://arxiv.org/abs/1802.00100",
  },
};

export function Cite({ k, sec }: { k: keyof typeof REFS | string; sec?: string }) {
  const [open, setOpen] = useState(false);
  const ref = REFS[k];
  if (!ref) return null;
  return (
    <span className="relative inline-block">
      <button
        onClick={() => setOpen((v) => !v)}
        onBlur={() => setTimeout(() => setOpen(false), 150)}
        className="text-signal/80 hover:text-signal text-[0.85em] font-mono cursor-pointer border-b border-dotted border-signal/40"
      >
        [{ref.label}
        {sec ? `, ${sec}` : ""}]
      </button>
      {open && (
        <span className="absolute z-50 bottom-full left-1/2 -translate-x-1/2 mb-2 w-72 rounded border border-line-bright bg-raise p-3 text-xs leading-relaxed text-ink-dim shadow-xl">
          {ref.detail}
          {ref.href && (
            <a
              href={ref.href}
              target="_blank"
              rel="noreferrer"
              className="block mt-1.5 text-signal hover:underline font-mono"
            >
              {ref.href.replace("https://", "").slice(0, 40)}…
            </a>
          )}
        </span>
      )}
    </span>
  );
}
