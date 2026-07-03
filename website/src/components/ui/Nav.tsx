"use client";

import { useEffect, useState } from "react";

const SECTIONS = [
  { id: "hero", label: "Start" },
  { id: "hawking", label: "01 Hawking" },
  { id: "origins", label: "02 Origins" },
  { id: "collider", label: "03 Collider" },
  { id: "telescopes", label: "04 Telescopes" },
  { id: "distance", label: "05 Distance" },
  { id: "fitting", label: "06 Fitting" },
  { id: "catalogue", label: "07 Catalogue" },
  { id: "sky", label: "08 Sky" },
  { id: "afterglow", label: "09 Afterglow" },
  { id: "transients", label: "10 Transients" },
  { id: "verdict", label: "11 Verdict" },
];

export function Nav() {
  const [active, setActive] = useState("hero");
  const [open, setOpen] = useState(false);

  useEffect(() => {
    const obs = new IntersectionObserver(
      (entries) => {
        for (const e of entries) {
          if (e.isIntersecting) setActive(e.target.id);
        }
      },
      { rootMargin: "-35% 0px -55% 0px" }
    );
    for (const s of SECTIONS) {
      const el = document.getElementById(s.id);
      if (el) obs.observe(el);
    }
    return () => obs.disconnect();
  }, []);

  return (
    <>
      {/* desktop rail */}
      <nav className="fixed left-4 top-1/2 z-40 hidden -translate-y-1/2 lg:block">
        <ul className="space-y-1">
          {SECTIONS.map((s) => (
            <li key={s.id}>
              <a
                href={`#${s.id}`}
                className={`group flex items-center gap-2 py-0.5 font-mono text-[10px] tracking-wider transition-colors ${
                  active === s.id ? "text-signal" : "text-ink-faint hover:text-ink-dim"
                }`}
              >
                <span
                  className={`h-px transition-all ${
                    active === s.id ? "w-6 bg-signal" : "w-3 bg-line-bright group-hover:w-5"
                  }`}
                />
                {s.label}
              </a>
            </li>
          ))}
        </ul>
      </nav>

      {/* mobile toggle */}
      <button
        onClick={() => setOpen((v) => !v)}
        className="fixed right-4 top-4 z-50 rounded border border-line-bright bg-surface/90 px-3 py-2 font-mono text-xs text-ink-dim backdrop-blur lg:hidden"
      >
        {open ? "close" : "menu"}
      </button>
      {open && (
        <nav className="fixed inset-0 z-40 bg-void/95 backdrop-blur lg:hidden">
          <ul className="flex h-full flex-col items-center justify-center gap-3">
            {SECTIONS.map((s) => (
              <li key={s.id}>
                <a
                  href={`#${s.id}`}
                  onClick={() => setOpen(false)}
                  className={`font-mono text-sm ${
                    active === s.id ? "text-signal" : "text-ink-dim"
                  }`}
                >
                  {s.label}
                </a>
              </li>
            ))}
          </ul>
        </nav>
      )}
    </>
  );
}
