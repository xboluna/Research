"use client";

import { useMemo, useState } from "react";
import candidates from "@/data/candidates.json";
import transients from "@/data/transients.json";
import { RepoLink } from "@/components/ui/RepoLink";
import { REPO_PATHS } from "@/lib/repo";
import { C } from "@/lib/colors";

/**
 * Mollweide sky map in galactic coordinates with GRB candidates and
 * LAT transients. Mirrors generate_mollweise() from BoresightSelection.py.
 */

function mollweide(lonDeg: number, latDeg: number): [number, number] {
  // lon in [-180,180], lat in [-90,90] → x in [-2√2,2√2], y in [-√2,√2]
  const lon = (lonDeg * Math.PI) / 180;
  const lat = (Math.max(-89.99, Math.min(89.99, latDeg)) * Math.PI) / 180;
  let theta = lat;
  for (let i = 0; i < 12; i++) {
    const denom = 2 + 2 * Math.cos(2 * theta);
    if (Math.abs(denom) < 1e-9) break; // at poles theta = lat already
    const dt = -(2 * theta + Math.sin(2 * theta) - Math.PI * Math.sin(lat)) / denom;
    theta += dt;
    if (Math.abs(dt) < 1e-8) break;
  }
  const x = (2 * Math.SQRT2 * lon * Math.cos(theta)) / Math.PI;
  const y = Math.SQRT2 * Math.sin(theta);
  return [x, y];
}

/** Convert equatorial RA/Dec (deg) to galactic l/b (deg). */
function eqToGal(raDeg: number, decDeg: number): [number, number] {
  const ra = (raDeg * Math.PI) / 180;
  const dec = (decDeg * Math.PI) / 180;
  // J2000 north galactic pole
  const raGP = (192.85948 * Math.PI) / 180;
  const decGP = (27.12825 * Math.PI) / 180;
  const lNCP = (122.93192 * Math.PI) / 180;
  const sinb =
    Math.sin(decGP) * Math.sin(dec) + Math.cos(decGP) * Math.cos(dec) * Math.cos(ra - raGP);
  const b = Math.asin(sinb);
  const y = Math.cos(dec) * Math.sin(ra - raGP);
  const x = Math.cos(decGP) * Math.sin(dec) - Math.sin(decGP) * Math.cos(dec) * Math.cos(ra - raGP);
  let l = lNCP - Math.atan2(y, x);
  l = ((l % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
  return [(l * 180) / Math.PI, (b * 180) / Math.PI];
}

const W = 680;
const H = 360;
const SX = W / (4 * Math.SQRT2 + 0.4);
const SY = H / (2 * Math.SQRT2 + 0.3);

function toPx(l: number, b: number): [number, number] {
  const lon = l > 180 ? l - 360 : l; // center on galactic center
  const [mx, my] = mollweide(-lon, b); // flip so +l is left (astronomy convention)
  return [W / 2 + mx * SX, H / 2 - my * SY];
}

export function SkyMap() {
  const [showGRBs, setShowGRBs] = useState(true);
  const [showTransients, setShowTransients] = useState(true);

  const grbPts = useMemo(
    () =>
      (candidates.sources as { name: string; ra: number; dec: number }[]).map((s) => {
        const [l, b] = eqToGal(s.ra, s.dec);
        return { name: s.name, l, b, px: toPx(l, b) };
      }),
    []
  );
  const trPts = useMemo(
    () =>
      (transients.sources as { name: string; lii: number; bii: number }[]).map((s) => ({
        name: s.name,
        l: s.lii,
        b: s.bii,
        px: toPx(s.lii, s.bii),
      })),
    []
  );

  // hemisphere balance (BoresightSelection.py caption logic)
  const stats = useMemo(() => {
    const pts = [
      ...(showGRBs ? grbPts : []),
      ...(showTransients ? trPts : []),
    ];
    const n = pts.filter((p) => p.b > 0).length;
    const s = pts.filter((p) => p.b < 0).length;
    const lonN = (l: number) => (l > 180 ? l - 360 : l);
    const e = pts.filter((p) => lonN(p.l) > 0).length;
    const w = pts.filter((p) => lonN(p.l) < 0).length;
    return { n, s, e, w };
  }, [grbPts, trPts, showGRBs, showTransients]);

  // graticule
  const grat = useMemo(() => {
    const lines: string[] = [];
    for (const b of [-60, -30, 0, 30, 60]) {
      let d = "";
      for (let l = -180; l <= 180; l += 5) {
        const [x, y] = toPx(((l + 360) % 360), b);
        d += (l === -180 ? "M" : "L") + x.toFixed(1) + "," + y.toFixed(1);
      }
      lines.push(d);
    }
    for (const l of [0, 60, 120, 180, 240, 300]) {
      let d = "";
      for (let b = -90; b <= 90; b += 5) {
        const [x, y] = toPx(l, b);
        d += (b === -90 ? "M" : "L") + x.toFixed(1) + "," + y.toFixed(1);
      }
      lines.push(d);
    }
    return lines;
  }, []);

  return (
    <div>
      <div className="mb-4 flex flex-wrap gap-x-6 gap-y-2">
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={showGRBs}
            onChange={(e) => setShowGRBs(e.target.checked)}
            className="accent-[#E8B45A]"
          />
          GRB candidates ({grbPts.length})
        </label>
        <label className="flex cursor-pointer items-center gap-2 font-mono text-[11px] text-ink-dim">
          <input
            type="checkbox"
            checked={showTransients}
            onChange={(e) => setShowTransients(e.target.checked)}
            className="accent-[#E86A6A]"
          />
          1FLT unassociated transients ({trPts.length})
        </label>
      </div>

      <div className="overflow-x-auto">
        <svg viewBox={`0 0 ${W} ${H}`} className="w-full min-w-[320px]" style={{ maxWidth: W }}>
          {/* ellipse boundary */}
          <ellipse
            cx={W / 2}
            cy={H / 2}
            rx={2 * Math.SQRT2 * SX}
            ry={Math.SQRT2 * SY}
            fill={C.void}
            stroke={C.lineBright}
          />
          {grat.map((d, i) => (
            <path key={i} d={d} fill="none" stroke={C.line} strokeWidth="0.6" />
          ))}
          {/* galactic plane emphasis */}
          <path d={grat[2]} fill="none" stroke={C.theoryDim} strokeWidth="1.2" opacity="0.7" />
          {/* galactic center */}
          <circle cx={W / 2} cy={H / 2} r="3" fill="none" stroke={C.theory} strokeWidth="1" />
          <text
            x={W / 2 + 7}
            y={H / 2 + 3.5}
            fontSize="8.5"
            fill={C.theoryDim}
            fontFamily="var(--font-plex-mono)"
          >
            GC
          </text>

          {showTransients &&
            trPts.map((p) => (
              <g key={p.name} transform={`translate(${p.px[0]} ${p.px[1]})`}>
                <line x1="-3.5" y1="-3.5" x2="3.5" y2="3.5" stroke={C.danger} strokeWidth="1.3" />
                <line x1="-3.5" y1="3.5" x2="3.5" y2="-3.5" stroke={C.danger} strokeWidth="1.3" />
              </g>
            ))}
          {showGRBs &&
            grbPts.map((p) => (
              <circle
                key={p.name}
                cx={p.px[0]}
                cy={p.px[1]}
                r="3.4"
                fill={C.burst}
                opacity="0.85"
              >
                <title>{p.name}</title>
              </circle>
            ))}
        </svg>
      </div>

      <div className="mt-3 grid gap-3 font-mono text-xs sm:grid-cols-2">
        <div className="rounded border border-line bg-void px-3 py-2.5">
          <span className="text-ink-faint">north / south of plane · </span>
          <span className="text-ink">
            {stats.n} / {stats.s}
          </span>
          <span className="text-ink-faint">
            {" "}
            (ratio {(stats.n / Math.max(stats.s, 1)).toFixed(2)})
          </span>
        </div>
        <div className="rounded border border-line bg-void px-3 py-2.5">
          <span className="text-ink-faint">toward / away from center · </span>
          <span className="text-ink">
            {stats.e} / {stats.w}
          </span>
          <span className="text-ink-faint">
            {" "}
            (ratio {(stats.e / Math.max(stats.w, 1)).toFixed(2)})
          </span>
        </div>
      </div>

      <p className="mt-3 max-w-2xl font-mono text-[11px] leading-relaxed text-ink-faint">
        Galactic coordinates, galactic center at the middle, plane along the equator. If these
        sources were neutron stars or other stellar remnants they would trace the disk; local PBH
        explosions should be isotropic. Sky positions from the{" "}
        <RepoLink path={REPO_PATHS.candidatesCsv} file>
          candidate CSV
        </RepoLink>{" "}
        and the 1FLT fitted-parameter table (cf.{" "}
        <RepoLink path={REPO_PATHS.boresightSelection} file>
          BoresightSelection.py
        </RepoLink>
        , thesis Fig. 3.5).
      </p>
    </div>
  );
}
