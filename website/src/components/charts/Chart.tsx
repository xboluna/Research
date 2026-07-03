"use client";

import { ReactNode, useId } from "react";
import { C } from "@/lib/colors";
import {
  linTicks,
  logTicks,
  makeScale,
  Scale,
  ScaleKind,
  tickLabel,
} from "./scale";

export interface ChartFrame {
  sx: Scale;
  sy: Scale;
  w: number;
  h: number;
  inner: { x0: number; x1: number; y0: number; y1: number };
}

/**
 * Dark, dithered SVG chart frame with log/linear axes.
 * Children receive scales via render prop.
 */
export function Chart({
  width = 640,
  height = 380,
  xDomain,
  yDomain,
  xKind = "log",
  yKind = "log",
  xLabel,
  yLabel,
  margin = { top: 14, right: 16, bottom: 44, left: 62 },
  children,
  overlay,
}: {
  width?: number;
  height?: number;
  xDomain: [number, number];
  yDomain: [number, number];
  xKind?: ScaleKind;
  yKind?: ScaleKind;
  xLabel?: string;
  yLabel?: string;
  margin?: { top: number; right: number; bottom: number; left: number };
  children: (f: ChartFrame) => ReactNode;
  overlay?: ReactNode;
}) {
  const uid = useId().replace(/[:]/g, "");
  const inner = {
    x0: margin.left,
    x1: width - margin.right,
    y0: height - margin.bottom,
    y1: margin.top,
  };
  const sx = makeScale(xKind, xDomain, [inner.x0, inner.x1]);
  const sy = makeScale(yKind, yDomain, [inner.y0, inner.y1]);
  const xt = xKind === "log" ? logTicks(xDomain) : linTicks(xDomain);
  const yt = yKind === "log" ? logTicks(yDomain, 7) : linTicks(yDomain);

  return (
    <div className="relative w-full overflow-x-auto">
      <svg
        viewBox={`0 0 ${width} ${height}`}
        className="w-full min-w-[320px] select-none"
        style={{ maxWidth: width }}
      >
        <defs>
          {/* Bayer 4x4 dither pattern for area fills */}
          <pattern id={`d4-${uid}`} width="4" height="4" patternUnits="userSpaceOnUse">
            {[
              [0, 0], [2, 2], [2, 0], [0, 2],
            ].map(([x, y], i) => (
              <rect key={i} x={x} y={y} width="1" height="1" fill="currentColor" opacity={i < 2 ? 0.5 : 0.25} />
            ))}
          </pattern>
          <clipPath id={`clip-${uid}`}>
            <rect
              x={inner.x0}
              y={inner.y1}
              width={inner.x1 - inner.x0}
              height={inner.y0 - inner.y1}
            />
          </clipPath>
        </defs>

        {/* frame + grid */}
        <rect
          x={inner.x0}
          y={inner.y1}
          width={inner.x1 - inner.x0}
          height={inner.y0 - inner.y1}
          fill={C.void}
          stroke={C.line}
        />
        {xt.map((v) => (
          <line
            key={`gx${v}`}
            x1={sx(v)}
            x2={sx(v)}
            y1={inner.y0}
            y2={inner.y1}
            stroke={C.line}
            strokeDasharray="1 3"
          />
        ))}
        {yt.map((v) => (
          <line
            key={`gy${v}`}
            x1={inner.x0}
            x2={inner.x1}
            y1={sy(v)}
            y2={sy(v)}
            stroke={C.line}
            strokeDasharray="1 3"
          />
        ))}

        {/* ticks */}
        {xt.map((v) => (
          <text
            key={`tx${v}`}
            x={sx(v)}
            y={inner.y0 + 16}
            textAnchor="middle"
            fontSize="10"
            fill={C.inkFaint}
            fontFamily="var(--font-plex-mono)"
          >
            {tickLabel(v, xKind)}
          </text>
        ))}
        {yt.map((v) => (
          <text
            key={`ty${v}`}
            x={inner.x0 - 8}
            y={sy(v) + 3}
            textAnchor="end"
            fontSize="10"
            fill={C.inkFaint}
            fontFamily="var(--font-plex-mono)"
          >
            {tickLabel(v, yKind)}
          </text>
        ))}

        {/* axis labels */}
        {xLabel && (
          <text
            x={(inner.x0 + inner.x1) / 2}
            y={height - 8}
            textAnchor="middle"
            fontSize="11"
            fill={C.inkDim}
            fontFamily="var(--font-plex-mono)"
          >
            {xLabel}
          </text>
        )}
        {yLabel && (
          <text
            transform={`translate(14 ${(inner.y0 + inner.y1) / 2}) rotate(-90)`}
            textAnchor="middle"
            fontSize="11"
            fill={C.inkDim}
            fontFamily="var(--font-plex-mono)"
          >
            {yLabel}
          </text>
        )}

        <g clipPath={`url(#clip-${uid})`}>{children({ sx, sy, w: width, h: height, inner })}</g>
        {overlay}
      </svg>
    </div>
  );
}

export function Legend({
  items,
}: {
  items: { label: string; color: string; dash?: string }[];
}) {
  return (
    <div className="mt-3 flex flex-wrap gap-x-5 gap-y-1.5 font-mono text-[11px] text-ink-dim">
      {items.map((it) => (
        <span key={it.label} className="inline-flex items-center gap-1.5">
          <svg width="18" height="6">
            <line
              x1="0"
              y1="3"
              x2="18"
              y2="3"
              stroke={it.color}
              strokeWidth="2"
              strokeDasharray={it.dash}
            />
          </svg>
          {it.label}
        </span>
      ))}
    </div>
  );
}
