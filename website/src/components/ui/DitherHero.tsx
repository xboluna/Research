"use client";

import { useEffect, useRef } from "react";

/**
 * Bayer-dithered evaporating black hole.
 *
 * Renders a low-resolution luminance field (photon ring + Hawking glow +
 * starfield) and quantizes it through a 4×4 ordered-dither matrix, upscaled
 * with nearest-neighbor for a chunky retro-print feel.
 */
const BAYER4 = [
  [0, 8, 2, 10],
  [12, 4, 14, 6],
  [3, 11, 1, 9],
  [15, 7, 13, 5],
];

// palette: void → dim blue → mint signal → warm burst core
const PALETTE: [number, number, number][] = [
  [6, 7, 11],
  [18, 22, 34],
  [42, 48, 71],
  [61, 122, 110],
  [127, 212, 193],
  [232, 180, 90],
];

export function DitherHero({ className }: { className?: string }) {
  const ref = useRef<HTMLCanvasElement>(null);
  const raf = useRef(0);

  useEffect(() => {
    const canvas = ref.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d", { alpha: false });
    if (!ctx) return;

    const SCALE = 3; // dither cell size in CSS px
    let W = 0;
    let H = 0;
    let buf: ImageData | null = null;
    let stars: { x: number; y: number; b: number; tw: number }[] = [];

    function resize() {
      if (!canvas || !ctx) return;
      const rect = canvas.getBoundingClientRect();
      W = Math.max(80, Math.floor(rect.width / SCALE));
      H = Math.max(60, Math.floor(rect.height / SCALE));
      canvas.width = W;
      canvas.height = H;
      ctx.imageSmoothingEnabled = false;
      buf = ctx.createImageData(W, H);
      const n = Math.floor((W * H) / 320);
      stars = Array.from({ length: n }, () => ({
        x: Math.random() * W,
        y: Math.random() * H,
        b: 0.25 + Math.random() * 0.5,
        tw: Math.random() * Math.PI * 2,
      }));
    }

    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(canvas);

    const t0 = performance.now();

    function frame(now: number) {
      if (!ctx || !buf) return;
      const t = (now - t0) / 1000;
      const cx = W * 0.5;
      const cy = H * 0.46;
      const R = Math.min(W, H) * 0.16; // horizon radius
      const ring = R * 1.45;

      const data = buf.data;
      for (let y = 0; y < H; y++) {
        for (let x = 0; x < W; x++) {
          const dx = x - cx;
          const dy = (y - cy) * 1.12; // slight squash
          const r = Math.sqrt(dx * dx + dy * dy);
          const th = Math.atan2(dy, dx);

          let L = 0.055; // base sky

          // photon ring: bright thin annulus with rotating hot spot
          const ringDist = Math.abs(r - ring);
          const hot = 0.55 + 0.45 * Math.sin(th * 1 + t * 0.6);
          L += Math.exp(-(ringDist * ringDist) / (R * 0.16)) * (0.55 + 0.35 * hot);

          // Hawking glow halo, breathing slowly
          const breathe = 0.8 + 0.2 * Math.sin(t * 0.45);
          L += Math.exp(-r / (R * 3.1)) * 0.34 * breathe;

          // evaporation flicker: radial noise streaks
          const streak =
            Math.sin(th * 9 + t * 2.0 + Math.sin(r * 0.14)) *
            Math.sin(r * 0.22 - t * 2.6);
          if (r > ring * 0.95) L += Math.max(0, streak) * 0.10 * Math.exp(-(r - ring) / (R * 2.2));

          // event horizon: hard shadow
          if (r < R) L = Math.max(0, 0.5 - (R - r) / (R * 0.35)) * 0.35;

          // dither quantization
          const levels = PALETTE.length - 1;
          const noise = (BAYER4[y & 3][x & 3] + 0.5) / 16;
          const v = Math.max(0, Math.min(0.999, L));
          const idx = Math.min(levels, Math.floor(v * levels + noise));
          const [pr, pg, pb] = PALETTE[idx];

          const o = (y * W + x) * 4;
          data[o] = pr;
          data[o + 1] = pg;
          data[o + 2] = pb;
          data[o + 3] = 255;
        }
      }

      // starfield on top (twinkling single px)
      for (const s of stars) {
        const tw = 0.5 + 0.5 * Math.sin(t * 1.4 + s.tw);
        const b = s.b * tw;
        if (b < 0.18) continue;
        const o = ((s.y | 0) * W + (s.x | 0)) * 4;
        const boost = b > 0.5 ? PALETTE[4] : PALETTE[2];
        data[o] = boost[0];
        data[o + 1] = boost[1];
        data[o + 2] = boost[2];
      }

      ctx.putImageData(buf, 0, 0);
      raf.current = requestAnimationFrame(frame);
    }

    raf.current = requestAnimationFrame(frame);
    return () => {
      cancelAnimationFrame(raf.current);
      ro.disconnect();
    };
  }, []);

  return (
    <canvas
      ref={ref}
      className={className}
      style={{ imageRendering: "pixelated", width: "100%", height: "100%" }}
      aria-hidden
    />
  );
}
