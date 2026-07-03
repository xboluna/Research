export const C = {
  void: "#06070B",
  surface: "#0C0E15",
  raise: "#12151F",
  line: "#1C2030",
  lineBright: "#2A3047",
  ink: "#E8EAF2",
  inkDim: "#8B91A7",
  inkFaint: "#565C72",
  signal: "#7FD4C1",
  signalDim: "#3D7A6E",
  burst: "#E8B45A",
  burstDim: "#8A6A35",
  theory: "#9D8CFF",
  theoryDim: "#5B519B",
  danger: "#E86A6A",
} as const;

/** Detector color assignments, consistent across all charts. */
export const DETECTOR_COLORS: Record<string, string> = {
  GBM: "#E8B45A",
  BATSE: "#B08D50",
  LAT: "#7FD4C1",
  HAWC: "#9D8CFF",
  VERITAS: "#E86A9C",
  LHAASO: "#6AAEE8",
  CTA: "#8BE86A",
};
