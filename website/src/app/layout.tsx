import type { Metadata } from "next";
import { IBM_Plex_Mono, Space_Grotesk } from "next/font/google";
import "./globals.css";

const spaceGrotesk = Space_Grotesk({
  subsets: ["latin"],
  variable: "--font-space-grotesk",
  display: "swap",
});

const plexMono = IBM_Plex_Mono({
  subsets: ["latin"],
  weight: ["400", "500", "600"],
  variable: "--font-plex-mono",
  display: "swap",
});

export const metadata: Metadata = {
  title: "Searching for Exploding Black Holes",
  description:
    "An interactive journey through the hunt for evaporating primordial black holes — Hawking radiation, gamma-ray telescopes, and the Fermi mission catalogs. Based on Boluna et al. (arXiv:2307.06467) and the UCSC BHRad research program.",
  keywords: [
    "primordial black holes",
    "Hawking radiation",
    "gamma-ray bursts",
    "Fermi LAT",
    "dark matter",
  ],
  openGraph: {
    title: "Searching for Exploding Black Holes",
    description:
      "Interactive visualizations of primordial black hole evaporation and the gamma-ray search for their explosions.",
    type: "website",
  },
};

export default function RootLayout({
  children,
}: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en" className={`${spaceGrotesk.variable} ${plexMono.variable}`}>
      <body>{children}</body>
    </html>
  );
}
