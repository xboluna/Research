import { Hero } from "@/components/sections/Hero";
import { Nav } from "@/components/ui/Nav";
import { Note, Panel, Prose, Section } from "@/components/ui/Section";
import { TeX } from "@/components/ui/Math";
import { Cite } from "@/components/ui/Cite";
import { MassDial } from "@/components/interactives/MassDial";
import { SpectrumExplorer } from "@/components/interactives/SpectrumExplorer";
import { BackwardsBurst } from "@/components/interactives/BackwardsBurst";
import { FormationTimeline } from "@/components/interactives/FormationTimeline";
import { MassFunctions } from "@/components/interactives/MassFunctions";
import { AlphaStaircase } from "@/components/interactives/AlphaStaircase";
import { SpectralIndex } from "@/components/interactives/SpectralIndex";
import { DetectorAtlas } from "@/components/interactives/DetectorAtlas";
import { ThreeMLPipeline } from "@/components/interactives/ThreeMLPipeline";
import { DistanceFrontier } from "@/components/interactives/DistanceFrontier";
import { Neighborhood } from "@/components/interactives/Neighborhood";
import { LightcurveComposer } from "@/components/interactives/LightcurveComposer";
import { CataloguePlayground } from "@/components/interactives/CataloguePlayground";
import { SkyMap } from "@/components/interactives/SkyMap";
import { AfterglowTimeline } from "@/components/interactives/AfterglowTimeline";
import { TransientFitter } from "@/components/interactives/TransientFitter";
import { EvidenceBoard } from "@/components/interactives/EvidenceBoard";
import { RepoLink } from "@/components/ui/RepoLink";
import { REPO_PATHS } from "@/lib/repo";

export default function Home() {
  return (
    <main>
      <Nav />
      <Hero />

      {/* ------------------------------------------------ 01 HAWKING */}
      <Section
        id="hawking"
        index="01"
        kicker="the physics of dying black holes"
        title="Black holes are not quite black"
      >
        <Prose>
          <p>
            In 1974, Hawking showed that quantum effects near the event horizon force black holes
            to radiate like thermal bodies <Cite k="hawking74" />. The temperature is set by
            nothing but the mass:
          </p>
          <TeX block math="k_B T_{\mathrm{BH}} = \frac{\hbar c^3}{8\pi G M} \;\approx\; \left(\frac{10^{16}\,\mathrm{g}}{M}\right)\,\mathrm{MeV}" />
          <p>
            The rule is merciless: <strong>smaller means hotter</strong>. A stellar-mass black
            hole radiates at a hundred-millionth of a kelvin — utterly invisible. But shrink one
            to the mass of a mountain and it glows in γ-rays. And because radiating mass makes it
            smaller still, evaporation <em>runs away</em>. The final moments are an explosion.
          </p>
        </Prose>

        <Panel title="dial a black hole" caption="Lifetime from τ ≈ 407 s × (M/10¹⁰ g)³ — Boluna et al. Eq. (4.4). The 'best detector' targets the band where emission peaks (≈ 4.3 k_BT for a photon blackbody-like spectrum).">
          <MassDial />
        </Panel>

        <Prose>
          <p>
            What does the light of a dying black hole look like? Two components{" "}
            <Cite k="ukwatta16" />
            <Cite k="macgibbon90" />: <em>direct</em> photons radiated straight off the horizon,
            and a <em>fragmentation</em> component — quarks and gluons that hadronize into pions,
            whose decays flood the spectrum with γ-rays. For any hole hotter than the QCD scale,
            fragmentation dominates.
          </p>
        </Prose>

        <Panel
          title="the photon spectrum, mass by mass"
          caption={
            <>
              Grid computed by the site&apos;s{" "}
              <RepoLink path={REPO_PATHS.dataPipeline} file>
                data pipeline
              </RepoLink>{" "}
              with the exact parameterization used in{" "}
              <RepoLink path={REPO_PATHS.analyticalModelling} file>
                Analytical_Modelling.ipynb
              </RepoLink>{" "}
              (Ukwatta et al. 2016, Eqs. 31–34), the same forms fit to BlackHawk output.
            </>
          }
        >
          <SpectrumExplorer />
        </Panel>

        <Prose>
          <p>
            The mass-loss rate integrates all of this over every available species:
          </p>
          <TeX block math="\frac{dM}{dt} = -\frac{\alpha(M)}{M^2} \quad\Rightarrow\quad \tau_{\mathrm{BH}} \sim \frac{M^3}{3\alpha}" />
          <p>
            This produces the strangest lightcurve in astrophysics. Every other transient rises
            and fades. A black hole explosion <strong>only ever brightens</strong>, following a
            universal curve fixed by known physics <Cite k="jcap" sec="§2" /> — right up to the
            instant the hole ceases to exist.
          </p>
        </Prose>

        <Panel
          title="the backwards burst"
          caption="Total photon output per second versus (inverse) remaining lifetime. Scrub toward the explosion: in the last second, output climbs five orders of magnitude."
        >
          <BackwardsBurst />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 02 ORIGINS */}
      <Section
        id="origins"
        index="02"
        kicker="where small black holes come from"
        title="Forged in the first second"
      >
        <Prose>
          <p>
            No star can make a small black hole — stellar collapse bottoms out near three solar
            masses. But the infant universe could. In its first fraction of a second, space was
            filled with dense plasma; any patch with density contrast{" "}
            <TeX math="\delta \gtrsim 0.45" /> would collapse the moment it entered the causal
            horizon <Cite k="carr20" />. The mass swallowed is the mass inside the horizon at
            that instant:
          </p>
          <TeX block math="M_H(t) \approx 10^{15}\,\mathrm{g} \times \left(\frac{t}{10^{-23}\,\mathrm{s}}\right)" />
          <p>
            So the cosmic clock doubles as a mass dial: <em>when</em> a primordial black hole
            (PBH) formed sets <em>how big</em> it was born — and therefore <em>when it dies</em>.
          </p>
        </Prose>

        <Panel
          title="the formation clock"
          caption="Slide through cosmic time. Expansion continuously grows the horizon, so each epoch mints a characteristic PBH mass."
        >
          <FormationTimeline />
        </Panel>

        <Prose>
          <p>
            Nature would not mint a single mass but a distribution ψ(M), whose shape encodes the
            formation mechanism. The explosion rate we could hope to see today is set by how much
            of that distribution sits at the critical mass{" "}
            <TeX math="M_U \simeq 5.1\times10^{14}\,\mathrm{g}" /> — holes whose lifetime equals
            the age of the universe <Cite k="jcap" sec="Eq. 3.8" />:
          </p>
          <TeX block math="\dot{n}_{\mathrm{PBH}} = \rho_{\mathrm{DM}}\,\frac{\psi_i(M_U)}{3\,t_U}" />
        </Prose>

        <Panel
          title="mass function gallery"
          caption={
            <>
              Three physically motivated families, normalized and evaluated by the{" "}
              <RepoLink path={REPO_PATHS.dataPipeline} file>
                data pipeline
              </RepoLink>{" "}
              (Boluna et al. Eqs. 3.9–3.14). The red line marks M_U: only the sliver of ψ crossing
              it explodes on our watch.
            </>
          }
        >
          <MassFunctions />
        </Panel>

        <Note>
          A subtlety worth keeping straight: cosmic expansion enables PBH{" "}
          <strong>formation</strong> by setting the horizon scale in the early universe. Holes
          are not being created by expansion today — the ones exploding now are fossils from
          t ≈ 10⁻²³ s, finishing a 13.8-billion-year evaporation.
        </Note>
      </Section>

      {/* ------------------------------------------------ 03 COLLIDER */}
      <Section
        id="collider"
        index="03"
        kicker="dark degrees of freedom"
        title="A collider built by gravity"
      >
        <Prose>
          <p>
            Here is the property that makes exploding black holes more than a curiosity: Hawking
            radiation is <strong>democratic</strong>. Gravity couples to everything, so a black
            hole radiates every particle species lighter than its temperature — quarks, neutrinos,
            gravitons, and <em>anything else that exists</em>, including particles that never
            touch our detectors <Cite k="thesis" sec="§2.8" />.
          </p>
          <p>
            The evaporation coefficient α(M) literally counts the particle content of nature. As
            the hole shrinks and heats, each new species unlocks like a threshold in a collider
            energy scan:
          </p>
        </Prose>

        <Panel
          title="the degrees-of-freedom staircase"
          caption="α(M) rises each time k_BT crosses a particle mass (thresholds marked). Toggle the dark sector: doubling the degrees of freedom roughly halves the lifetime and dilutes the photon share of the luminosity — an observable, falsifiable shift."
        >
          <AlphaStaircase />
        </Panel>

        <Prose>
          <p>
            In its final seconds a PBH reaches temperatures beyond any accelerator — a
            <em> cosmic-scale collider</em> whose luminosity, spectrum, and duration all depend
            on the full particle spectrum of nature. A dark sector with N copies of the Standard
            Model would shorten the final TeV burst by ~N and dim its photons by the same
            factor <Cite k="thesis" sec="§2.8" />. Watching one explosion is a census of
            everything that can exist.
          </p>
          <p>
            Even without new physics, evaporation carries a unique fingerprint. The measured
            power-law slope of the spectrum evolves as the hole shrinks, converging to a
            universal value that no astrophysical source reproduces:
          </p>
        </Prose>

        <Panel
          title="the universal spectral index"
          caption={
            <>
              Regenerated from the{" "}
              <RepoLink path={REPO_PATHS.dataPipeline} file>
                data pipeline
              </RepoLink>{" "}
              (Boluna et al. Eq. 4.8, cf. Fig. 6). At 1 GeV the slope sweeps through γ ≈ 2–3
              (quasar-like territory) before locking onto γ → 1.5 near death.
            </>
          }
        >
          <SpectralIndex />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 04 TELESCOPES */}
      <Section
        id="telescopes"
        index="04"
        kicker="the instrument fleet"
        title="Six eyes on the γ-ray sky"
      >
        <Prose>
          <p>
            No single telescope covers the evaporation story. The spectrum spans nine decades of
            energy as the hole heats from MeV to TeV, so the search is inherently{" "}
            <strong>multi-mission</strong> <Cite k="thesis" sec="§2.7" />: space telescopes
            (Fermi's GBM and LAT) own the keV–GeV band with huge fields of view; ground arrays
            (VERITAS, HAWC, LHAASO) own the TeV band with vast effective areas.
          </p>
        </Prose>

        <Panel
          title="effective area atlas"
          caption={
            <>
              Digitized instrument response curves from{" "}
              <RepoLink path={REPO_PATHS.effectiveAreas}>EffectiveAreas/*.dat</RepoLink> in the
              BHRad repo (cf. Boluna et al. Fig. 3). Overlay the PBH spectrum and slide its mass:
              watch which instrument &apos;owns&apos; each stage of evaporation.
            </>
          }
        >
          <DetectorAtlas />
        </Panel>

        <Prose>
          <p>
            Combining instruments honestly is its own craft. The analysis behind this site uses{" "}
            <strong>threeML</strong> <Cite k="threeml" />, the multi-mission maximum-likelihood
            framework: every instrument keeps its own response and background model, while a
            single physical source model is fit jointly across all of them{" "}
            <Cite k="thesis" sec="Fig. 2.7" />.
          </p>
        </Prose>

        <Panel
          title="anatomy of a joint fit"
          caption={
            <>
              The{" "}
              <RepoLink path={REPO_PATHS.fittingPipeline} file>
                threeML fitting pipeline
              </RepoLink>{" "}
              used for every candidate in §06 — click through each stage.
            </>
          }
        >
          <ThreeMLPipeline />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 05 DISTANCE */}
      <Section
        id="distance"
        index="05"
        kicker="the brutal arithmetic of flux"
        title="Only the nearest explosions count"
      >
        <Prose>
          <p>
            Here is the sobering part. A PBH explosion is intrinsically <em>identical</em> every
            time — same mass, same lightcurve, same luminosity. That makes the detection
            criterion pure geometry <Cite k="jcap" sec="§4" />:
          </p>
          <TeX block math="N_S = \int\!\!\int \frac{d^2N}{dE\,dt}\,\frac{A_{\mathrm{eff}}(E)}{4\pi d^2}\,dE\,dt \;\geq\; 10 \quad\text{or}\quad \frac{N_S}{\sqrt{N_B}} \geq 5" />
          <p>
            Run that requirement through every instrument and you get a{" "}
            <strong>sensitivity frontier</strong>: the farthest distance each telescope could see
            a hole of a given mass.
          </p>
        </Prose>

        <Panel
          title="the sensitivity frontier"
          caption={
            <>
              Regenerated end-to-end by this site&apos;s{" "}
              <RepoLink path={REPO_PATHS.dataPipeline} file>
                data pipeline
              </RepoLink>{" "}
              from{" "}
              <RepoLink path={REPO_PATHS.analyticalModelling} file>
                Analytical_Modelling.ipynb
              </RepoLink>{" "}
              detectability code — same effective areas, backgrounds, and detection criteria (cf.
              Boluna et al. Fig. 4). Everything below a curve is visible to that instrument.
            </>
          }
        >
          <DistanceFrontier />
        </Panel>

        <Prose>
          <p>
            The peak reach is a fraction of a parsec — thousands of times closer than the nearest
            star. And nearness cuts both ways: anything that close should visibly{" "}
            <em>move</em>. At galactic speeds (~220 km/s), a source at 0.01 pc drifts about a
            degree per year <Cite k="jcap" sec="Eq. 4.9" /> — a smoking-gun streak, but also a
            reason catalogs might discard the very sources we want.
          </p>
        </Prose>

        <Panel
          title="our cosmic backyard"
          caption="Log-scale map from Earth outward. The detection horizons sit deep inside the Oort cloud's neighborhood — set the search radius and see how many exploding holes the abundance limits allow inside it."
        >
          <Neighborhood />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 06 FITTING */}
      <Section
        id="fitting"
        index="06"
        kicker="template vs. reality"
        title="Fitting the flash"
      >
        <Prose>
          <p>
            Suppose a short burst trips the GBM. Is it a neutron-star merger at a billion
            parsecs, or a black hole dying at a thousandth of one? The lightcurve holds the
            answer. The PBH template rises as <TeX math="(\tau - t)^{-0.52}" /> toward the
            explosion epoch τ — time-reversed compared to every conventional burst — optionally
            followed by an afterglow from the ejecta shell <Cite k="thesis" sec="§2.9–2.10" />.
          </p>
        </Prose>

        <Panel
          title="fit real bursts by hand"
          caption={
            <>
              Background-subtracted 100 ms lightcurves for seven candidate GRBs, exported directly
              from{" "}
              <RepoLink path={REPO_PATHS.lightcurveData}>
                Lightcurve_Fitting/~100ms_Source_Data
              </RepoLink>{" "}
              in the BHRad repo. Drive the template parameters the Bayesian sampler explores.
            </>
          }
        >
          <LightcurveComposer />
        </Panel>

        <Prose>
          <p>
            In the production pipeline these fits run through threeML with ultranest's nested
            sampling — hundreds of live points exploring normalization, decay index, and
            afterglow timing against every GBM detector simultaneously, returning posterior
            distributions and Bayesian evidence for template comparison.
          </p>
        </Prose>
      </Section>

      {/* ------------------------------------------------ 07 CATALOGUE */}
      <Section
        id="catalogue"
        index="07"
        kicker="mining fourteen years of Fermi data"
        title="The candidate hunt"
      >
        <Prose>
          <p>
            The Fermi GBM has logged thousands of bursts since 2008. Almost all are ordinary.
            The search strategy <Cite k="thesis" sec="§3.2" /> applies three cuts inherited from
            the physics: the burst must be <strong>short</strong> (T90 within 0.2–5 s),{" "}
            <strong>hard</strong> (more high-energy fluence than low), and{" "}
            <strong>local</strong> — no measured redshift, since a real PBH burst comes from
            inside our stellar neighborhood.
          </p>
          <p>
            The idea traces back to Cline's BATSE analyses in the 1990s{" "}
            <Cite k="cline97" sec="§3" />, which found a curious subpopulation of very short,
            anomalously hard bursts. The{" "}
            <RepoLink path={REPO_PATHS.catalogSearch} file>
              repo&apos;s Fermi catalog search
            </RepoLink>{" "}
            modernizes this on Fermi&apos;s much deeper catalog.
          </p>
        </Prose>

        <Panel
          title="filter playground — real candidates"
          caption={
            <>
              The 36 sources that survived the repo&apos;s full selection (
              <RepoLink path={REPO_PATHS.candidatesCsv} file>
                H&gt;1_T90[0.2-5]_RS=0.csv
              </RepoLink>
              ), with live re-filtering. Hardness here is LAT fluence over GBM fluence; Cline&apos;s
              original hardness–duration anticorrelation used BATSE bands.
            </>
          }
        >
          <CataloguePlayground />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 08 SKY */}
      <Section
        id="sky"
        index="08"
        kicker="the isotropy test"
        title="Where do they point?"
      >
        <Prose>
          <p>
            Geometry offers a free hypothesis test. Neutron stars, magnetars, X-ray binaries —
            everything born of stars traces the Milky Way's disk. But local PBH explosions
            sample only our tiny corner of the halo, so their sky map should be{" "}
            <strong>isotropic</strong>, indifferent to the galactic plane{" "}
            <Cite k="cline97" sec="Fig. 4" /> <Cite k="thesis" sec="§2.3" />.
          </p>
        </Prose>

        <Panel
          title="the candidate sky, in galactic coordinates"
          caption={
            <>
              Mollweide projection with the galactic center at the middle and the plane along the
              equator — the same visualization produced by{" "}
              <RepoLink path={REPO_PATHS.boresightSelection} file>
                BoresightSelection.py
              </RepoLink>{" "}
              in the BHRad repo. Equatorial coordinates converted to galactic on the fly.
            </>
          }
        >
          <SkyMap />
        </Panel>

        <Prose>
          <p>
            Cline additionally reported an anomalous cluster of very short bursts in one octant
            of the sky <Cite k="cline02" /> — never confirmed, never fully refuted. A larger
            candidate sample with Fermi-era statistics is exactly what could settle it.
          </p>
        </Prose>
      </Section>

      {/* ------------------------------------------------ 09 AFTERGLOW */}
      <Section
        id="afterglow"
        index="09"
        kicker="the multi-messenger encore"
        title="Echoes weeks after the flash"
      >
        <Prose>
          <p>
            The γ-ray flash may not be the last word. The explosion dumps ~10²⁵ J of
            electron-positron pairs into the interstellar medium; Rees noted in 1977 that this
            conducting fireball plows the ambient magnetic field into a coherent low-frequency
            radio pulse <Cite k="rees77" /> <Cite k="blandford77" /> — potentially detectable
            across the entire galaxy, far beyond the γ-ray horizon.
          </p>
          <p>
            And blazar physics adds a slower channel: adiabatically expanding ejecta re-radiate
            at progressively lower frequencies, with radio peaking 40–140 days after the γ-ray
            flare in observed jets. Would a spherical PBH shell do the same? That is an open
            question the{" "}
            <RepoLink path={REPO_PATHS.afterglowNotes} file>
              repo&apos;s afterglow notes
            </RepoLink>{" "}
            flag explicitly <Cite k="thesis" sec="§2.10" />.
          </p>
        </Prose>

        <Panel
          title="one explosion, two messengers"
          caption="Click along the timeline. The γ-ray burst and the radio afterglow are separated by up to nine orders of magnitude in time — a coincidence search across archives is the discovery strategy."
        >
          <AfterglowTimeline />
        </Panel>
      </Section>

      {/* ------------------------------------------------ 10 TRANSIENTS */}
      <Section
        id="transients"
        index="10"
        kicker="the slow-burn search"
        title="Sources that only brighten"
      >
        <Prose>
          <p>
            There's a second way to catch a dying black hole: <em>before</em> the explosion.
            For months to years, a nearby PBH would appear as a faint GeV point source with a
            steadily rising flux <Cite k="jcap" sec="Eq. 5.3" />:
          </p>
          <TeX block math="F_\gamma(t) \simeq 2.7\times10^{-8}\left(\frac{\mathrm{pc}}{d}\right)^{2}\left(\frac{\tau - t}{\mathrm{s}}\right)^{-0.533}\,\mathrm{cm^{-2}\,s^{-1}}" />
          <p>
            The Fermi LAT Transient Catalog lists 35 sources with no known counterpart at any
            wavelength. The{" "}
            <RepoLink path={REPO_PATHS.transientFits} file>
              thesis transient-fitting pipeline
            </RepoLink>{" "}
            fit every one of them with this two-parameter model{" "}
            <Cite k="thesis" sec="§3.1" />.
          </p>
        </Prose>

        <Panel
          title="fit a rising transient"
          caption="A mock decade-long LAT lightcurve with realistic scatter. Trade τ against d and feel the degeneracy that dominates the real fits — then check the frontier plot in §05 with the real fitted sources overlaid."
        >
          <TransientFitter />
        </Panel>

        <Prose>
          <p>
            The verdict so far: the fits succeed — two clusters of solutions, one near 10¹² g at
            milliparsec distances, one near 10¹⁴⁻¹⁵ g further out — but both imply proper motion
            of degrees per year that the LAT catalog does not observe. The simplest reading is
            that these transients are something else. The method, though, now exists and sharpens
            with every year of data.
          </p>
        </Prose>
      </Section>

      {/* ------------------------------------------------ 11 VERDICT */}
      <Section
        id="verdict"
        index="11"
        kicker="scientific honesty"
        title="Where the evidence stands"
      >
        <Prose>
          <p>
            No exploding black hole has been found. This site would be dishonest to imply
            otherwise — and dishonest to imply the search is hopeless. Here is the current
            ledger, claim by claim. Click to expand.
          </p>
        </Prose>

        <Panel title="the evidence board">
          <EvidenceBoard />
        </Panel>

        <Prose>
          <p>
            The deeper reason to keep looking: a single confirmed detection would
            simultaneously prove black holes evaporate (quantum gravity's only accessible
            prediction), demonstrate dark matter physics beyond WIMPs, and census every particle
            degree of freedom in nature — visible or dark. Few observations in physics carry
            that much payload.
          </p>
          <p>
            The tools are improving on schedule: CTA will push the sensitivity frontier past
            every current instrument, LAT keeps accumulating exposure, and the template methods
            built in this research — universal lightcurves, spectral-index screening, transient
            fitting, proper-motion vetoes — are ready for the data.
          </p>
        </Prose>

        <footer className="mt-20 border-t border-line pt-8 pb-4">
          <div className="grid gap-6 sm:grid-cols-2">
            <div>
              <div className="font-mono text-[10px] uppercase tracking-widest text-ink-faint">
                primary sources
              </div>
              <ul className="mt-2 space-y-1 text-[12px] text-ink-dim">
                <li>
                  <a
                    href="https://arxiv.org/abs/2307.06467"
                    className="text-signal hover:underline"
                    target="_blank"
                    rel="noreferrer"
                  >
                    Boluna, Profumo, Blé & Hennings — Searching for Exploding Black Holes (JCAP)
                  </a>
                </li>
                <li>
                  <a
                    href="https://www.xboluna.com/media/documents/xboluna_UCSC_thesis.pdf"
                    className="text-signal hover:underline"
                    target="_blank"
                    rel="noreferrer"
                  >
                    Boluna — Detection Methods for Discovering Evaporating PBHs (UCSC thesis)
                  </a>
                </li>
                <li>Cline, Sanders & Hong — ApJ 486, 169 (1997)</li>
                <li>Ukwatta et al. — Astropart. Phys. 80, 90 (2016)</li>
              </ul>
            </div>
            <div>
              <div className="font-mono text-[10px] uppercase tracking-widest text-ink-faint">
                about this site
              </div>
              <p className="mt-2 text-[12px] leading-relaxed text-ink-dim">
                Every chart is generated from the{" "}
                <RepoLink path={REPO_PATHS.bhrad}>BHRad research repository</RepoLink> — the same
                effective areas, spectra parameterizations, candidate catalogs, and fitted parameters
                used in the published analysis. Regenerated by the{" "}
                <RepoLink path={REPO_PATHS.dataPipeline} file>
                  website data pipeline
                </RepoLink>
                .
              </p>
              <a
                href="https://github.com/xboluna/Research"
                target="_blank"
                rel="noreferrer"
                className="mt-4 inline-flex w-fit items-center gap-2 rounded border border-signal/50 bg-signal/10 px-4 py-2 font-mono text-[11px] uppercase tracking-widest text-signal transition-colors hover:bg-signal/20"
              >
                view on GitHub ↗
              </a>
            </div>
          </div>
        </footer>
      </Section>
    </main>
  );
}
