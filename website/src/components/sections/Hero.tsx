import { DitherHero } from "@/components/ui/DitherHero";
import { RepoLink } from "@/components/ui/RepoLink";
import { REPO_PATHS } from "@/lib/repo";

export function Hero() {
  return (
    <section id="hero" className="relative flex min-h-svh flex-col">
      <div className="absolute inset-0">
        <DitherHero className="h-full w-full" />
        <div className="absolute inset-0 bg-gradient-to-b from-transparent via-transparent to-void" />
      </div>

      <div className="relative z-10 mx-auto flex w-full max-w-5xl flex-1 flex-col justify-end px-5 pb-20 sm:px-8">
        <p className="font-mono text-[11px] uppercase tracking-[0.3em] text-signal">
          an interactive field guide
        </p>
        <h1 className="mt-4 max-w-3xl text-4xl font-medium leading-[1.08] tracking-tight text-ink sm:text-6xl">
          Searching for
          <br />
          Exploding Black Holes
        </h1>
        <p className="mt-6 max-w-xl rounded-lg border border-line/80 bg-void/85 px-5 py-4 text-[15px] leading-relaxed text-ink shadow-[0_12px_40px_rgba(0,0,0,0.55)] backdrop-blur-md sm:text-base">
          Half a century ago Stephen Hawking predicted that black holes evaporate — and that the
          smallest ones end their lives in an explosion brighter than anything gravity has built
          since. Somewhere in the data of our γ-ray telescopes, that flash may already be waiting.
        </p>
        <div className="mt-8 flex flex-wrap items-center gap-x-6 gap-y-2 font-mono text-[11px] text-ink-faint">
          <span>
            based on{" "}
            <a
              href="https://arxiv.org/abs/2307.06467"
              target="_blank"
              rel="noreferrer"
              className="text-signal hover:underline"
            >
              Boluna, Profumo, Blé & Hennings (JCAP)
            </a>
          </span>
          <span className="hidden sm:inline">·</span>
          <span>
            and the{" "}
            <RepoLink path={REPO_PATHS.bhrad}>UCSC BHRad research program</RepoLink>
          </span>
        </div>
        <a
          href="#hawking"
          className="mt-10 inline-flex w-fit items-center gap-2 rounded border border-signal/50 bg-signal/10 px-5 py-2.5 font-mono text-xs uppercase tracking-widest text-signal transition-colors hover:bg-signal/20"
        >
          begin ↓
        </a>
      </div>
    </section>
  );
}
