import type { ReactNode } from "react";

export function Section({
  id,
  index,
  kicker,
  title,
  children,
}: {
  id: string;
  index: string;
  kicker: string;
  title: string;
  children: ReactNode;
}) {
  return (
    <section id={id} className="relative mx-auto max-w-5xl px-5 py-24 sm:px-8">
      <header className="mb-10">
        <div className="flex items-baseline gap-3 font-mono text-xs tracking-[0.2em] text-ink-faint uppercase">
          <span className="text-signal">{index}</span>
          <span>{kicker}</span>
        </div>
        <h2 className="mt-3 text-3xl font-medium tracking-tight text-ink sm:text-4xl">
          {title}
        </h2>
      </header>
      {children}
    </section>
  );
}

export function Prose({ children }: { children: ReactNode }) {
  return (
    <div className="max-w-3xl space-y-4 text-[15px] leading-relaxed text-ink-dim [&_strong]:text-ink [&_em]:text-ink/90">
      {children}
    </div>
  );
}

export function Panel({
  title,
  caption,
  children,
}: {
  title?: string;
  caption?: ReactNode;
  children: ReactNode;
}) {
  return (
    <figure className="my-10 overflow-hidden rounded-lg border border-line bg-surface">
      {title && (
        <figcaption className="border-b border-line px-4 py-2.5 font-mono text-[11px] uppercase tracking-[0.15em] text-ink-faint">
          {title}
        </figcaption>
      )}
      <div className="p-4 sm:p-5">{children}</div>
      {caption && (
        <figcaption className="border-t border-line px-4 py-3 text-xs leading-relaxed text-ink-faint">
          {caption}
        </figcaption>
      )}
    </figure>
  );
}

export function Note({ children }: { children: ReactNode }) {
  return (
    <aside className="my-6 max-w-3xl rounded border-l-2 border-theory/60 bg-theory/5 px-4 py-3 text-sm leading-relaxed text-ink-dim">
      {children}
    </aside>
  );
}
