import katex from "katex";

/** Server-rendered KaTeX. Use $-free TeX source. */
export function TeX({ math, block = false }: { math: string; block?: boolean }) {
  const html = katex.renderToString(math, {
    displayMode: block,
    throwOnError: false,
    strict: false,
  });
  return (
    <span
      className={block ? "block my-3" : "inline"}
      dangerouslySetInnerHTML={{ __html: html }}
    />
  );
}
