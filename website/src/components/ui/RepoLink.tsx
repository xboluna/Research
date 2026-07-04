import type { ReactNode } from "react";
import { REPO_ROOT, repoFileUrl, repoUrl } from "@/lib/repo";

type RepoLinkProps = {
  /** Path inside the repo, e.g. `BHRad/Modelling/Analytical_Modelling.ipynb`. Omit for repo root. */
  path?: string;
  /** Use `/blob/` instead of `/tree/` — better for direct file links. */
  file?: boolean;
  children: ReactNode;
  className?: string;
};

export function RepoLink({
  path,
  file = false,
  children,
  className = "text-signal hover:underline",
}: RepoLinkProps) {
  const href = path ? (file ? repoFileUrl(path) : repoUrl(path)) : REPO_ROOT;

  return (
    <a href={href} target="_blank" rel="noreferrer" className={className}>
      {children}
    </a>
  );
}
