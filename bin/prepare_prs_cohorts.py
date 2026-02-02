#!/usr/bin/env python3
"""
prepare_prs_cohorts.py

Builds cohort TSVs for PRS→PheWAS by:
  - reading PRS matrix (FID/IID + selected PRS trait columns)
  - applying Blom rank-based inverse normal transform (INT) to each selected PRS column
  - merging INT PRS columns into an existing cohort-with-covariates TSV
  - writing:
      * traits_manifest.tsv
      * cohort_with_covariates.prs.pooled.tsv
      * optionally cohort_with_covariates.prs.<ancestry>.tsv (filtered by ancestry_select)

Inputs:
  --prs-matrix: TSV with columns FID, IID, <traits...>
  --cohort-covariates: TSV with person_id, age*, sex_at_birth, pcs, genetic_ancestry, ...
"""

import argparse
import re
from pathlib import Path
from typing import List, Optional, Set

import numpy as np
import pandas as pd
from scipy.stats import norm


def parse_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    s = str(x).strip().lower()
    return s in ("1", "true", "t", "yes", "y", "on")


def safe_tag(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name))


def rank_int(series: pd.Series, rng: np.random.Generator, jitter_sd: float) -> pd.Series:
    """Blom rank-based inverse normal transform with optional tiny jitter."""
    x = pd.to_numeric(series, errors="coerce").astype(float)
    mask = x.notna()

    x2 = x.copy()
    if jitter_sd and jitter_sd > 0:
        x2.loc[mask] = x2.loc[mask] + rng.normal(0.0, float(jitter_sd), size=int(mask.sum()))

    ranks = x2.loc[mask].rank(method="average")
    n = int(mask.sum())
    if n == 0:
        return pd.Series(np.nan, index=series.index, dtype=float)

    p = (ranks - 0.375) / (n + 0.25)
    p = p.clip(1e-12, 1 - 1e-12)

    out = pd.Series(np.nan, index=series.index, dtype=float)
    out.loc[mask] = norm.ppf(p)
    return out


def resolve_traits(prs_path: Path, prs_traits_arg: str) -> List[str]:
    """
    Resolve trait list from:
      - "ALL"
      - comma-separated list
      - file path with one trait per line (comments allowed with #)
    """
    hdr = pd.read_csv(prs_path, sep="\t", nrows=0)
    cols = hdr.columns.tolist()

    id_cols = {"FID", "IID"}
    trait_cols = [c for c in cols if c not in id_cols]

    arg = (prs_traits_arg or "").strip()
    if not arg:
        return []

    if arg.upper() == "ALL":
        return trait_cols

    p = Path(arg)
    if p.exists() and p.is_file():
        traits = [
            ln.strip()
            for ln in p.read_text().splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        return traits

    traits = [t.strip() for t in arg.split(",") if t.strip()]
    return traits


def parse_ancestry_select(sel_raw: str) -> str:
    """
    Normalize ancestry selection:
      - ALL
      - POOLED
      - LIST (comma list)
    """
    s = (sel_raw or "ALL").strip()
    su = s.upper()
    if su in ("ALL",):
        return "ALL"
    if su in ("POOLED", "POOL", "POOLED_ONLY", "POOLED-ONLY", "POOLEDONLY"):
        return "POOLED"
    return "LIST"


def parse_selected_tags(sel_raw: str) -> Set[str]:
    """For LIST mode: return normalized lower-case tags, excluding pooled."""
    tags = set()
    for tok in (sel_raw or "").split(","):
        t = tok.strip()
        if not t:
            continue
        t = safe_tag(t).lower()
        tags.add(t)
    tags.discard("pooled")
    return tags


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prs-matrix", required=True, help="PRS matrix TSV (FID/IID + trait cols)")
    ap.add_argument("--cohort-covariates", required=True, help="Cohort covariates TSV (person_id + covariates)")
    ap.add_argument("--prs-traits", required=True, help="ALL | comma list | file path (one per line)")
    ap.add_argument("--out-prefix", default="cohort_with_covariates.prs")
    ap.add_argument("--stratify-by-ancestry", default="true")
    ap.add_argument("--ancestry-col", default="genetic_ancestry")
    ap.add_argument("--person-id-col", default="person_id")

    # NEW: choose which ancestries to output
    ap.add_argument(
        "--ancestry-select",
        default="ALL",
        help="ALL | POOLED | comma list (e.g. eur,afr). Case-insensitive. "
             "POOLED writes pooled cohort only.",
    )

    # INT settings
    ap.add_argument("--int-seed", type=int, default=1)
    ap.add_argument("--int-jitter-sd", type=float, default=1e-12)

    # logging
    ap.add_argument("--log", default="prep.log")
    args = ap.parse_args()

    prs_path = Path(args.prs_matrix)
    cov_path = Path(args.cohort_covariates)

    rng = np.random.default_rng(int(args.int_seed))
    logp = Path(args.log)

    def log(msg: str):
        msg = str(msg)
        if logp.exists():
            logp.write_text(logp.read_text() + msg + "\n")
        else:
            logp.write_text(msg + "\n")

    log("=== PREPARE_PRS_COHORTS ===")
    log(f"PRS matrix: {prs_path}")
    log(f"Cohort covariates: {cov_path}")
    log(f"prs_traits: {args.prs_traits}")
    log(f"stratify_by_ancestry: {args.stratify_by_ancestry}")
    log(f"ancestry_col: {args.ancestry_col}")
    log(f"person_id_col: {args.person_id_col}")
    log(f"ancestry_select: {args.ancestry_select}")
    log(f"int_seed: {args.int_seed}")
    log(f"int_jitter_sd: {args.int_jitter_sd}")

    if not prs_path.exists():
        raise SystemExit(f"PRS matrix not found: {prs_path}")
    if not cov_path.exists():
        raise SystemExit(f"Cohort covariates not found: {cov_path}")

    # --- traits ---
    traits = resolve_traits(prs_path, args.prs_traits)
    if not traits:
        raise SystemExit("No PRS traits resolved from --prs-traits (empty).")

    hdr = pd.read_csv(prs_path, sep="\t", nrows=0)
    cols = hdr.columns.tolist()
    missing_traits = [t for t in traits if t not in cols]
    if missing_traits:
        raise SystemExit(
            f"These PRS traits were not found in PRS matrix header: "
            f"{missing_traits[:20]}{'...' if len(missing_traits) > 20 else ''}"
        )

    # --- load PRS subset (IID + selected traits) ---
    usecols = ["IID"] + traits
    prs = pd.read_csv(prs_path, sep="\t", usecols=usecols)
    prs = prs.rename(columns={"IID": args.person_id_col})
    prs[args.person_id_col] = prs[args.person_id_col].astype("string")

    # --- INT each trait and build manifest ---
    manifest_rows = []
    int_cols = []
    for t in traits:
        ttag = safe_tag(t)
        iv_col = f"int_{ttag}"
        prs[iv_col] = rank_int(prs[t], rng=rng, jitter_sd=float(args.int_jitter_sd))
        manifest_rows.append({"trait": t, "trait_tag": ttag, "iv_col": iv_col})
        int_cols.append(iv_col)

    prs_int = prs[[args.person_id_col] + int_cols].copy()

    # --- load cohort covariates ---
    cohort = pd.read_csv(cov_path, sep="\t", dtype={args.person_id_col: "string"})

    # detect age column
    if "age_at_last_ehr_event" in cohort.columns:
        age_col = "age_at_last_ehr_event"
    elif "age_at_last_event" in cohort.columns:
        age_col = "age_at_last_event"
    else:
        raise SystemExit(
            "Could not find age column in cohort covariates "
            "(expected age_at_last_ehr_event or age_at_last_event)."
        )

    pc_cols = [f"pc{i}" for i in range(1, 11)]
    required = [args.person_id_col, age_col, "sex_at_birth"] + pc_cols

    stratify = parse_bool(args.stratify_by_ancestry)
    if stratify:
        required.append(args.ancestry_col)

    missing_cov = [c for c in required if c not in cohort.columns]
    if missing_cov:
        raise SystemExit(f"Missing required covariate columns in cohort file: {missing_cov}")

    cohort_base = cohort[required].copy()

    # --- merge ---
    merged = cohort_base.merge(prs_int, on=args.person_id_col, how="inner")
    merged = merged.dropna(subset=required + int_cols)

    log(f"After merge+dropna: N={merged.shape[0]}")

    # --- write outputs ---
    out_prefix = args.out_prefix

    # manifest
    pd.DataFrame(manifest_rows).to_csv("traits_manifest.tsv", sep="\t", index=False)
    log("Wrote traits_manifest.tsv")

    # pooled always
    pooled_path = f"{out_prefix}.pooled.tsv"
    merged.to_csv(pooled_path, sep="\t", index=False)
    log(f"Wrote pooled cohort: {pooled_path}")

    # ancestry selection mode
    mode = parse_ancestry_select(args.ancestry_select)
    if not stratify:
        log("Stratify disabled -> skipping ancestry cohorts.")
        log("Done.")
        return

    if mode == "POOLED":
        log("ancestry_select=POOLED -> skipping ancestry cohorts.")
        log("Done.")
        return

    selected_tags: Optional[Set[str]] = None
    if mode == "LIST":
        selected_tags = parse_selected_tags(args.ancestry_select)
        if not selected_tags:
            # User gave empty/pooled-only-like list; treat as pooled
            log(f"ancestry_select list resolved empty -> skipping ancestry cohorts.")
            log("Done.")
            return
        log(f"ancestry_select tags: {sorted(selected_tags)}")

    ancestries = sorted(merged[args.ancestry_col].astype(str).unique())
    log(f"Ancestries observed: {ancestries}")

    wrote_any = 0
    for anc in ancestries:
        tag = safe_tag(anc).lower()

        if selected_tags is not None and tag not in selected_tags:
            continue

        out_path = f"{out_prefix}.{tag}.tsv"
        merged.loc[merged[args.ancestry_col].astype(str) == anc].to_csv(out_path, sep="\t", index=False)
        wrote_any += 1
        log(f"Wrote ancestry cohort: {out_path}")

    if selected_tags is not None and wrote_any == 0:
        log("WARNING: ancestry_select was a list but none matched observed ancestries (check naming).")

    log("Done.")


if __name__ == "__main__":
    main()

