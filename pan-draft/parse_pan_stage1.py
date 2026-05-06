#!/usr/bin/env python3
"""parse_pan_stage1.py
Aggregate ORF coverage and unique gene counts from N stage-1 stdouts.
Usage: python3 parse_pan_stage1.py <pan_id> <log_dir> <out_tsv>
"""
import sys, re, statistics
from pathlib import Path

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: parse_pan_stage1.py <pan_id> <log_dir> <out_tsv>")
    pan_id, log_dir, out_tsv = sys.argv[1:]

    log_files = list(Path(log_dir).glob("stdout_step1_*.txt"))
    if not log_files:
        sys.exit(f"No stdout_step1_*.txt in {log_dir}")

    orf_covs, n_genes = [], []
    for f in log_files:
        text = f.read_text()
        cov = re.search(r"ORF coverage:\s+([\d.]+)", text)
        gen = re.search(r"Creating Gene-Reaction list\.\.\.\s+(\d+)\s+unique genes", text)
        if cov: orf_covs.append(float(cov.group(1)))
        if gen: n_genes.append(int(gen.group(1)))

    def stats(xs):
        if not xs: return ("NA",)*4
        return (f"{statistics.median(xs):.2f}", f"{min(xs):.2f}",
                f"{max(xs):.2f}", f"{statistics.mean(xs):.2f}")

    cov = stats(orf_covs); gen = stats(n_genes)
    fields = ["pan_id","n_mags_parsed",
              "orf_cov_median","orf_cov_min","orf_cov_max","orf_cov_mean",
              "n_genes_median","n_genes_min","n_genes_max","n_genes_mean"]

    out = Path(out_tsv); write_header = not out.exists()
    with out.open("a") as f:
        if write_header: f.write("\t".join(fields) + "\n")
        f.write("\t".join([pan_id, str(len(log_files)), *cov, *gen]) + "\n")
    print(f"Aggregated {len(log_files)} MAGs for {pan_id} -> {out_tsv}")

if __name__ == "__main__":
    main()
