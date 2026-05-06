#!/usr/bin/env python3
"""
parse_gapseq_report.py
Parse one gapseq fill stdout into a single TSV row.

Usage:
    python3 parse_gapseq_report.py <organism_id> <report.txt> <out_tsv>

First call creates the file with header. Subsequent calls append rows.
"""
import sys, re
from pathlib import Path

def grab(pattern, text, group=1, default="NA"):
    m = re.search(pattern, text)
    return m.group(group).strip() if m else default

def section(text, start_pat, end_pat):
    s = re.search(start_pat, text)
    if not s: return ""
    e = re.search(end_pat, text[s.end():])
    return text[s.end():s.end() + e.start()] if e else text[s.end():]

def parse(report_path):
    text = Path(report_path).read_text()

    orf_cov = grab(r"ORF coverage:\s+([\d.]+)\s*%", text)
    n_genes = grab(r"Creating Gene-Reaction list\.\.\.\s+(\d+)\s+unique genes", text)

    s1 = section(text, r"1\.\s*Initial gapfilling", r"2\.\s*Biomass gapfilling")
    s1_added = grab(r"Added reactions:\s+(\d+)", s1)
    s1_core  = grab(r"Added core reactions:\s+(\d+)", s1)
    s1_grow  = grab(r"Final growth rate:\s+([\d.]+)", s1)

    s3 = section(text, r"3\.\s*Energy source gapfilling", r"4\.\s*Checking")
    s3_filled    = grab(r"Filled components:\s+(\d+)", s3)
    s3_added     = grab(r"Added reactions:\s+(\d+)", s3)
    s3_grow      = grab(r"Final growth rate:\s+([\d.]+)", s3)
    s3_components = grab(r"Filled components:\s+\d+\s+\(\s*(.+?)\s*\)", s3)

    s4 = section(text, r"4\.\s*Checking for potential", r"Uptake at limit")
    s4_filled = grab(r"Filled components:\s+(\d+)", s4)
    s4_added  = grab(r"Added reactions:\s+(\d+)", s4)
    s4_grow   = grab(r"Final growth rate:\s+([\d.]+)", s4)

    uptake   = grab(r"Uptake at limit:\s*\n([^\n]+)", text)
    products = grab(r"Top 10 produced metabolites[^\n]*\n([^\n]+)", text)

    return {
        "orf_coverage_pct": orf_cov,
        "n_unique_genes":   n_genes,
        "s1_added":         s1_added,
        "s1_core":          s1_core,
        "s1_growth":        s1_grow,
        "s3_filled":        s3_filled,
        "s3_added":         s3_added,
        "s3_growth":        s3_grow,
        "s3_components":    s3_components,
        "s4_filled":        s4_filled,
        "s4_added":         s4_added,
        "final_growth":     s4_grow,
        "uptake_at_limit":  uptake,
        "top_products":     products,
    }

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: parse_gapseq_report.py <organism_id> <report.txt> <out_tsv>")
    org_id, report_path, out_tsv = sys.argv[1:]

    fields = ["organism", "orf_coverage_pct", "n_unique_genes",
              "s1_added", "s1_core", "s1_growth",
              "s3_filled", "s3_added", "s3_growth", "s3_components",
              "s4_filled", "s4_added", "final_growth",
              "uptake_at_limit", "top_products"]

    data = parse(report_path)
    data["organism"] = org_id

    out = Path(out_tsv)
    write_header = not out.exists()
    with out.open("a") as f:
        if write_header:
            f.write("\t".join(fields) + "\n")
        f.write("\t".join(data[k] for k in fields) + "\n")
    print(f"Parsed {org_id} -> {out_tsv}")

if __name__ == "__main__":
    main()
