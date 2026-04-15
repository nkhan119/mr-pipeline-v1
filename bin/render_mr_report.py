#!/usr/bin/env python3
"""
render_mr_report.py — MR Pipeline CDC 1.0.0  Stage B6
Reads all TSV results staged into the work directory, injects JSON
into the HTML template, and writes MR_Report.html.
"""

import argparse
import csv
import glob
import json
import os
from datetime import datetime
from pathlib import Path


def load_tsvs(pattern):
    rows = []
    for f in sorted(glob.glob(pattern)):
        try:
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    rows.append(row)
        except Exception as e:
            print(f"  [WARN] {f}: {e}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template",    required=True)
    ap.add_argument("--author",      default="")
    ap.add_argument("--affiliation", default="")
    ap.add_argument("--institute",   default="")
    ap.add_argument("--github",      default="")
    args = ap.parse_args()

    date_str = datetime.now().strftime("%B %d, %Y · %H:%M")

    mr_rows   = load_tsvs("*_mr.tsv")
    wald_rows = load_tsvs("*_wald.tsv")
    het_rows  = load_tsvs("*_het.tsv")
    loo_rows  = load_tsvs("*_loo.tsv")
    snp_rows  = load_tsvs("*_single.tsv")
    mvmr_rows = load_tsvs("*_mvmr.tsv")

    print(f"  uvMR estimates  : {len(mr_rows)}")
    print(f"  Instruments     : {len(wald_rows)}")
    print(f"  Het tests       : {len(het_rows)}")
    print(f"  LOO rows        : {len(loo_rows)}")
    print(f"  Single-SNP rows : {len(snp_rows)}")
    print(f"  MVMR estimates  : {len(mvmr_rows)}")

    tmpl = Path(args.template).read_text()

    replacements = {
        "__MR_DATA__":    json.dumps(mr_rows),
        "__WALD_DATA__":  json.dumps(wald_rows),
        "__HET_DATA__":   json.dumps(het_rows),
        "__LOO_DATA__":   json.dumps(loo_rows),
        "__SNP_DATA__":   json.dumps(snp_rows),
        "__MVMR_DATA__":  json.dumps(mvmr_rows),
        "__AUTHOR__":     args.author,
        "__AFFILIATION__":args.affiliation,
        "__INSTITUTE__":  args.institute,
        "__GITHUB__":     args.github,
        "__DATE__":       date_str,
    }
    for k, v in replacements.items():
        tmpl = tmpl.replace(k, v)

    Path("MR_Report.html").write_text(tmpl)
    print(f"\n  MR_Report.html written  ({len(tmpl)//1024} KB)")


if __name__ == "__main__":
    main()
