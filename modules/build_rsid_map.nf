// ================================================================
//  modules/build_rsid_map.nf — Stage B1
// ================================================================
process BUILD_RSID_MAP {
    label 'small'
    publishDir "${params.out_dir}/refs", mode: 'copy'

    input:
    path bim_file

    output:
    path "rsid_map.tsv.gz", emit: rsid_map

    script:
    """
    python3 << 'PYEOF'
import gzip, csv

seen  = set()
n_out = 0
print(f"[RSID] Reading {bim_file} ...")

with gzip.open("rsid_map.tsv.gz", "wt") as fout:
    w = csv.writer(fout, delimiter="\\t")
    w.writerow(["CHR", "BP", "rsID"])
    with open("${bim_file}") as fin:
        for line in fin:
            p = line.strip().split()
            if len(p) < 4: continue
            chrom, rsid, _, bp = p[0], p[1], p[2], p[3]
            if not rsid.startswith("rs"): continue
            try: chrom, bp = int(chrom), int(bp)
            except: continue
            if (chrom, bp) in seen: continue
            seen.add((chrom, bp))
            w.writerow([chrom, bp, rsid])
            n_out += 1

print(f"[RSID] Wrote {n_out:,} entries to rsid_map.tsv.gz")
PYEOF
    """
}
