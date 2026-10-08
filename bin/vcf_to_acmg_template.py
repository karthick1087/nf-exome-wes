#!/usr/bin/env python3
"""
Generate an ACMG/AMP manual curation template from an annotated VCF.

Columns follow ACMG 2015 criteria groups:
  PVS1  PS1-4  PM1-6  PP1-5  (Pathogenic)
  BA1   BS1-4  BP1-7  (Benign)

Usage:
    python3 vcf_to_acmg_template.py \\
        --vcf sample.snpeff.vcf.gz \\
        --bcftools-vcf sample.fully_annotated.vcf.gz \\
        --sample SAMPLE001 \\
        --out acmg_curation.tsv
"""

import argparse
import gzip
from pathlib import Path

from snpeff_parser import parse_ann_entries, pick_ann_entry


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(s):
    d = {}
    for part in s.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def auto_flags(ann, info, zyg, vaf):
    """Suggest criteria flags for curator review (not final classification)."""
    flags = []
    impact = ann.get("IMPACT", "")
    annotation = ann.get("ANNOTATION", "")
    af = info.get("AF", info.get("AF_sas", ""))
    clnsig = info.get("CLNSIG", "")

    # PM2 — absent or extremely low in population
    try:
        if af and float(af) < 0.0001:
            flags.append("PM2_suggested")
    except ValueError:
        flags.append("PM2_suggested")

    # PVS1 — null variant (LOF)
    if any(x in annotation for x in ("frameshift", "stop_gained", "splice_donor", "splice_acceptor")):
        flags.append("PVS1_review")

    # PS3/PP3 — computational (manual confirmation needed)
    if "missense" in annotation:
        flags.append("PP3_review")

    # PM3 — detected in trans (needs parental testing)
    if zyg == "Heterozygous" and "recessive" in clnsig.lower():
        flags.append("PM3_needs_parental")

    # BA1 — common variant
    try:
        if af and float(af) > 0.05:
            flags.append("BA1_suggested")
    except ValueError:
        pass

    return ";".join(flags) if flags else "none_suggested"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bcftools-vcf", default="")
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    bcf_lookup = {}
    if args.bcftools_vcf and Path(args.bcftools_vcf).exists():
        with open_text(args.bcftools_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                c = line.strip().split("\t")
                bcf_lookup[(c[0], int(c[1]), c[3], c[4])] = parse_info(c[7])

    acmg_pathogenic = ["PVS1", "PS1", "PS2", "PS3", "PS4",
                       "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
                       "PP1", "PP2", "PP3", "PP4", "PP5"]
    acmg_benign = ["BA1", "BS1", "BS2", "BS3", "BS4",
                   "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]

    header = [
        "sample", "gene", "hgvs_c", "hgvs_p", "chrom", "pos",
        "zygosity", "vaf_pct", "impact", "clinvar_sig", "gnomad_af",
        "auto_suggested_flags",
        "final_classification",   # Pathogenic / Likely pathogenic / VUS / Benign
        "curator", "review_date", "notes",
    ] + acmg_pathogenic + acmg_benign

    rows = []
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
            info = parse_info(cols[7])
            ann = pick_ann_entry(parse_ann_entries(info.get("ANN", "")), alt=alt)

            fmt = cols[8].split(":")
            sval = cols[9].split(":")
            fmap = dict(zip(fmt, sval))
            gt = fmap.get("GT", "./.")
            if gt in ("0/0", "0|0", "./.", ".|."):
                continue

            zyg = {"0/1": "Het", "1/0": "Het", "0|1": "Het",
                   "1/1": "Hom", "1|1": "Hom"}.get(gt, gt)

            ad = fmap.get("AD", "")
            vaf = ""
            if ad and ad != ".":
                p = [int(x) for x in ad.split(",")]
                vaf = f"{p[1]/sum(p)*100:.1f}" if sum(p) else ""

            bcf = bcf_lookup.get((chrom, pos, ref, alt), {})
            merged = {**bcf, **info}
            flags = auto_flags(ann, merged, zyg, vaf)

            row = {
                "sample": args.sample,
                "gene": ann.get("GENE", "."),
                "hgvs_c": ann.get("HGVS_C", "."),
                "hgvs_p": ann.get("HGVS_P", "."),
                "chrom": chrom, "pos": pos,
                "zygosity": zyg, "vaf_pct": vaf,
                "impact": ann.get("IMPACT", "."),
                "clinvar_sig": merged.get("CLNSIG", "."),
                "gnomad_af": merged.get("AF", merged.get("AF_sas", ".")),
                "auto_suggested_flags": flags,
                "final_classification": "PENDING",
                "curator": "", "review_date": "", "notes": "",
            }
            for c in acmg_pathogenic + acmg_benign:
                row[c] = ""   # curator fills: Y / N / NA
            rows.append(row)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(k, "")) for k in header) + "\n")

    print(f"ACMG template: {len(rows)} variants → {out}")
    print("Open in Excel/Calc and fill criteria columns + final_classification.")


if __name__ == "__main__":
    main()