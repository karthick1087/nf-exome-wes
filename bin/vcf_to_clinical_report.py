#!/usr/bin/env python3
"""
Build a clinical-style variant table from annotated VCF files.
Parses SnpEff ANN field, ClinVar INFO, gnomAD AF, genotype, VAF.

Usage:
    python3 vcf_to_clinical_report.py \\
        --vcf sample.snpeff.vcf.gz \\
        --bcftools-vcf sample.fully_annotated.vcf.gz \\
        --clinvar-ref /path/to/clinvar.vcf.gz \\
        --sample SAMPLE001 \\
        --out report.tsv
"""

import argparse
import gzip
from pathlib import Path

from clinical_disorder import (
    fetch_clinvar_nearby,
    find_nearby_clinvar,
    format_disorder,
    get_omim_lookup,
    needs_clinvar_enrichment,
    resolve_disorder_inheritance,
)
from snpeff_parser import parse_ann_entries, pick_ann_entry


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(info_str):
    d = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def parse_gt(fmt, sample_val):
    fmt_keys = fmt.split(":")
    sample_parts = sample_val.split(":")
    fmap = dict(zip(fmt_keys, sample_parts))
    gt = fmap.get("GT", "./.")
    dp = int(fmap.get("DP", 0) or 0)
    ad = fmap.get("AD", "")
    vaf = 0.0
    if ad and ad != ".":
        parts = [int(x) for x in ad.split(",")]
        total = sum(parts)
        alt = parts[1] if len(parts) > 1 else 0
        vaf = (alt / total * 100) if total else 0.0
    zyg = {"0/1": "Heterozygous", "1/0": "Heterozygous",
           "0|1": "Heterozygous", "1|0": "Heterozygous",
           "1/1": "Homozygous",  "1|1": "Homozygous"}.get(gt, gt)
    return zyg, dp, vaf


def load_bcftools_info(bcf_path):
    """Load ClinVar + gnomAD INFO keyed by (chrom, pos, ref, alt)."""
    lookup = {}
    if not bcf_path or not Path(bcf_path).exists():
        return lookup
    with open_text(bcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            key = (cols[0], int(cols[1]), cols[3], cols[4])
            info = parse_info(cols[7])
            lookup[key] = info
    return lookup


def enrich_clinvar(gene, clnsig, clndn, clndisdb, nearby, omim):
    """Apply nearby ClinVar match and resolve disorder/inheritance."""
    if nearby:
        if not clnsig or clnsig == ".":
            clnsig = nearby.get("CLNSIG", clnsig) or clnsig
        if format_disorder(clndn) == ".":
            clndn = nearby.get("CLNDN", clndn) or clndn
        if not clndisdb or clndisdb == ".":
            clndisdb = nearby.get("CLNDISDB", clndisdb) or clndisdb

    disorder, inheritance = resolve_disorder_inheritance(
        gene, clnsig, clndn, clndisdb, omim=omim
    )
    return clnsig, clndn, disorder, inheritance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bcftools-vcf", default="")
    ap.add_argument("--clinvar-ref", default="")
    ap.add_argument("--omim-table", default="")
    ap.add_argument("--annovar", default="")
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    args = ap.parse_args()

    omim = get_omim_lookup(args.omim_table or None)
    bcf_info = load_bcftools_info(args.bcftools_vcf)
    raw_variants = []

    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, pos, vid, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            if alt == ".":
                continue

            info = parse_info(cols[7])
            ann = pick_ann_entry(parse_ann_entries(info.get("ANN", "")), alt=alt)
            zyg, dp, vaf = parse_gt(cols[8], cols[9])

            if zyg in ("Reference", "0/0", "0|0"):
                continue

            bcf = bcf_info.get((chrom, pos, ref, alt), {})
            clnsig = bcf.get("CLNSIG", info.get("CLNSIG", "."))
            clndn = bcf.get("CLNDN", info.get("CLNDN", "."))
            clndisdb = bcf.get("CLNDISDB", info.get("CLNDISDB", "."))
            af = bcf.get("AF", bcf.get("AF_sas", info.get("AF", ".")))
            gene = ann.get("GENE", ".")

            disorder, inheritance = resolve_disorder_inheritance(
                gene, clnsig, clndn, clndisdb, omim=omim
            )

            raw_variants.append({
                "sample": args.sample,
                "gene": gene,
                "transcript": ann.get("TRANSCRIPT", ".") or ".",
                "variant_type": ann.get("ANNOTATION", "."),
                "hgvs_c": ann.get("HGVS_C", "."),
                "hgvs_p": ann.get("HGVS_P", "."),
                "chrom": chrom,
                "pos": pos,
                "rsid": vid if vid != "." else ".",
                "ref": ref,
                "alt": alt,
                "zygosity": zyg,
                "depth": dp,
                "vaf_pct": f"{vaf:.2f}",
                "impact": ann.get("IMPACT", "."),
                "clinvar_sig": clnsig,
                "clinvar_disease": clndn,
                "gnomad_af": af,
                "classification": "REVIEW",
                "_clndisdb": clndisdb,
                "_disorder": disorder,
                "_inheritance": inheritance,
                "_needs_enrichment": needs_clinvar_enrichment(
                    clnsig, clndn, disorder, inheritance
                ),
            })

    enrich_positions = [
        (v["chrom"], v["pos"])
        for v in raw_variants
        if v["_needs_enrichment"] and v["gene"] != "."
    ]
    clinvar_index = fetch_clinvar_nearby(
        args.clinvar_ref, enrich_positions, bcftools=args.bcftools
    )

    rows = []
    filled_disorder = 0
    filled_inheritance = 0
    for v in raw_variants:
        disorder = v["_disorder"]
        inheritance = v["_inheritance"]
        nearby = {}
        if v["_needs_enrichment"]:
            nearby = find_nearby_clinvar(
                v["chrom"], v["pos"], v["gene"], clinvar_index
            )
        clnsig, clndn, disorder, inheritance = enrich_clinvar(
            v["gene"], v["clinvar_sig"], v["clinvar_disease"],
            v["_clndisdb"], nearby, omim,
        )
        if disorder != ".":
            filled_disorder += 1
        if inheritance != ".":
            filled_inheritance += 1
        v["clinvar_sig"] = clnsig
        v["clinvar_disease"] = clndn
        v["disorder"] = disorder
        v["inheritance"] = inheritance
        del v["_clndisdb"]
        del v["_disorder"]
        del v["_inheritance"]
        del v["_needs_enrichment"]
        rows.append(v)

    header = [
        "sample", "gene", "transcript", "variant_type", "hgvs_c", "hgvs_p",
        "chrom", "pos", "rsid", "ref", "alt", "zygosity", "depth", "vaf_pct",
        "impact", "clinvar_sig", "clinvar_disease", "disorder", "inheritance",
        "gnomad_af", "classification",
    ]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in header) + "\n")

    print(
        f"Clinical report: {len(rows)} variants → {out}\n"
        f"  disorder filled: {filled_disorder}/{len(rows)}\n"
        f"  inheritance filled: {filled_inheritance}/{len(rows)}"
    )


if __name__ == "__main__":
    main()