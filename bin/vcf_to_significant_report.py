#!/usr/bin/env python3
"""
Build a significant-findings table for pathogenic / likely pathogenic variants.

Output columns:
  Gene, Transcript, Variant Type, hgvs_c, hgvs_p, genomic_coord, rsid, region,
  Zygosity, Coverage, VAF, Disorder, Inheritance, Classification, Attributes

Usage:
    python3 vcf_to_significant_report.py \\
        --vcf sample.snpeff.vcf.gz \\
        --bcftools-vcf sample.fully_annotated.vcf.gz \\
        --clinvar-ref clinvar.vcf.gz \\
        --sample SAMPLE001 \\
        --out significant_findings.tsv
"""

import argparse
import gzip
import re
import subprocess
from pathlib import Path

from clinical_disorder import (
    fetch_clinvar_nearby,
    find_nearby_clinvar,
    format_disorder,
    get_clinvar_disease_lookup,
    get_omim_lookup,
    needs_clinvar_enrichment,
    resolve_disorder_inheritance,
)
from snpeff_parser import format_region, humanize_annotation, parse_ann_entries, pick_ann_entry


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
    zyg = {
        "0/1": "Heterozygous", "1/0": "Heterozygous",
        "0|1": "Heterozygous", "1|0": "Heterozygous",
        "1/1": "Homozygous", "1|1": "Homozygous",
    }.get(gt, gt)
    return zyg, dp, vaf


def normalize_classification(clnsig):
    if not clnsig or clnsig == ".":
        return ".", 99
    best = (".", 99)
    for part in re.split(r"[|,/]", clnsig):
        low = part.lower().replace("_", " ")
        if "benign" in low:
            continue
        if "likely" in low and "pathogenic" in low:
            cand = ("Likely pathogenic", 1)
        elif "pathogenic" in low:
            cand = ("Pathogenic", 0)
        elif "conflicting" in low:
            cand = ("Conflicting interpretations of pathogenicity", 2)
        elif "risk factor" in low or "risk_factor" in part.lower():
            cand = ("Risk factor", 3)
        else:
            continue
        if cand[1] < best[1]:
            best = cand
    return best


def is_significant(clnsig):
    label, rank = normalize_classification(clnsig)
    return label != ".", label, rank


def suggest_attributes(ann, clnsig, zyg, vaf_pct, inheritance, gnomad_af):
    """Suggest ACMG-style evidence tags (PDF format)."""
    attrs = []
    annotation = ann.get("ANNOTATION", "")
    impact = ann.get("IMPACT", "")
    clow = (clnsig or "").lower()

    # PM2 — absent or extremely rare in population
    try:
        if not gnomad_af or gnomad_af == "." or float(gnomad_af) < 0.0001:
            attrs.append("PM2")
    except ValueError:
        attrs.append("PM2")

    # PVS1 — null variant (LOF)
    if any(x in annotation for x in (
        "splice_acceptor", "splice_donor", "frameshift", "stop_gained", "start_lost"
    )):
        attrs.append("PVS1")

    # PM1 — missense in functionally important domain
    if impact in ("HIGH", "MODERATE") and "missense" in annotation:
        attrs.append("PM1")

    # PM3 — recessive disorder (compound het or homozygous)
    if inheritance and "recessive" in inheritance.lower():
        attrs.append("PM3")

    # PP2 — missense in gene where missense is common mechanism
    if "missense" in annotation:
        attrs.append("PP2")

    # PP3 — computational evidence for missense
    if "missense" in annotation:
        attrs.append("PP3")

    # PP5 — reputable source / ClinVar pathogenic
    if "pathogenic" in clow:
        attrs.append("PP5")

    # PP1_Moderate — segregation evidence (recessive / unclear phase; manual confirm)
    if "pathogenic" in clow and inheritance and inheritance != ".":
        if "recessive" in inheritance.lower() or zyg == "Homozygous":
            attrs.append("PP1_Moderate")

    # De-duplicate preserving order
    seen = set()
    ordered = []
    for a in attrs:
        if a not in seen:
            seen.add(a)
            ordered.append(a)
    return ", ".join(ordered) if ordered else "."


def load_bcftools_info(bcf_path):
    lookup = {}
    if not bcf_path or not Path(bcf_path).exists():
        return lookup
    with open_text(bcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            key = (cols[0], int(cols[1]), cols[3], cols[4])
            lookup[key] = parse_info(cols[7])
    return lookup


def load_rsid_supplement(path=None):
    """Optional position-based rsID overrides (newer dbSNP entries)."""
    path = Path(path or Path(__file__).resolve().parent / "resources" / "rsid_supplement_hg38.tsv")
    rsmap = {}
    if not path.exists():
        return rsmap
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            try:
                chrom, pos, ref, alt, rsid = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            except ValueError:
                continue
            if rsid and rsid != ".":
                rsmap[(chrom.replace("chr", ""), pos, ref, alt)] = rsid
    return rsmap


def fetch_rsid_map(dbsnp_vcfs, positions, bcftools="bcftools"):
    """Map (chrom, pos, ref, alt) -> rsID from one or more dbSNP VCFs."""
    if not dbsnp_vcfs or not positions:
        return {}
    if isinstance(dbsnp_vcfs, (str, Path)):
        dbsnp_vcfs = [dbsnp_vcfs]

    by_chrom = {}
    for chrom, pos, ref, alt in positions:
        by_chrom.setdefault(chrom, set()).add(int(pos))

    rsmap = {}
    fmt = "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n"
    for dbsnp_vcf in dbsnp_vcfs:
        if not dbsnp_vcf or not Path(dbsnp_vcf).exists():
            continue
        for chrom, pos_set in by_chrom.items():
            min_pos = max(1, min(pos_set) - 2)
            max_pos = max(pos_set) + 2
            for region in (f"{chrom}:{min_pos}-{max_pos}", f"chr{chrom}:{min_pos}-{max_pos}"):
                try:
                    proc = subprocess.run(
                        [bcftools, "query", "-f", fmt, str(dbsnp_vcf), "-r", region],
                        check=True, capture_output=True, text=True,
                    )
                except (subprocess.CalledProcessError, FileNotFoundError):
                    continue
                for line in proc.stdout.splitlines():
                    cols = line.split("\t")
                    if len(cols) < 5:
                        continue
                    c, p, rid, r, a = cols[0], cols[1], cols[2], cols[3], cols[4]
                    c = c.replace("chr", "")
                    if rid and rid != ".":
                        rsmap[(c, int(p), r, a)] = rid
                        if "," in a:
                            for alt_allele in a.split(","):
                                rsmap[(c, int(p), r, alt_allele)] = rid
    return rsmap


def enrich_clinvar(gene, clnsig, clndn, clndisdb, nearby, omim):
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
    ap.add_argument("--dbsnp-vcf", default="")
    ap.add_argument("--omim-table", default="")
    ap.add_argument("--sample", default="SAMPLE")
    ap.add_argument("--out", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    args = ap.parse_args()

    omim = get_omim_lookup(args.omim_table or None)
    get_clinvar_disease_lookup()
    bcf_info = load_bcftools_info(args.bcftools_vcf)

    candidates = []
    with open_text(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            chrom, pos, vid, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]
            if alt == ".":
                continue

            info = parse_info(cols[7])
            zyg, dp, vaf = parse_gt(cols[8], cols[9])
            if zyg in ("Reference", "0/0", "0|0", "./.", ".|."):
                continue

            bcf = bcf_info.get((chrom, pos, ref, alt), {})
            clnsig = bcf.get("CLNSIG", info.get("CLNSIG", "."))
            clndn = bcf.get("CLNDN", info.get("CLNDN", "."))
            clndisdb = bcf.get("CLNDISDB", info.get("CLNDISDB", "."))

            entries = parse_ann_entries(info.get("ANN", ""))
            ann = pick_ann_entry(entries, alt)
            gene = ann.get("GENE", ".")
            impact = ann.get("IMPACT", ".")

            sig_ok, classification, sig_rank = is_significant(clnsig)
            if not sig_ok:
                clow = (clnsig or "").lower()
                if "benign" in clow or gene == ".":
                    continue
                if impact not in ("HIGH", "MODERATE"):
                    continue

            disorder, inheritance = resolve_disorder_inheritance(
                gene, clnsig, clndn, clndisdb, omim=omim
            )
            needs_nearby = needs_clinvar_enrichment(clnsig, clndn, disorder, inheritance)

            candidates.append({
                "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                "vid": vid, "zyg": zyg, "dp": dp, "vaf": vaf,
                "clnsig": clnsig, "clndn": clndn, "clndisdb": clndisdb,
                "classification": classification, "sig_rank": sig_rank,
                "ann": ann, "gene": gene,
                "disorder": disorder, "inheritance": inheritance,
                "needs_nearby": needs_nearby,
                "gnomad_af": bcf.get("AF", bcf.get("AF_sas", info.get("AF", "."))),
            })

    nearby_positions = [
        (c["chrom"], c["pos"]) for c in candidates if c["needs_nearby"] and c["gene"] != "."
    ]
    clinvar_index = fetch_clinvar_nearby(
        args.clinvar_ref, nearby_positions, bcftools=args.bcftools
    )

    rs_positions = [(c["chrom"], c["pos"], c["ref"], c["alt"]) for c in candidates]
    dbsnp_sources = [args.dbsnp_vcf] if args.dbsnp_vcf else []
    rsmap = fetch_rsid_map(dbsnp_sources, rs_positions, bcftools=args.bcftools)
    rsmap.update(load_rsid_supplement())

    rows = []
    for c in candidates:
        nearby = {}
        if c["needs_nearby"]:
            nearby = find_nearby_clinvar(c["chrom"], c["pos"], c["gene"], clinvar_index)
        clnsig, clndn, disorder, inheritance = enrich_clinvar(
            c["gene"], c["clnsig"], c["clndn"], c["clndisdb"], nearby, omim
        )
        sig_ok, classification, sig_rank = is_significant(clnsig)
        if not sig_ok:
            continue

        ann = c["ann"]
        rsid = c["vid"] if c["vid"] not in (".", "") else "."
        rsid = rsmap.get((c["chrom"], c["pos"], c["ref"], c["alt"]), rsid)

        hgvs_c = ann.get("HGVS_C", ".")
        hgvs_p = ann.get("HGVS_P", ".")
        if hgvs_p in (".", ""):
            hgvs_p = "."

        rows.append({
            "sample": args.sample,
            "gene": c["gene"],
            "transcript": ann.get("TRANSCRIPT", "."),
            "variant_type": humanize_annotation(ann.get("ANNOTATION", ".")),
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "genomic_coord": f"{c['chrom']}:{c['pos']}",
            "rsid": rsid,
            "region": format_region(ann.get("ANNOTATION", ""), ann.get("RANK", ".")),
            "zygosity": c["zyg"],
            "coverage": c["dp"],
            "vaf_pct": f"{c['vaf']:.2f}%",
            "disorder": disorder,
            "inheritance": inheritance,
            "classification": classification,
            "attributes": suggest_attributes(
                ann, clnsig, c["zyg"], c["vaf"], inheritance, c["gnomad_af"]
            ),
            "clinvar_sig": clnsig,
            "impact": ann.get("IMPACT", "."),
            "_sig_rank": sig_rank,
            "_gene_sort": c["gene"],
        })

    rows.sort(key=lambda r: (r["_sig_rank"], r["_gene_sort"], r["genomic_coord"]))

    header = [
        "sample", "gene", "transcript", "variant_type", "hgvs_c", "hgvs_p",
        "genomic_coord", "rsid", "region", "zygosity", "coverage", "vaf_pct",
        "disorder", "inheritance", "classification", "attributes",
        "clinvar_sig", "impact",
    ]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in header) + "\n")

    print(f"Significant findings: {len(rows)} variants → {out}")


if __name__ == "__main__":
    main()