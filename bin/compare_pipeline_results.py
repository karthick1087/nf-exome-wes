#!/usr/bin/env python3
"""Compare bcftools vs GATK pipeline results for a sample."""

import gzip
import sys
from pathlib import Path


TARGETS = [
    ("PAH", "12", 102894920, "c.169-2A>G splice", "Homozygous"),
    ("TBX5", "12", 114399514, "c.361T>C p.Trp121Arg", "Heterozygous"),
]


def vcf_variants(vcf_path: Path, chrom: str, pos: int, window: int = 0):
    if not vcf_path.exists():
        return []
    lo, hi = pos - window, pos + window
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    hits = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if cols[0] != chrom:
                continue
            p = int(cols[1])
            if window and not (lo <= p <= hi):
                if window == 0 and p != pos:
                    continue
                elif window and not (lo <= p <= hi):
                    continue
            if not window and p != pos:
                continue
            fmt = cols[8].split(":")
            sample = cols[9].split(":")
            fmap = dict(zip(fmt, sample))
            gt = fmap.get("GT", "./.")
            ad = fmap.get("AD", "")
            dp = fmap.get("DP", "0")
            vaf = ""
            if ad:
                parts = [int(x) for x in ad.split(",")]
                total = sum(parts)
                if total and len(parts) > 1:
                    vaf = f"{parts[1] / total * 100:.1f}%"
            hits.append({
                "pos": p, "ref": cols[3], "alt": cols[4], "gt": gt,
                "dp": dp, "vaf": vaf, "qual": cols[5], "filter": cols[6],
            })
    return hits


def count_variants(vcf_path: Path) -> int:
    if not vcf_path.exists():
        return -1
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    n = 0
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                n += 1
    return n


def pick_vcf(outdir: Path, sample_id: str, pipeline: str):
    variants = outdir / "variants"
    if pipeline == "gatk":
        for name in (f"{sample_id}.pass.vcf.gz", f"{sample_id}.filtered.vcf.gz"):
            p = variants / name
            if p.exists():
                return p
    else:
        p = variants / f"{sample_id}.filtered.vcf.gz"
        if p.exists():
            return p
    return variants / f"{sample_id}.filtered.vcf.gz"


def main():
    if len(sys.argv) < 3:
        print("Usage: compare_pipeline_results.py <sample_dir> <sample_id>", file=sys.stderr)
        sys.exit(2)
    sample_dir = Path(sys.argv[1])
    sample_id = sys.argv[2]

    bcftools_dir = sample_dir / "results-bcftools"
    gatk_dir = sample_dir / "results-gatk"

    bcftools_vcf = pick_vcf(bcftools_dir, sample_id, "bcftools")
    gatk_vcf = pick_vcf(gatk_dir, sample_id, "gatk")

    out = sample_dir / "pipeline_comparison.tsv"
    lines = [
        "section\tfield\tbcftools\tgatk\tnotes",
        f"summary\tvcf_path\t{bcftools_vcf}\t{gatk_vcf}\t",
        f"summary\tvariant_count\t{count_variants(bcftools_vcf)}\t{count_variants(gatk_vcf)}\t",
        f"summary\tclinical_report\t{bcftools_dir / 'annotation/reports' / f'{sample_id}_clinical_report.tsv'}\t"
        f"{gatk_dir / 'annotation/reports' / f'{sample_id}_clinical_report.tsv'}\t",
    ]

    for gene, chrom, pos, desc, exp_zyg in TARGETS:
        b_hits = vcf_variants(bcftools_vcf, chrom, pos)
        g_hits = vcf_variants(gatk_vcf, chrom, pos)
        b_str = "NOT_FOUND"
        g_str = "NOT_FOUND"
        if b_hits:
            h = b_hits[0]
            b_str = f"{h['ref']}>{h['alt']} GT={h['gt']} DP={h['dp']} VAF={h['vaf']} FILTER={h['filter']}"
        if g_hits:
            h = g_hits[0]
            g_str = f"{h['ref']}>{h['alt']} GT={h['gt']} DP={h['dp']} VAF={h['vaf']} FILTER={h['filter']}"
        match = "MATCH" if b_hits and g_hits and b_hits[0]["alt"] == g_hits[0]["alt"] else (
            "BCF_ONLY" if b_hits and not g_hits else "GATK_ONLY" if g_hits and not b_hits else "DIFF/NONE"
        )
        lines.append(f"target\t{gene}_{chrom}:{pos}\t{b_str}\t{g_str}\t{desc}; expected {exp_zyg}; {match}")

    text = "\n".join(lines) + "\n"
    out.write_text(text)
    print(text)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()