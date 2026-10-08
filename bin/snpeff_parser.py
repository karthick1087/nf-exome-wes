"""Correct SnpEff ANN field parsing (SnpEff 4.x pipe-delimited format)."""

# SnpEff 4/5 ANN is 16 pipe fields. cDNA, CDS, and AA are each one "pos/len" field.
SNPEFF_KEYS = [
    "ALLELE", "ANNOTATION", "IMPACT", "GENE", "GENE_ID",
    "FEATURE_TYPE", "TRANSCRIPT", "BIOTYPE", "RANK",
    "HGVS_C", "HGVS_P", "CDNA", "CDS", "AA", "DISTANCE", "ERRORS",
]


def parse_ann_entries(ann_field):
    if not ann_field:
        return []
    entries = []
    for entry in ann_field.split(","):
        parts = entry.split("|")
        entries.append(dict(zip(SNPEFF_KEYS, parts + [""] * (len(SNPEFF_KEYS) - len(parts)))))
    return entries


def pick_ann_entry(entries, alt=None):
    """Prefer protein_coding transcript with HGVS.c."""
    if not entries:
        return {}
    if alt:
        for e in entries:
            if e.get("ALLELE") == alt and e.get("HGVS_C", "").startswith("c."):
                return e
    for e in entries:
        if e.get("BIOTYPE") == "protein_coding" and e.get("HGVS_C", "").startswith("c."):
            return e
    for e in entries:
        if e.get("BIOTYPE") == "protein_coding":
            return e
    return entries[0]


def humanize_annotation(annotation):
    if not annotation or annotation == ".":
        return "."
    primary = annotation.split("&")[0].replace("_variant", " variant")
    return primary.replace("_", " ")


def format_region(annotation, rank):
    if not rank or rank == "." or "/" not in rank:
        return "."
    idx = rank.split("/")[0]
    ann = (annotation or "").lower()
    if "intron" in ann:
        return f"Intron {idx}"
    if "utr" in ann:
        if "3_prime" in ann:
            return "3' UTR"
        if "5_prime" in ann:
            return "5' UTR"
        return f"UTR {idx}"
    if any(x in ann for x in ("missense", "synonymous", "nonsense", "frameshift", "splice", "exon")):
        return f"Exon {idx}"
    return rank