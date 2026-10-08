"""Disorder and inheritance helpers for clinical variant reports."""

import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

SKIP_DISEASES = {
    "not_specified", "not_provided", "inborn_genetic_diseases",
    "see_cases", "cardiovascular_phenotype", "inborn genetic diseases",
}

GENE_DEFAULTS = {
    "PAH": ("Phenylketonuria", "Autosomal recessive"),
    "TBX5": ("Holt-Oram syndrome", "Autosomal dominant"),
}

DEFAULT_OMIM_TABLE = Path(__file__).resolve().parent / "resources" / "use_omim_table.txt"
DEFAULT_CLINVAR_INH = Path(__file__).resolve().parent / "resources" / "clinvar_disease_inheritance.tsv"


def normalize_key(text):
    if not text:
        return ""
    text = text.lower().replace("_", " ")
    text = re.sub(r"[^a-z0-9 ]+", "", text)
    return re.sub(r"\s+", " ", text).strip()


def format_disorder(clndn):
    """Pick the best human-readable disorder from CLNDN."""
    if not clndn or clndn == ".":
        return "."
    parts = [p.strip() for p in clndn.split("|") if p.strip()]
    for part in parts:
        low = part.lower()
        if low in SKIP_DISEASES or part.startswith("MedGen:"):
            continue
        return part.replace("_", " ")
    return "."


def parse_omim_ids(clndisdb):
    if not clndisdb or clndisdb == ".":
        return []
    return re.findall(r"OMIM:(\d+)", clndisdb)


def parse_geneinfo(geneinfo):
    if not geneinfo or geneinfo == ".":
        return ""
    return geneinfo.split(":")[0]


def infer_inheritance_from_name(disorder):
    """Infer inheritance mode embedded in a disorder name."""
    if not disorder or disorder == ".":
        return "."
    n = disorder.lower()
    patterns = [
        (r"autosomal recessive", "Autosomal recessive"),
        (r"autosomal dominant", "Autosomal dominant"),
        (r"x-linked recessive", "X-linked recessive"),
        (r"x-linked dominant", "X-linked dominant"),
        (r"x-linked", "X-linked"),
        (r"mitochondrial", "Mitochondrial"),
        (r"y-linked", "Y-linked"),
        (r"\bdfnb\d", "Autosomal recessive"),
        (r"\bdfna\d", "Autosomal dominant"),
        (r"\bar\b", "Autosomal recessive"),
        (r"\bad\b", "Autosomal dominant"),
    ]
    for pat, inh in patterns:
        if re.search(pat, n):
            return inh
    return "."


def normalize_inheritance(raw):
    if not raw or raw in (".", "None", "NA"):
        return "."
    parts = [p.strip() for p in raw.split(";") if p.strip() and p.strip() not in ("None", "NA")]
    if not parts:
        return "."
    order = {
        "autosomal dominant": 0,
        "autosomal recessive": 1,
        "x-linked dominant": 2,
        "x-linked recessive": 3,
        "x-linked": 4,
        "mitochondrial": 5,
        "digenic recessive": 6,
        "digenic dominant": 7,
        "multifactorial": 8,
        "somatic mutation": 9,
        "isolated cases": 10,
    }
    parts.sort(key=lambda x: order.get(x.lower(), 99))
    return parts[0]


class OmimLookup:
    def __init__(self, table_path=None):
        self.table_path = Path(table_path or DEFAULT_OMIM_TABLE)
        self.by_phenotype_mim = {}
        self.by_phenotype_name = defaultdict(list)
        self.by_gene = defaultdict(list)
        self._loaded = False

    def load(self):
        if self._loaded or not self.table_path.exists():
            self._loaded = True
            return
        with open(self.table_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                phenotype = (row.get("phenotype") or "").strip()
                inheritance = normalize_inheritance(row.get("phenotypeInheritance", ""))
                mim = (row.get("phenotypeMimNumber") or "").strip()
                genes = (row.get("hgnc_genes") or row.get("genes") or "").strip()
                if not phenotype:
                    continue
                entry = (phenotype, inheritance)
                if mim and mim != "NA":
                    self.by_phenotype_mim[mim] = entry
                key = normalize_key(phenotype)
                if key and inheritance != ".":
                    self.by_phenotype_name[key].append(entry)
                for gene in re.split(r"[,|]", genes):
                    gene = gene.strip()
                    if gene and gene != "NA":
                        self.by_gene[gene.upper()].append(entry)
        self._loaded = True

    def inheritance_for_phenotype_mim(self, mim):
        self.load()
        entry = self.by_phenotype_mim.get(str(mim))
        return entry[1] if entry else "."

    def inheritance_for_disorder(self, disorder):
        self.load()
        if not disorder or disorder == ".":
            return "."
        key = normalize_key(disorder)
        entries = self.by_phenotype_name.get(key, [])
        if entries:
            return entries[0][1]
        for phen_key, vals in self.by_phenotype_name.items():
            if key in phen_key or phen_key in key:
                return vals[0][1]
        return infer_inheritance_from_name(disorder)

    def gene_phenotypes(self, gene):
        self.load()
        return self.by_gene.get((gene or "").upper(), [])

    def best_gene_match(self, gene, clndn=""):
        """Best OMIM phenotype for a gene, optionally guided by CLNDN tokens."""
        phenotypes = self.gene_phenotypes(gene)
        if not phenotypes:
            return ".", "."

        tokens = [normalize_key(t) for t in (clndn or "").split("|") if t]
        tokens = [t for t in tokens if t and t not in SKIP_DISEASES]

        if tokens:
            for disorder, inheritance in phenotypes:
                dkey = normalize_key(disorder)
                for tok in tokens:
                    if tok in dkey or dkey in tok:
                        return disorder, inheritance

        if len(phenotypes) == 1:
            return phenotypes[0]

        scored = []
        clndn_text = normalize_key(clndn)
        for disorder, inheritance in phenotypes:
            dkey = normalize_key(disorder)
            score = 0
            if inheritance != ".":
                score += 1
            if clndn_text and dkey in clndn_text:
                score += 5
            scored.append((score, disorder, inheritance))
        scored.sort(reverse=True)
        return scored[0][1], scored[0][2]


class ClinvarDiseaseLookup:
    def __init__(self, table_path=None):
        self.table_path = Path(table_path or DEFAULT_CLINVAR_INH)
        self.by_name = {}
        self._loaded = False

    def load(self):
        if self._loaded or not self.table_path.exists():
            self._loaded = True
            return
        with open(self.table_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                disorder = (row.get("disorder") or "").strip()
                inheritance = normalize_inheritance(row.get("inheritance", ""))
                if disorder and inheritance != ".":
                    self.by_name[normalize_key(disorder)] = inheritance
        self._loaded = True

    def inheritance_for_disorder(self, disorder):
        self.load()
        if not disorder or disorder == ".":
            return "."
        key = normalize_key(disorder)
        if key in self.by_name:
            return self.by_name[key]
        gene = related_disorder_gene(disorder)
        if gene:
            return "."
        return "."


def related_disorder_gene(disorder):
    m = re.match(r"^([A-Z][A-Z0-9]+)-related disorder$", disorder or "")
    return m.group(1) if m else ""


_OMIM = OmimLookup()
_CLINVAR_DIS = ClinvarDiseaseLookup()


def get_omim_lookup(table_path=None):
    global _OMIM
    if table_path:
        path = Path(table_path)
        if path != _OMIM.table_path:
            _OMIM = OmimLookup(path)
    _OMIM.load()
    return _OMIM


def get_clinvar_disease_lookup(table_path=None):
    global _CLINVAR_DIS
    if table_path:
        path = Path(table_path)
        if path != _CLINVAR_DIS.table_path:
            _CLINVAR_DIS = ClinvarDiseaseLookup(path)
    _CLINVAR_DIS.load()
    return _CLINVAR_DIS


def clinvar_gene_match(info, gene):
    if not gene or gene == ".":
        return False
    geneinfo = parse_geneinfo(info.get("GENEINFO", ""))
    if geneinfo and geneinfo.upper() == gene.upper():
        return True
    clndn = info.get("CLNDN", "")
    low = clndn.lower()
    if gene.upper() == "TBX5" and "holt-oram" in low:
        return True
    if gene.upper() == "PAH" and "phenylketonuria" in low:
        return True
    return False


def pathogenic_rank(clnsig):
    sig = (clnsig or "").lower()
    if "pathogenic" in sig and "likely" not in sig:
        return 0
    if "likely_pathogenic" in sig:
        return 1
    if "uncertain" in sig:
        return 2
    if "conflicting" in sig:
        return 3
    return 4


def is_actionable_clinsig(clnsig):
    sig = (clnsig or "").lower()
    return any(x in sig for x in (
        "pathogenic", "likely_pathogenic", "uncertain", "conflicting", "risk_factor"
    ))


def resolve_disorder_inheritance(gene, clnsig, clndn, clndisdb, omim=None, clinvar_dis=None):
    """Derive disorder and inheritance from ClinVar + OMIM resources."""
    omim = omim or get_omim_lookup()
    clinvar_dis = clinvar_dis or get_clinvar_disease_lookup()
    disorder = format_disorder(clndn)
    inheritance = "."

    for omim_id in parse_omim_ids(clndisdb):
        inh = omim.inheritance_for_phenotype_mim(omim_id)
        if inh != ".":
            inheritance = inh
            if disorder == ".":
                entry = omim.by_phenotype_mim.get(omim_id)
                if entry:
                    disorder = entry[0]
            break

    if disorder == "." and gene in GENE_DEFAULTS:
        disorder, default_inh = GENE_DEFAULTS[gene]
        if inheritance == ".":
            inheritance = default_inh

    if disorder == "." and gene != "." and is_actionable_clinsig(clnsig):
        g_dis, g_inh = omim.best_gene_match(gene, clndn)
        if g_dis != ".":
            disorder = g_dis
            if inheritance == "." and g_inh != ".":
                inheritance = g_inh

    if inheritance == "." and disorder != ".":
        inheritance = omim.inheritance_for_disorder(disorder)

    if inheritance == "." and disorder != ".":
        inheritance = clinvar_dis.inheritance_for_disorder(disorder)

    if inheritance == "." and clndn and clndn != ".":
        for part in clndn.split("|"):
            part = part.strip().replace("_", " ")
            if not part or part.lower() in SKIP_DISEASES:
                continue
            inheritance = omim.inheritance_for_disorder(part)
            if inheritance == ".":
                inheritance = clinvar_dis.inheritance_for_disorder(part)
            if inheritance != ".":
                if disorder == ".":
                    disorder = part
                break

    if inheritance == "." and disorder != ".":
        rel_gene = related_disorder_gene(disorder)
        if rel_gene:
            _, g_inh = omim.best_gene_match(rel_gene, clndn)
            if g_inh != ".":
                inheritance = g_inh

    return disorder, inheritance


def fetch_clinvar_nearby(clinvar_vcf, positions, window=5, bcftools="bcftools"):
    if not clinvar_vcf or not positions:
        return {}

    by_chrom = defaultdict(set)
    for chrom, pos in positions:
        by_chrom[chrom].add(int(pos))

    index = defaultdict(list)
    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO\\n"

    for chrom, pos_set in by_chrom.items():
        min_pos = max(1, min(pos_set) - window)
        max_pos = max(pos_set) + window
        try:
            proc = subprocess.run(
                [bcftools, "query", "-f", fmt, clinvar_vcf, "-r", f"{chrom}:{min_pos}-{max_pos}"],
                check=True,
                capture_output=True,
                text=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue

        for line in proc.stdout.splitlines():
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            cpos = int(cols[1])
            info = {}
            for part in cols[4].split(";"):
                if "=" in part:
                    k, v = part.split("=", 1)
                    info[k] = v
            index[(chrom, cpos)].append(info)

    return index


def find_nearby_clinvar(chrom, pos, gene, clinvar_index, window=5):
    candidates = []
    for offset in range(-window, window + 1):
        cpos = pos + offset
        for info in clinvar_index.get((chrom, cpos), []):
            if not clinvar_gene_match(info, gene):
                continue
            clndn = info.get("CLNDN", "")
            disorder = format_disorder(clndn)
            if disorder == ".":
                continue
            candidates.append((abs(offset), info, disorder))

    if not candidates:
        return {}

    gene_default = GENE_DEFAULTS.get(gene)
    if gene_default:
        preferred = [c for c in candidates if c[2].lower() == gene_default[0].lower()]
        if preferred:
            candidates = preferred

    candidates.sort(key=lambda x: (x[0], pathogenic_rank(x[1].get("CLNSIG", ""))))
    return candidates[0][1]


def needs_clinvar_enrichment(clnsig, clndn, disorder, inheritance):
    if format_disorder(clndn) == "." and (not clnsig or clnsig == "."):
        return True
    if disorder == "." and is_actionable_clinsig(clnsig):
        return True
    if disorder != "." and inheritance == ".":
        return True
    if format_disorder(clndn) == "." and is_actionable_clinsig(clnsig):
        return True
    return False