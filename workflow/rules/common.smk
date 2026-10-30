import re
from pathlib import Path

import pandas as pd
from snakemake.exceptions import WorkflowError

# Output files each module can build, by target name (see `targets` in config)
TARGETS = {}

SAFE_NAME = re.compile(r"^[A-Za-z0-9._-]+$")


# ------------------------------------------------------------------------------
# Sample sheet: one row per read group (unit), the pool columns repeated on
# every row of the pool
# ------------------------------------------------------------------------------

SHEET = config["samples"]
POOL_COLUMNS = ["type", "ploidy"]

units = pd.read_csv(SHEET, sep="\t", dtype=str, comment="#").fillna("")
_missing = {"sample", "unit", "fq1", "fq2", *POOL_COLUMNS} - set(units.columns)
if _missing:
    raise WorkflowError(f"{SHEET}: missing column(s) {sorted(_missing)}")

# One library per pool unless stated otherwise
if "library" not in units:
    units["library"] = ""
units["library"] = units["library"].where(units["library"] != "", units["sample"])

bad_names = [n for n in [*units["sample"], *units["unit"]] if not SAFE_NAME.match(n)]
if bad_names:
    raise WorkflowError(f"Sample/unit names may only use A-Z a-z 0-9 . _ -: {bad_names}")
if "all" in set(units["sample"]):
    raise WorkflowError(f"{SHEET}: 'all' is reserved (results/all/)")
if units.duplicated(["sample", "unit"]).any():
    raise WorkflowError(f"{SHEET}: duplicated (sample, unit) pairs")

# No effect on the pipeline yet, kept for downstream steps
POOL_TYPES = {"drones", "workers"}
if not set(units["type"]) <= POOL_TYPES:
    raise WorkflowError(
        f"{SHEET}: unknown type(s) {sorted(set(units['type']) - POOL_TYPES)}, "
        f"expected one of {sorted(POOL_TYPES)}"
    )
if not units["ploidy"].str.fullmatch(r"[1-9]\d*").all():
    raise WorkflowError(f"{SHEET}: ploidy must be a positive integer")

_values = units.groupby("sample")[POOL_COLUMNS].nunique()
_differ = _values.index[(_values > 1).any(axis=1)].tolist()
if _differ:
    raise WorkflowError(
        f"{SHEET}: {' or '.join(POOL_COLUMNS)} differ between the rows of {_differ}"
    )

units["ploidy"] = units["ploidy"].astype(int)
samples = (
    units.drop_duplicates("sample")[["sample", *POOL_COLUMNS]]
    .set_index("sample", drop=False)
)
units = units.set_index(["sample", "unit"], drop=False)

SAMPLES = samples.index.tolist()
ROUNDS = int(config["bqsr"]["rounds"])


def get_unit(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit)]


def sample_units(sample):
    return units.loc[units["sample"] == sample, "unit"].tolist()


# ------------------------------------------------------------------------------
# Reference (indexes must be provided, the pipeline does not build them)
# ------------------------------------------------------------------------------

REF = config["reference"]["fasta"]
REF_FAI = f"{REF}.fai"
REF_STEM = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", "", REF)
REF_DICT = f"{REF_STEM}.dict"
REF_BWA = [f"{REF}{ext}" for ext in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]

_missing = [f for f in [REF, REF_FAI, REF_DICT, *REF_BWA] if not Path(f).exists()]
if _missing:
    raise WorkflowError(
        "Reference files missing:\n  "
        + "\n  ".join(_missing)
        + f"\nBuild them with:\n  bwa-mem2 index {REF}\n  samtools faidx {REF}"
        f"\n  gatk CreateSequenceDictionary -R {REF}"
    )


# ------------------------------------------------------------------------------
# Regions for parallel calling
# ------------------------------------------------------------------------------


def make_regions(fai, min_contig_size, include, exclude):
    """Regions for parallel calling, each a list of (contig, start, end) 0-based.

    Only contigs in include (all if empty) and not in exclude are kept. Contigs
    >= min_contig_size get one region each; all smaller contigs (unplaced
    scaffolds, mitochondrion) share one last region.
    """
    regions, small = [], []
    with open(fai) as f:
        for line in f:
            contig, length = line.split("\t")[:2]
            if (include and contig not in include) or contig in exclude:
                continue
            interval = (contig, 0, int(length))
            if int(length) >= min_contig_size:
                regions.append([interval])
            else:
                small.append(interval)
    if small:
        regions.append(small)
    return regions


_include = set(config["regions"]["include"])
_exclude = set(config["regions"]["exclude"])
_contigs = {line.split("\t")[0] for line in open(REF_FAI)}
for _key, _names in [("include", _include), ("exclude", _exclude)]:
    if not _names <= _contigs:
        raise WorkflowError(
            f"regions.{_key}: contigs not in {REF_FAI}: {sorted(_names - _contigs)}"
        )

REGIONS = {
    f"{i:04d}": region
    for i, region in enumerate(
        make_regions(
            REF_FAI,
            int(config["regions"]["min_contig_size"]),
            _include,
            _exclude,
        ),
        start=1,
    )
}


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------


def pool_ploidy(sample):
    """HaplotypeCaller ploidy of the bootstrap calls and final gVCF (samples.tsv)."""
    return int(samples.loc[sample, "ploidy"])


def java_opts(wildcards, resources):
    """Leave 20% of the job memory to the JVM overhead."""
    return f"-Xmx{int(resources.mem_mb * 0.8)}m"


def bam_after_round(sample, r):
    """BAM after BQSR round r (round 0 = duplicate-marked BAM)."""
    if r == 0:
        prefix = f"results/{sample}/mapping/dedup"
    else:
        prefix = f"results/{sample}/bqsr/round{r}.recal"
    return {"bam": f"{prefix}.bam", "bai": f"{prefix}.bai"}


def round_input(wildcards):
    """BAM entering the BQSR round of the wildcards."""
    return bam_after_round(wildcards.sample, int(wildcards.round) - 1)


wildcard_constraints:
    sample="|".join(map(re.escape, SAMPLES)),
    unit="|".join(map(re.escape, units["unit"].unique())),
    region=r"\d+",
    round=r"\d+",
