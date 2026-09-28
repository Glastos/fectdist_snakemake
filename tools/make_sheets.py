#!/usr/bin/env python3
"""Write config/samples.tsv and config/units.tsv from FASTQ file names.

FASTQ names must follow the provider's pattern:
    {DNA barcode}_{library}_{flowcell}_{lane}_{1|2}.fq.gz
    e.g. GPS263943_MKDN260040585-1A_23M5MYLT4_L4_1.fq.gz
The DNA barcode is looked up in the "CB ech ADN" column of the sample table to
get the pool type ("Type Matrice") and a readable name ("Nom Ext").

Usage:
    tools/make_sheets.py --table ../Table_echantillons_FecDist.xlsm /data/*/*.fq.gz
    find /data -name "*.fq.gz" > files.txt
    tools/make_sheets.py --table ../Table_echantillons_FecDist.xlsm --list files.txt

Nothing is written if any file or sample has a problem; all problems are listed.
"""

import argparse
import os
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

import pandas as pd

FASTQ_NAME = re.compile(
    r"(?P<sample>GPS\d+)_(?P<library>[^_]+)_(?P<flowcell>[^_]+)_(?P<lane>L\d+)"
    r"_(?P<read>[12])\.f(ast)?q\.gz"
)
# "Type Matrice" without accents and case -> pipeline pool type
POOL_TYPES = {"pool male": "drones", "pool ouvriere": "workers"}


def normalise(text):
    text = unicodedata.normalize("NFKD", str(text))
    return "".join(c for c in text if not unicodedata.combining(c)).strip().lower()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("fastq", nargs="*", help="FASTQ files (R1 and R2)")
    parser.add_argument("--list", type=Path, help="text file with one FASTQ path per line")
    parser.add_argument("--table", type=Path, required=True, help="Table_echantillons_FecDist.xlsm")
    parser.add_argument("--samples", type=Path, default=Path("config/samples.tsv"))
    parser.add_argument("--units", type=Path, default=Path("config/units.tsv"))
    parser.add_argument("--drone-ploidy", type=int, default=2)
    parser.add_argument("--worker-ploidy", type=int, default=50)
    parser.add_argument("--force", action="store_true", help="overwrite existing sheets")
    return parser.parse_args()


def main():
    args = parse_args()
    paths = [Path(p) for p in args.fastq]
    if args.list:
        paths += [Path(line.strip()) for line in args.list.read_text().splitlines() if line.strip()]
    if not paths:
        sys.exit("No FASTQ given")

    problems = []

    # (sample, library, flowcell, lane) -> {"1": path, "2": path}
    units = defaultdict(dict)
    for path in paths:
        m = FASTQ_NAME.fullmatch(path.name)
        if not m:
            problems.append(f"name does not match {{barcode}}_{{library}}_{{flowcell}}_{{lane}}_{{1|2}}.fq.gz: {path}")
            continue
        key = (m["sample"], m["library"], m["flowcell"], m["lane"])
        if m["read"] in units[key]:
            problems.append(f"duplicate R{m['read']} for {'_'.join(key)}: {units[key][m['read']]} and {path}")
        units[key][m["read"]] = os.path.abspath(path)
        if not path.exists():
            problems.append(f"file not found: {path}")
    for key, reads in units.items():
        for read in {"1", "2"} - set(reads):
            problems.append(f"missing R{read} for {'_'.join(key)}")
    unit_names = defaultdict(list)
    for sample, library, flowcell, lane in units:
        unit_names[(sample, f"{flowcell}_{lane}")].append(library)
    for (sample, unit), libraries in unit_names.items():
        if len(libraries) > 1:
            problems.append(f"{sample}: several libraries on {unit} {libraries}, unit names would clash")

    table = pd.read_excel(args.table, dtype=str)
    table = table.dropna(subset=["CB ech ADN"]).set_index("CB ech ADN")
    duplicated = set(table.index[table.index.duplicated()])
    ploidy = {"drones": args.drone_ploidy, "workers": args.worker_ploidy}
    samples = {}
    for sample in sorted({key[0] for key in units}):
        if sample not in table.index:
            problems.append(f"{sample}: not in the 'CB ech ADN' column of {args.table}")
            continue
        if sample in duplicated:
            problems.append(f"{sample}: several rows in {args.table}")
            continue
        matrix = table.loc[sample, "Type Matrice"]
        pool_type = POOL_TYPES.get(normalise(matrix))
        if pool_type is None:
            problems.append(f"{sample}: unknown Type Matrice {matrix!r}")
            continue
        samples[sample] = (pool_type, ploidy[pool_type], table.loc[sample, "Nom Ext"])

    for sheet in (args.samples, args.units):
        if sheet.exists() and not args.force:
            problems.append(f"{sheet} exists (use --force to overwrite)")

    if problems:
        sys.exit("Nothing written:\n  " + "\n  ".join(problems))

    args.samples.parent.mkdir(parents=True, exist_ok=True)
    with open(args.samples, "w") as out:
        out.write("sample\ttype\tploidy\tname\n")
        for sample, (pool_type, pool_ploidy, name) in samples.items():
            out.write(f"{sample}\t{pool_type}\t{pool_ploidy}\t{name}\n")

    args.units.parent.mkdir(parents=True, exist_ok=True)
    with open(args.units, "w") as out:
        out.write("sample\tunit\tfq1\tfq2\tlibrary\n")
        for (sample, library, flowcell, lane), reads in sorted(units.items()):
            out.write(f"{sample}\t{flowcell}_{lane}\t{reads['1']}\t{reads['2']}\t{library}\n")

    types = pd.Series([t for t, _, _ in samples.values()]).value_counts().to_dict()
    print(f"{args.samples}: {len(samples)} pools {types}")
    print(f"{args.units}: {len(units)} units")


if __name__ == "__main__":
    main()
