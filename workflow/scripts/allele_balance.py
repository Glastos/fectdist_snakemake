"""Allele balance of a drone pool at the SNPs called 0/1 and 1/1 in its gVCF.

Keeps biallelic SNPs (first ALT allele, <NON_REF> ignored), with a depth
(ref + alt reads) >= min_depth and <= max_depth_ratio x the median depth of all
these SNPs, then writes a random sample of at most max_sites per class.
"""

import gzip
import random
import statistics
import subprocess

params = snakemake.params
query = subprocess.Popen(
    [
        "bcftools", "query",
        "-i", 'GT="het" || GT="AA"',
        "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%AD]\n",
        snakemake.input.gvcf,
    ],
    stdout=subprocess.PIPE,
    text=True,
)

sites = {"het": [], "hom": []}
for line in query.stdout:
    chrom, pos, ref, alt, gt, ad = line.rstrip("\n").split("\t")
    alt = alt.split(",")[0]
    if len(ref) != 1 or alt not in ("A", "C", "G", "T"):
        continue
    gt = gt.replace("|", "/")
    if gt == "0/1":
        genotype = "het"
    elif gt == "1/1":
        genotype = "hom"
    else:
        continue
    counts = ad.split(",")
    ref_ad, alt_ad = int(counts[0]), int(counts[1])
    if ref_ad + alt_ad > 0:
        sites[genotype].append((chrom, int(pos), ref_ad, alt_ad))
if query.wait() != 0:
    raise RuntimeError(f"bcftools query failed on {snakemake.input.gvcf}")

depths = [r + a for rows in sites.values() for _, _, r, a in rows]
median_depth = float(statistics.median(depths)) if depths else 0.0
max_depth = params.max_depth_ratio * median_depth
rng = random.Random(params.seed)

summary = {"sample": snakemake.wildcards.sample, "median_depth": median_depth,
           "min_depth": params.min_depth, "max_depth": max_depth}
with gzip.open(snakemake.output.sites, "wt") as out:
    out.write("genotype\tchrom\tpos\tref_ad\talt_ad\n")
    for genotype, rows in sites.items():
        kept = [s for s in rows if params.min_depth <= s[2] + s[3] <= max_depth]
        sampled = sorted(rng.sample(kept, min(len(kept), params.max_sites)))
        summary |= {f"n_{genotype}": len(rows), f"n_{genotype}_kept": len(kept),
                    f"n_{genotype}_sampled": len(sampled)}
        for chrom, pos, ref_ad, alt_ad in sampled:
            out.write(f"{genotype}\t{chrom}\t{pos}\t{ref_ad}\t{alt_ad}\n")

with open(snakemake.output.summary, "w") as out:
    out.write("\t".join(summary) + "\n")
    out.write("\t".join(str(v) for v in summary.values()) + "\n")
