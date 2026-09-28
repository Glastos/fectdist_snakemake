"""Build the bwa read group of one unit from the first read name of its R1 FASTQ.

PU ({flowcell}.{lane}.{barcode}) is parsed from Illumina (Casava >= 1.8) names:
  @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<barcode>
Any other format stops the job.
"""

import gzip
import re

ILLUMINA = re.compile(
    r"@[^:\s]+:\d+:(?P<flowcell>[^:\s]+):(?P<lane>\d+):\d+:\d+:\d+(:\S+)?"
    r"\s+\d:[YN]:\d+:(?P<barcode>\S+)"
)

fastq = snakemake.input[0]
sample = snakemake.wildcards.sample

with gzip.open(fastq, "rt") as f:
    header = f.readline().strip()

with open(snakemake.log[0], "w") as log:
    log.write(f"First read name: {header}\n")

m = ILLUMINA.fullmatch(header)
if not m:
    raise ValueError(f"{fastq}: first read name is not in Illumina (Casava >= 1.8) format: {header}")

# Literal \t: bwa mem -R expands them
with open(snakemake.output[0], "w") as out:
    out.write(
        f"@RG\\tID:{sample}.{snakemake.wildcards.unit}\\tSM:{sample}\\tPL:ILLUMINA"
        f"\\tLB:{snakemake.params.library}"
        f"\\tPU:{m['flowcell']}.{m['lane']}.{m['barcode']}\n"
    )
