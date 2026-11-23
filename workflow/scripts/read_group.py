"""Build the bwa read group of one unit from the first read name of its R1 FASTQ.

PU ({flowcell}.{lane}.{barcode}) is parsed from Illumina-style names:
  @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<barcode>
as written by Illumina and by AVITI
Any other format stops the job. PL is the platform of the unit in the sample sheet.
"""

import gzip
import re

READ_NAME = re.compile(
    r"@[^:\s]+:[^:\s]+:(?P<flowcell>[^:\s]+):(?P<lane>\d+):\d+:\d+:\d+(:\S+)?"
    r"\s+\d:[YN]:\d+:(?P<barcode>\S+)"
)

fastq = snakemake.input[0]
sample = snakemake.wildcards.sample

with gzip.open(fastq, "rt") as f:
    header = f.readline().strip()

with open(snakemake.log[0], "w") as log:
    log.write(f"First read name: {header}\n")

m = READ_NAME.fullmatch(header)
if not m:
    raise ValueError(
        f"{fastq}: first read name is not in Illumina-style format (@<instrument>:<run>:"
        f"<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<barcode>): {header}"
    )

# Literal \t: bwa mem -R expands them
with open(snakemake.output[0], "w") as out:
    out.write(
        f"@RG\\tID:{sample}.{snakemake.wildcards.unit}\\tSM:{sample}"
        f"\\tPL:{snakemake.params.platform}"
        f"\\tLB:{snakemake.params.library}"
        f"\\tPU:{m['flowcell']}.{m['lane']}.{m['barcode']}\n"
    )
