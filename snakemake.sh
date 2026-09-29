#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --partition=unlimitq
#SBATCH --output=snakemake.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guilhem.huau@inrae.fr

# A preambule to keep track of when the job has run
echo "Job start: $(date '+%Y-%m-%d %R:%S.%N %Z')"

# We load snakemake module
module load bioinfo/Snakemake/8.20.3

# Dispay snakemake version
snakemake --version

# Run snakemake workflow
snakemake --profile slurm --dry-run

# If the date is displayed, the job reachs the end
echo "Job end: $(date '+%Y-%m-%d %R:%S.%N %Z')"