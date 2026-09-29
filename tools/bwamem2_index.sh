#!/bin/sh
#SBATCH --cpus-per-task=8 
#SBATCH --mem=16392

module purge
module load bioinfo/bwa-mem2/2.2.1
module load bioinfo/samtools/1.23
module load devel/python/Python-3.11.1
module load bioinfo/GATK/4.6.2.0


bwa-mem2 index resources/reference_HAv3_1/GCF_003254395.2_Amel_HAv3.1_genomic.fna

samtools faidx resources/reference_HAv3_1/GCF_003254395.2_Amel_HAv3.1_genomic.fna

gatk CreateSequenceDictionary -R resources/reference_HAv3_1/GCF_003254395.2_Amel_HAv3.1_genomic.fna
