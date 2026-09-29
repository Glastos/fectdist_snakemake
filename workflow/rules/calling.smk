# ------------------------------------------------------------------------------
# Per-pool gVCF, ploidy from samples.tsv
# ------------------------------------------------------------------------------


rule haplotype_caller_gvcf:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bai",
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bed="results/regions/all.bed",
    output:
        gvcf="results/gvcf/{sample}.g.vcf.gz",
        tbi="results/gvcf/{sample}.g.vcf.gz.tbi",
    params:
        java=java_opts,
        ploidy=lambda wildcards: pool_ploidy(wildcards.sample),
    log:
        "logs/haplotype_caller_gvcf/{sample}.log",
    threads: 2
    resources:
        mem_mb=32000,
        runtime=4320,
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options "{params.java}" HaplotypeCaller'
        " -R {input.ref} -I {input.bam} -L {input.bed} -ERC GVCF"
        " -ploidy {params.ploidy} --native-pair-hmm-threads {threads}"
        " --tmp-dir {resources.tmpdir} -O {output.gvcf} 2> {log}"


TARGETS["gvcf"] = expand("results/gvcf/{sample}.g.vcf.gz", sample=SAMPLES) + BQSR_CHECKS
