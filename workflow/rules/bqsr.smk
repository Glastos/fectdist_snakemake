# ------------------------------------------------------------------------------
# Bootstrap BQSR per pool (no known variant set for the species)
#
# Round r, starting from the BAM of round r-1 (round 0 = duplicate-marked BAM):
#   HaplotypeCaller per region -> hard filters -> PASS sites = known sites
#   -> BaseRecalibrator (whole pool, one model per read group) -> ApplyBQSR
# ------------------------------------------------------------------------------


rule bootstrap_call:
    input:
        unpack(round_input),
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bed="results/regions/{region}.bed",
    output:
        vcf=temp("results/bqsr/{sample}/round{round}/{region}.vcf.gz"),
        tbi=temp("results/bqsr/{sample}/round{round}/{region}.vcf.gz.tbi"),
    params:
        java=java_opts,
        ploidy=lambda wildcards: pool_ploidy(wildcards.sample),
    log:
        "logs/bootstrap_call/{sample}/round{round}/{region}.log",
    threads: 2
    resources:
        mem_mb=10000,
        runtime=1440,
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options "{params.java}" HaplotypeCaller'
        " -R {input.ref} -I {input.bam} -L {input.bed} -ploidy {params.ploidy}"
        " --native-pair-hmm-threads {threads} --tmp-dir {resources.tmpdir}"
        " -O {output.vcf} 2> {log}"


rule bootstrap_filter:
    input:
        vcf="results/bqsr/{sample}/round{round}/{region}.vcf.gz",
        tbi="results/bqsr/{sample}/round{round}/{region}.vcf.gz.tbi",
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
    output:
        flagged=temp("results/bqsr/{sample}/round{round}/{region}.flagged.vcf.gz"),
        flagged_tbi=temp(
            "results/bqsr/{sample}/round{round}/{region}.flagged.vcf.gz.tbi"
        ),
        vcf=temp("results/bqsr/{sample}/round{round}/{region}.pass.vcf.gz"),
        tbi=temp("results/bqsr/{sample}/round{round}/{region}.pass.vcf.gz.tbi"),
    params:
        java=java_opts,
        site_filter=config["bqsr"]["site_filter"],
        genotype_filter=config["bqsr"]["genotype_filter"],
    log:
        "logs/bootstrap_filter/{sample}/round{round}/{region}.log",
    resources:
        mem_mb=4000,
        runtime=60,
    conda:
        "../envs/gatk.yaml"
    shell:
        '(gatk --java-options "{params.java}" VariantFiltration'
        " -R {input.ref} -V {input.vcf}"
        ' --filter-expression "{params.site_filter}" --filter-name HQ_fail'
        ' --genotype-filter-expression "{params.genotype_filter}"'
        " --genotype-filter-name GQ_fail"
        " -O {output.flagged}"
        ' && gatk --java-options "{params.java}" SelectVariants'
        " -R {input.ref} -V {output.flagged} --exclude-filtered"
        " -O {output.vcf}) 2> {log}"


rule base_recalibrator:
    input:
        unpack(round_input),
        known=lambda wildcards: expand(
            "results/bqsr/{sample}/round{round}/{region}.pass.vcf.gz",
            sample=wildcards.sample,
            round=wildcards.round,
            region=REGIONS,
        ),
        known_tbi=lambda wildcards: expand(
            "results/bqsr/{sample}/round{round}/{region}.pass.vcf.gz.tbi",
            sample=wildcards.sample,
            round=wildcards.round,
            region=REGIONS,
        ),
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
        bed="results/regions/all.bed",
    output:
        "results/bqsr/{sample}/round{round}.table",
    params:
        java=java_opts,
        known=lambda wildcards, input: " ".join(
            f"--known-sites {vcf}" for vcf in input.known
        ),
    log:
        "logs/base_recalibrator/{sample}/round{round}.log",
    resources:
        mem_mb=10000,
        runtime=1440,
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options "{params.java}" BaseRecalibrator'
        " -R {input.ref} -I {input.bam} -L {input.bed} {params.known}"
        " --tmp-dir {resources.tmpdir} -O {output} 2> {log}"


# Applied once per round of bootstrap.
rule apply_bqsr:
    input:
        unpack(round_input),
        table="results/bqsr/{sample}/round{round}.table",
        ref=REF,
        fai=REF_FAI,
        dict=REF_DICT,
    output:
        bam=temp("results/bqsr/{sample}/round{round}.recal.bam"),
        bai=temp("results/bqsr/{sample}/round{round}.recal.bai"),
    params:
        java=java_opts,
    log:
        "logs/apply_bqsr/{sample}/round{round}.log",
    resources:
        mem_mb=10000,
        runtime=1440,
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options "{params.java}" ApplyBQSR'
        " -R {input.ref} -I {input.bam} --bqsr-recal-file {input.table}"
        " --tmp-dir {resources.tmpdir} -O {output.bam} 2> {log}"


# Hard link keeps the last round once the temp intermediates are removed
rule final_bam:
    input:
        unpack(lambda wildcards: bam_after_round(wildcards.sample, ROUNDS)),
    output:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bai",
    shell:
        "ln -f {input.bam} {output.bam} && ln -f {input.bai} {output.bai}"


TARGETS["bam"] = expand("results/bam/{sample}.bam", sample=SAMPLES)
