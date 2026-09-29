# ------------------------------------------------------------------------------
# Per read group (unit): bwa-mem2 mem | samtools sort
# Per pool: MarkDuplicates on all its read groups at once (merge + dedup)
#
# Follows GATK "How should I pre-process data from multiplexed sequencing and
# multi-library designs?" (article 360035889471)
# ------------------------------------------------------------------------------


rule read_group:
    input:
        lambda wildcards: get_unit(wildcards)["fq1"],
    output:
        "results/mapping/{sample}/{unit}.rg.txt",
    params:
        library=lambda wildcards: get_unit(wildcards)["library"],
    log:
        "logs/read_group/{sample}/{unit}.log",
    script:
        "../scripts/read_group.py"


rule bwa_mem:
    input:
        fq1=lambda wildcards: get_unit(wildcards)["fq1"],
        fq2=lambda wildcards: get_unit(wildcards)["fq2"],
        rg="results/mapping/{sample}/{unit}.rg.txt",
        ref=REF,
        idx=REF_BWA,
    output:
        temp("results/mapping/{sample}/{unit}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}/{unit}.log",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=720,
    conda:
        "../envs/mapping.yaml"
    shell:
        '(bwa-mem2 mem -M -t {threads} -R "$(cat {input.rg})" {input.ref} {input.fq1} {input.fq2}'
        " | samtools sort -@ 2 -m 1G -T {output}.tmp -o {output} -) 2> {log}"


rule mark_duplicates:
    input:
        lambda wildcards: expand(
            "results/mapping/{sample}/{unit}.sorted.bam",
            sample=wildcards.sample,
            unit=sample_units(wildcards.sample),
        ),
    output:
        bam=temp("results/mapping/{sample}.dedup.bam"),
        bai=temp("results/mapping/{sample}.dedup.bai"),
        metrics="results/qc/mark_duplicates/{sample}.metrics.txt",
    params:
        java=java_opts,
        inputs=lambda wildcards, input: " ".join(f"-I {bam}" for bam in input),
    log:
        "logs/mark_duplicates/{sample}.log",
    resources:
        mem_mb=16000,
        runtime=720,
    conda:
        "../envs/gatk.yaml"
    shell:
        'gatk --java-options "{params.java}" MarkDuplicates {params.inputs}'
        " -O {output.bam} -M {output.metrics} --CREATE_INDEX true"
        " --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 2040 --TMP_DIR {resources.tmpdir}"
        " 2> {log}"
