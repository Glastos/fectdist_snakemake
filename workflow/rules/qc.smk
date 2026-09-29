# ------------------------------------------------------------------------------
# QC per pool, displayed in results/qc/qc_report.html (workflow/report/qc_report.qmd)
#   samtools stats and mosdepth on the final BAM, BQSR tables (bqsr.smk),
#   allele balance of the drone pools from their gVCF
# ------------------------------------------------------------------------------

DRONES = samples.loc[samples["type"] == "drones", "sample"].tolist()


rule samtools_stats:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bai",
    output:
        "results/qc/samtools_stats/{sample}.stats",
    log:
        "logs/samtools_stats/{sample}.log",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../envs/qc.yaml"
    shell:
        "samtools stats -@ {threads} {input.bam} > {output} 2> {log}"


rule mosdepth:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bai",
    output:
        summary="results/qc/mosdepth/{sample}.mosdepth.summary.txt",
        dist="results/qc/mosdepth/{sample}.mosdepth.global.dist.txt",
        windows="results/qc/mosdepth/{sample}.regions.bed.gz",
    params:
        prefix="results/qc/mosdepth/{sample}",
        min_mapq=config["qc"]["mosdepth_min_mapq"],
        window=config["qc"]["mosdepth_window"],
    log:
        "logs/mosdepth/{sample}.log",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../envs/qc.yaml"
    shell:
        "mosdepth --threads {threads} --no-per-base --mapq {params.min_mapq}"
        " --by {params.window} {params.prefix} {input.bam} 2> {log}"


rule drone_allele_balance:
    input:
        gvcf="results/gvcf/{sample}.g.vcf.gz",
        tbi="results/gvcf/{sample}.g.vcf.gz.tbi",
    output:
        sites="results/qc/allele_balance/{sample}.sites.tsv.gz",
        summary="results/qc/allele_balance/{sample}.summary.tsv",
    params:
        min_depth=config["qc"]["allele_balance"]["min_depth"],
        max_depth_ratio=config["qc"]["allele_balance"]["max_depth_ratio"],
        max_sites=config["qc"]["allele_balance"]["max_sites"],
        seed=config["qc"]["allele_balance"]["seed"],
    log:
        "logs/drone_allele_balance/{sample}.log",
    resources:
        mem_mb=4000,
        runtime=120,
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/allele_balance.py"


# Pools and their type, read by the report
rule qc_samples:
    output:
        "results/qc/samples.tsv",
    params:
        rows=samples[["sample", "type", "ploidy"]].values.tolist(),
    run:
        with open(output[0], "w") as out:
            out.write("sample\ttype\tploidy\n")
            for sample, pool_type, ploidy in params.rows:
                out.write(f"{sample}\t{pool_type}\t{ploidy}\n")


# Rendered from a copy in results/qc/, so the HTML lands next to the QC files
rule qc_report:
    input:
        qmd=f"{workflow.basedir}/report/qc_report.qmd",
        samples="results/qc/samples.tsv",
        stats=expand("results/qc/samtools_stats/{sample}.stats", sample=SAMPLES),
        depth=expand(
            "results/qc/mosdepth/{sample}.{file}",
            sample=SAMPLES,
            file=["mosdepth.summary.txt", "mosdepth.global.dist.txt", "regions.bed.gz"],
        ),
        bqsr=expand(
            "results/bqsr/{sample}/round{round}.table",
            sample=SAMPLES,
            round=range(1, ROUNDS + 1),
        )
        + BQSR_CHECKS,
        allele_balance=expand(
            "results/qc/allele_balance/{sample}.{file}",
            sample=DRONES,
            file=["sites.tsv.gz", "summary.tsv"],
        ),
    output:
        "results/qc/qc_report.html",
    log:
        "logs/qc_report/qc_report.log",
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/report.yaml"
    shell:
        "cp {input.qmd} results/qc/qc_report.qmd"
        " && quarto render results/qc/qc_report.qmd --to html > {log} 2>&1"


TARGETS["qc"] = ["results/qc/qc_report.html"]
