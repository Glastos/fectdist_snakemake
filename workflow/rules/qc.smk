# ------------------------------------------------------------------------------
# QC per pool (results/{sample}/qc/), displayed in results/all/qc/qc_report.html
# (workflow/report/qc_report.qmd)
#   samtools stats and mosdepth on the final BAM, BQSR tables (bqsr.smk),
#   allele balance of the drone pools from their gVCF (with the gvcf target),
#   pool type and ploidy
# ------------------------------------------------------------------------------

DRONES = samples.loc[samples["type"] == "drones", "sample"].tolist()

# Allele balance of the drone pools only with the gvcf target: `targets: [qc]`
# alone reports on the pools up to their BAM, without calling the gVCFs
ALLELE_BALANCE = []
if "gvcf" in config["targets"]:
    ALLELE_BALANCE = expand(
        "results/{sample}/qc/{sample}.allele_balance.{file}",
        sample=DRONES,
        file=["sites.tsv.gz", "summary.tsv"],
    )


rule samtools_stats:
    input:
        bam="results/{sample}/{sample}.bam",
        bai="results/{sample}/{sample}.bai",
    output:
        "results/{sample}/qc/{sample}.samtools.stats",
    log:
        "logs/{sample}/samtools_stats/{sample}.log",
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
        bam="results/{sample}/{sample}.bam",
        bai="results/{sample}/{sample}.bai",
    output:
        summary="results/{sample}/qc/{sample}.mosdepth.summary.txt",
        dist="results/{sample}/qc/{sample}.mosdepth.global.dist.txt",
        windows="results/{sample}/qc/{sample}.regions.bed.gz",
    params:
        prefix="results/{sample}/qc/{sample}",
        min_mapq=config["qc"]["mosdepth_min_mapq"],
        window=config["qc"]["mosdepth_window"],
    log:
        "logs/{sample}/mosdepth/{sample}.log",
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
        gvcf="results/{sample}/{sample}.g.vcf.gz",
        tbi="results/{sample}/{sample}.g.vcf.gz.tbi",
    output:
        sites="results/{sample}/qc/{sample}.allele_balance.sites.tsv.gz",
        summary="results/{sample}/qc/{sample}.allele_balance.summary.tsv",
    params:
        min_depth=config["qc"]["allele_balance"]["min_depth"],
        max_depth_ratio=config["qc"]["allele_balance"]["max_depth_ratio"],
        max_sites=config["qc"]["allele_balance"]["max_sites"],
        seed=config["qc"]["allele_balance"]["seed"],
    log:
        "logs/{sample}/drone_allele_balance/{sample}.log",
    resources:
        mem_mb=4000,
        runtime=120,
    conda:
        "../envs/qc.yaml"
    script:
        "../scripts/allele_balance.py"


# Type and ploidy of the pool, kept in its folder: the report reads them there,
# for the pools of any sample sheet (rewritten if the sheet changes them)
rule pool_info:
    output:
        "results/{sample}/qc/{sample}.pool.tsv",
    params:
        type=lambda w: samples.loc[w.sample, "type"],
        ploidy=lambda w: samples.loc[w.sample, "ploidy"],
    run:
        with open(output[0], "w") as out:
            out.write("sample\ttype\tploidy\n")
            out.write(f"{wildcards.sample}\t{params.type}\t{params.ploidy}\n")


# Rendered from a copy in results/all/qc/, so the HTML lands next to the QC files.
# Builds the QC of the pools of the sample sheet, but the report shows every pool
# folder of results/ with QC files, from earlier runs with other sheets too.
rule qc_report:
    input:
        qmd=f"{workflow.basedir}/report/qc_report.qmd",
        pool=expand("results/{sample}/qc/{sample}.pool.tsv", sample=SAMPLES),
        stats=expand("results/{sample}/qc/{sample}.samtools.stats", sample=SAMPLES),
        depth=expand(
            "results/{sample}/qc/{sample}.{file}",
            sample=SAMPLES,
            file=["mosdepth.summary.txt", "mosdepth.global.dist.txt", "regions.bed.gz"],
        ),
        bqsr=expand(
            "results/{sample}/bqsr/round{round}.table",
            sample=SAMPLES,
            round=range(1, ROUNDS + 1),
        )
        + BQSR_CHECKS,
        allele_balance=ALLELE_BALANCE,
    output:
        "results/all/qc/qc_report.html",
    log:
        "logs/all/qc_report/qc_report.log",
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/report.yaml"
    shell:
        "cp {input.qmd} results/all/qc/qc_report.qmd"
        " && quarto render results/all/qc/qc_report.qmd --to html > {log} 2>&1"


TARGETS["qc"] = ["results/all/qc/qc_report.html"]
