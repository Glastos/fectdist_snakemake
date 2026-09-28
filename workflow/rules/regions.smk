# ------------------------------------------------------------------------------
# BED files of the regions (see make_regions in common.smk)
# ------------------------------------------------------------------------------


rule region_bed:
    output:
        "results/regions/{region}.bed",
    params:
        intervals=lambda wildcards: REGIONS[wildcards.region],
    run:
        with open(output[0], "w") as bed:
            for contig, start, end in params.intervals:
                bed.write(f"{contig}\t{start}\t{end}\n")


rule all_regions_bed:
    output:
        "results/regions/all.bed",
    params:
        intervals=[interval for region in REGIONS.values() for interval in region],
    run:
        with open(output[0], "w") as bed:
            for contig, start, end in params.intervals:
                bed.write(f"{contig}\t{start}\t{end}\n")
