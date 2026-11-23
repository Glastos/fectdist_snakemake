#!/usr/bin/env Rscript
# Write the sample sheet of the pipeline from AVITI (Element Biosciences) FASTQ
# file names, for individual (not pooled) diploid samples.
#
# FASTQ names must follow the provider's pattern:
#     {ID}_{index}-{flowcell}_{lane}_{R1|R2}.fastq.gz
# {index} is the index pair (or a single index). The names give no library ID,
# so the index stands for the library: the lanes of a library share its index,
# MarkDuplicates looks for duplicates within each library.
# No sample table is read: every sample gets the type, ploidy and platform set
# below. The sheet has one row per read group (sample x library x lane).
#
# Edit the settings below, then run from the pipeline directory:
#     tools/make_sheets_AVITI.R    or, in R:  source("tools/make_sheets_AVITI.R")
# Nothing is written if any file has a problem; all problems are listed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
})

# ------------------------------------------------------------------------------
# Settings (paths relative to the pipeline directory)
# ------------------------------------------------------------------------------

# One FASTQ (R1 or R2) per line, e.g. find /data -name "*.fastq.gz" > tools/list_fastq_aviti.txt
# The paths are written as they are in the sheet
FASTQ_LIST <- "tools/list_fastq_aviti.txt"
OUTPUT <- "config/samples_aviti.tsv"
OVERWRITE <- FALSE      # replace an existing OUTPUT
TYPE <- "individuals"
PLOIDY <- 2L
PLATFORM <- "ELEMENT"   # PL of the read groups (SAM specification value for AVITI)

FASTQ_NAME <- "^([^_]+)_([ACGTN]+(?:-[ACGTN]+)*)-([^_-]+)_(L[0-9]+)_R([12])\\.f(ast)?q\\.gz$"

# ------------------------------------------------------------------------------

paths <- trimws(readLines(FASTQ_LIST, warn = FALSE))
paths <- paths[paths != ""]
if (length(paths) == 0) stop(FASTQ_LIST, ": no FASTQ", call. = FALSE)

# One row per FASTQ; id = sample_index_flowcell_lane, shared by R1 and R2
m <- str_match(basename(paths), FASTQ_NAME)
fastq <- tibble(path = paths, sample = m[, 2], library = m[, 3],
                unit = str_c(m[, 4], "_", m[, 5]), read = factor(m[, 6], c("1", "2")),
                id = str_c(sample, "_", library, "_", unit))
unmatched <- filter(fastq, is.na(sample))
fastq <- filter(fastq, !is.na(sample))
duplicates <- fastq |>
  summarise(n = n(), paths = str_flatten(path, " and "), .by = c(id, read)) |>
  filter(n > 1)

# One row per read group, with its R1 and R2
units <- fastq |>
  distinct(id, read, .keep_all = TRUE) |>
  pivot_wider(id_cols = c(id, sample, library, unit), names_from = read,
              names_prefix = "fq", names_expand = TRUE, values_from = path)

# Read groups are named {sample}.{unit} (unit = {flowcell}_{lane}): two
# libraries of one sample on the same lane would share a name
shared_lanes <- units |>
  summarise(n = n(), libraries = str_flatten(library, ", "), .by = c(sample, unit)) |>
  filter(n > 1)

problems <- c(
  sprintf("name does not match {ID}_{index}-{flowcell}_{lane}_{R1|R2}.fastq.gz: %s",
          unmatched$path),
  sprintf("duplicate R%s for %s: %s", duplicates$read, duplicates$id, duplicates$paths),
  sprintf("missing R1 for %s", units$id[is.na(units$fq1)]),
  sprintf("missing R2 for %s", units$id[is.na(units$fq2)]),
  sprintf("%s: several libraries on %s (%s)", shared_lanes$sample, shared_lanes$unit,
          shared_lanes$libraries),
  sprintf("file not found: %s", fastq$path[!file.exists(fastq$path)]),
  if (file.exists(OUTPUT) && !OVERWRITE) sprintf("%s exists (set OVERWRITE <- TRUE)", OUTPUT)
)
# message(): stop() would truncate a long list
if (length(problems) > 0) {
  message(str_flatten(problems, "\n"))
  stop("Nothing written, ", length(problems), " problem(s) above", call. = FALSE)
}

sheet <- units |>
  mutate(type = TYPE, ploidy = PLOIDY, platform = PLATFORM) |>
  arrange(sample, library, unit) |>
  select(sample, type, ploidy, unit, fq1, fq2, library, platform)
write_tsv(sheet, OUTPUT, na = "")

cat(sprintf("%s: %d samples (%s, ploidy %d, %s), %d read groups\n", OUTPUT,
            n_distinct(sheet$sample), TYPE, PLOIDY, PLATFORM, nrow(sheet)))
