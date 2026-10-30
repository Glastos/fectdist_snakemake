#!/usr/bin/env Rscript
# Write the sample sheet of the pipeline from FASTQ file names.
#
# FASTQ names must follow the provider's pattern:
#     {DNA barcode}_{library}_{flowcell}_{lane}_{1|2}.fq.gz
#     e.g. GPS263943_MKDN260040585-1A_23M5MYLT4_L4_1.fq.gz
# The DNA barcode is looked up in the "CB ech ADN" column of the sample table to
# get the pool type ("Type Matrice") and a readable name ("Nom Ext").
# The sheet has one row per read group (pool x library x lane), the pool columns
# (type, ploidy, name) repeated on every row of the pool.
#
# Edit the settings below, then run from the pipeline directory:
#     tools/make_sheets.R          or, in R:  source("tools/make_sheets.R")
# Nothing is written if any file or pool has a problem; all problems are listed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(readr)
})

# ------------------------------------------------------------------------------
# Settings (paths relative to the pipeline directory)
# ------------------------------------------------------------------------------

TABLE <- "../Table_echantillons_FecDist.xlsm"
# One FASTQ (R1 or R2) per line, e.g. find /data -name "*.fq.gz" > config/files.txt
# The paths are written as they are in the sheet
FASTQ_LIST <- "config/files.txt"
OUTPUT <- "config/samples.tsv"
OVERWRITE <- FALSE      # replace an existing OUTPUT
SKIP_UNKNOWN <- FALSE   # leave out (with a warning) pools absent from TABLE instead of stopping
PLOIDY <- c(drones = 2L, workers = 50L)

# "Type Matrice" in lower case -> pool type (pool male, pool mâle, pool ouvrière;
# accents as \u escapes so that they match in any locale)
POOL_TYPES <- setNames(c("drones", "drones", "workers"),
                       c("pool male", "pool m\u00e2le", "pool ouvri\u00e8re"))
FASTQ_NAME <- "^(GPS[0-9]+)_([^_]+)_([^_]+)_(L[0-9]+)_([12])\\.f(ast)?q\\.gz$"

# ------------------------------------------------------------------------------

paths <- trimws(readLines(FASTQ_LIST, warn = FALSE))
paths <- paths[paths != ""]
if (length(paths) == 0) stop(FASTQ_LIST, ": no FASTQ", call. = FALSE)

# One row per FASTQ; id = sample_library_flowcell_lane, shared by R1 and R2
m <- str_match(basename(paths), FASTQ_NAME)
fastq <- tibble(path = paths, sample = m[, 2], library = m[, 3],
                unit = str_c(m[, 4], "_", m[, 5]), read = factor(m[, 6], c("1", "2")),
                id = str_c(sample, "_", library, "_", unit))
unmatched <- filter(fastq, is.na(sample))
fastq <- filter(fastq, !is.na(sample))
duplicates <- fastq |>
  summarise(n = n(), paths = str_flatten(path, " and "), .by = c(id, read)) |>
  filter(n > 1)

pools <- suppressMessages(read_excel(TABLE, col_types = "text")) |>
  filter(!is.na(`CB ech ADN`)) |>
  transmute(sample = `CB ech ADN`, matrix = `Type Matrice`,
            type = unname(POOL_TYPES[str_to_lower(str_trim(matrix))]),
            ploidy = unname(PLOIDY[type]), name = coalesce(`Nom Ext`, ""))
unknown <- sort(setdiff(fastq$sample, pools$sample))
if (SKIP_UNKNOWN) fastq <- filter(fastq, !sample %in% unknown)

# One row per read group, with its R1 and R2
units <- fastq |>
  distinct(id, read, .keep_all = TRUE) |>
  pivot_wider(id_cols = c(id, sample, library, unit), names_from = read,
              names_prefix = "fq", names_expand = TRUE, values_from = path)

used <- filter(pools, sample %in% units$sample)
several_rows <- unique(used$sample[duplicated(used$sample)])
untyped <- filter(used, is.na(type), !sample %in% several_rows)

problems <- c(
  sprintf("name does not match {barcode}_{library}_{flowcell}_{lane}_{1|2}.fq.gz: %s",
          unmatched$path),
  sprintf("duplicate R%s for %s: %s", duplicates$read, duplicates$id, duplicates$paths),
  sprintf("missing R1 for %s", units$id[is.na(units$fq1)]),
  sprintf("missing R2 for %s", units$id[is.na(units$fq2)]),
  sprintf("file not found: %s", fastq$path[!file.exists(fastq$path)]),
  if (!SKIP_UNKNOWN) sprintf("%s: not in the 'CB ech ADN' column of %s", unknown, TABLE),
  sprintf("%s: several rows in %s", several_rows, TABLE),
  sprintf("%s: unknown Type Matrice '%s'", untyped$sample, untyped$matrix),
  if (file.exists(OUTPUT) && !OVERWRITE) sprintf("%s exists (set OVERWRITE <- TRUE)", OUTPUT)
)
# message(): stop() would truncate a long list
if (length(problems) > 0) {
  message(str_flatten(problems, "\n"))
  stop("Nothing written, ", length(problems), " problem(s) above", call. = FALSE)
}

sheet <- units |>
  inner_join(pools, by = "sample", relationship = "many-to-one") |>
  arrange(sample, library, unit) |>
  select(sample, type, ploidy, name, unit, fq1, fq2, library)
write_tsv(sheet, OUTPUT, na = "")

types <- count(distinct(sheet, sample, type), type)
cat(sprintf("%s: %d pools (%s), %d read groups\n", OUTPUT, n_distinct(sheet$sample),
            str_flatten(str_c(types$type, ": ", types$n), ", "), nrow(sheet)))
if (SKIP_UNKNOWN && length(unknown) > 0) {
  message(sprintf("Warning: %d pools not in %s, left out: %s", length(unknown), TABLE,
                  str_flatten(unknown, ", ")))
}
