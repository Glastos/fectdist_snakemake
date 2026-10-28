#!/usr/bin/env Rscript
# R port of tools/make_sheets.py: same options, checks and output.
# Needs the readxl and stringi packages.

USAGE <- "Write config/samples.tsv and config/units.tsv from FASTQ file names.

FASTQ names must follow the provider's pattern:
    {DNA barcode}_{library}_{flowcell}_{lane}_{1|2}.fq.gz
    e.g. GPS263943_MKDN260040585-1A_23M5MYLT4_L4_1.fq.gz
The DNA barcode is looked up in the \"CB ech ADN\" column of the sample table to
get the pool type (\"Type Matrice\") and a readable name (\"Nom Ext\").

Usage:
    tools/make_sheets.R --table ../Table_echantillons_FecDist.xlsm /data/*/*.fq.gz
    find /data -name \"*.fq.gz\" > files.txt
    tools/make_sheets.R --table ../Table_echantillons_FecDist.xlsm --list files.txt

Nothing is written if any file or sample has a problem; all problems are listed.
With --skip-unknown, pools absent from the table are left out (with a warning)
instead of being a problem.

Options:
  FASTQ ...            FASTQ files (R1 and R2)
  --list FILE          text file with one FASTQ path per line
  --table FILE         Table_echantillons_FecDist.xlsm (required)
  --samples FILE       default: config/samples.tsv
  --units FILE         default: config/units.tsv
  --drone-ploidy N     default: 2
  --worker-ploidy N    default: 50
  --force              overwrite existing sheets
  --skip-unknown       leave out pools absent from the table instead of stopping
  -h, --help           show this help
"

FASTQ_NAME <- "^(GPS[0-9]+)_([^_]+)_([^_]+)_(L[0-9]+)_([12])\\.f(ast)?q\\.gz$"
# "Type Matrice" without accents and case -> pipeline pool type
POOL_TYPES <- c("pool male" = "drones", "pool ouvriere" = "workers")


die <- function(...) {
  message(...)
  quit(status = 1)
}


normalise <- function(text) {
  trimws(stringi::stri_trans_general(text, "NFKD; [:Mn:] Remove; Lower"))
}


abspath <- function(path) {
  path <- ifelse(startsWith(path, "/"), path, file.path(getwd(), path))
  vapply(strsplit(path, "/"), function(parts) {
    kept <- character()
    for (part in parts[!parts %in% c("", ".")]) {
      kept <- if (part == "..") head(kept, -1) else c(kept, part)
    }
    paste0("/", paste(kept, collapse = "/"))
  }, character(1))
}


# useBytes keeps accents in "Nom Ext" as UTF-8 whatever the locale
write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(paste(names(df), collapse = "\t"), do.call(paste, c(df, sep = "\t")))
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
}


parse_args <- function(argv) {
  args <- list(fastq = character(), list = NULL, table = NULL,
               samples = "config/samples.tsv", units = "config/units.tsv",
               drone_ploidy = "2", worker_ploidy = "50",
               force = FALSE, skip_unknown = FALSE)
  flags <- c("--force", "--skip-unknown")
  options <- c("--list", "--table", "--samples", "--units", "--drone-ploidy", "--worker-ploidy")
  i <- 1
  while (i <= length(argv)) {
    arg <- argv[i]
    name <- gsub("-", "_", sub("^--", "", arg))
    if (arg %in% c("-h", "--help")) {
      cat(USAGE)
      quit(status = 0)
    } else if (arg %in% flags) {
      args[[name]] <- TRUE
    } else if (arg %in% options) {
      if (i == length(argv)) die(arg, " needs a value (see --help)")
      i <- i + 1
      args[[name]] <- argv[i]
    } else if (startsWith(arg, "-")) {
      die("unknown option ", arg, " (see --help)")
    } else {
      args$fastq <- c(args$fastq, arg)
    }
    i <- i + 1
  }
  if (is.null(args$table)) die("--table is required (see --help)")
  for (name in c("drone_ploidy", "worker_ploidy")) {
    if (!grepl("^[0-9]+$", args[[name]])) die("--", gsub("_", "-", name), " must be an integer")
    args[[name]] <- as.integer(args[[name]])
  }
  args
}


main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  paths <- args$fastq
  if (!is.null(args$list)) {
    lines <- trimws(readLines(args$list, warn = FALSE))
    paths <- c(paths, lines[lines != ""])
  }
  if (length(paths) == 0) die("No FASTQ given")

  problems <- character()

  # One row per FASTQ: sample, library, flowcell, lane, read, path
  m <- regmatches(basename(paths), regexec(FASTQ_NAME, basename(paths)))
  ok <- lengths(m) > 0
  problems <- c(problems, sprintf(
    "name does not match {barcode}_{library}_{flowcell}_{lane}_{1|2}.fq.gz: %s", paths[!ok]))
  fastq <- as.data.frame(matrix(
    as.character(unlist(lapply(m[ok], `[`, 2:6))), ncol = 5, byrow = TRUE,
    dimnames = list(NULL, c("sample", "library", "flowcell", "lane", "read"))))
  fastq$path <- paths[ok]

  # id = sample_library_flowcell_lane, shared by the R1 and R2 of a unit
  fastq$id <- paste(fastq$sample, fastq$library, fastq$flowcell, fastq$lane, sep = "_")
  id_read <- paste(fastq$id, fastq$read)
  dup <- duplicated(id_read)
  problems <- c(problems, sprintf(
    "duplicate R%s for %s: %s and %s", fastq$read[dup], fastq$id[dup],
    fastq$path[match(id_read[dup], id_read)], fastq$path[dup]))
  fastq <- fastq[!dup, ]

  excel <- as.data.frame(suppressMessages(
    readxl::read_excel(args$table, col_types = "text", trim_ws = FALSE)))
  missing <- setdiff(c("CB ech ADN", "Type Matrice", "Nom Ext"), names(excel))
  if (length(missing) > 0) die(args$table, ": missing column(s) ", paste(missing, collapse = ", "))
  excel <- excel[!is.na(excel[["CB ech ADN"]]), ]
  barcodes <- excel[["CB ech ADN"]]
  unknown <- sort(setdiff(fastq$sample, barcodes), method = "radix")
  if (args$skip_unknown) {
    fastq <- fastq[!fastq$sample %in% unknown, ]
  } else {
    problems <- c(problems, sprintf("%s: not in the 'CB ech ADN' column of %s", unknown, args$table))
  }

  # One row per unit with its R1 and R2
  units <- unique(fastq[c("id", "sample", "library", "flowcell", "lane")])
  for (read in c("1", "2")) {
    mates <- fastq[fastq$read == read, ]
    path <- mates$path[match(units$id, mates$id)]
    problems <- c(problems, sprintf("missing R%s for %s", read, units$id[is.na(path)]))
    units[[paste0("fq", read)]] <- path
  }
  problems <- c(problems, sprintf("file not found: %s", fastq$path[!file.exists(fastq$path)]))
  units$unit <- paste(units$flowcell, units$lane, sep = "_")
  sample_unit <- paste(units$sample, units$unit)
  for (clash in unique(sample_unit[duplicated(sample_unit)])) {
    rows <- units[sample_unit == clash, ]
    problems <- c(problems, sprintf(
      "%s: several libraries on %s [%s], unit names would clash",
      rows$sample[1], rows$unit[1], paste(rows$library, collapse = ", ")))
  }

  several_rows <- unique(barcodes[duplicated(barcodes)])
  ploidy <- c(drones = args$drone_ploidy, workers = args$worker_ploidy)
  samples <- data.frame(sample = sort(setdiff(units$sample, unknown), method = "radix"))
  row <- match(samples$sample, barcodes)
  samples$matrix <- excel[["Type Matrice"]][row]
  samples$type <- unname(POOL_TYPES[normalise(samples$matrix)])
  samples$ploidy <- unname(ploidy[samples$type])
  samples$name <- excel[["Nom Ext"]][row]
  samples$name[is.na(samples$name)] <- ""
  in_several_rows <- samples$sample %in% several_rows
  untyped <- !in_several_rows & is.na(samples$type)
  problems <- c(problems,
                sprintf("%s: several rows in %s", samples$sample[in_several_rows], args$table),
                sprintf("%s: unknown Type Matrice '%s'", samples$sample[untyped], samples$matrix[untyped]))
  samples <- samples[!in_several_rows & !untyped, ]

  for (sheet in c(args$samples, args$units)) {
    if (file.exists(sheet) && !args$force) {
      problems <- c(problems, sprintf("%s exists (use --force to overwrite)", sheet))
    }
  }

  if (length(problems) > 0) die("Nothing written:\n  ", paste(problems, collapse = "\n  "))

  write_tsv(samples[c("sample", "type", "ploidy", "name")], args$samples)

  # radix = C-locale order, like Python's sorted()
  units <- units[order(units$sample, units$library, units$flowcell, units$lane, method = "radix"), ]
  units$fq1 <- abspath(units$fq1)
  units$fq2 <- abspath(units$fq2)
  write_tsv(units[c("sample", "unit", "fq1", "fq2", "library")], args$units)

  counts <- sort(table(samples$type), decreasing = TRUE)
  cat(sprintf("%s: %d pools {%s}\n", args$samples, nrow(samples),
              paste(names(counts), counts, sep = ": ", collapse = ", ")))
  cat(sprintf("%s: %d units\n", args$units, nrow(units)))
  if (args$skip_unknown && length(unknown) > 0) {
    message(sprintf("Warning: %d pools not in %s, left out: %s",
                    length(unknown), args$table, paste(unknown, collapse = ", ")))
  }
}


main()
