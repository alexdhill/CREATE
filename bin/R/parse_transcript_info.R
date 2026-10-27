## REQUIRED NOTICE: Copyright (c) 2020-2023, Regents of the University of California
## All rights reserved. https://polyformproject.org/licenses/noncommercial/1.0.0
##
## This software was developed by the Daniel Kim lab at the University of California, Santa Cruz.
## Authors: Roman E. Reggiardo, Vikas Peddu, Alex D. Hill
##
## The licensor grants you a copyright license for the software to do everything you might do with
## the software that would otherwise infringe the licensor’s copyright in it for any permitted
## purpose.
##
## As far as the law allows, the software comes as is, without any warranty or condition, and the
## licensor will not be liable to you for any damages arising out of these terms or the use or
## nature of the software, under any kind of legal claim.

attributes <-
    c(
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_biotype",
        "source_gene",
        "source_transcript"
    )

log_ranks <- c(FATAL = 1L, ERROR = 1L, WARN = 2L, INFO = 3L, DEBUG = 4L)

.log <- function(level, message) {
    stopifnot(level %in% names(log_ranks))
    if (log_ranks[[level]] <= getOption("create.verbosity", 1L)) {
        stamp <- format(Sys.time(), "%d-%m-%Y %H:%M:%S-%Z")
        cat(sprintf("[%s] [%s] %s\n", stamp, level, message), file = stderr())
    }
    invisible(NULL)
}

.read_options <- function(args) {
    stopifnot(is.character(args))
    verbosity <- 1L
    paths <- character(0)
    consumed <- FALSE
    for (index in seq_along(args)) {
        token <- args[index]
        level <- if (index < length(args)) args[index + 1L] else ""
        if (consumed) {
            consumed <- FALSE
        } else if (token == "-v" || token == "--verbose") {
            consumed <- grepl("^[1-4]$", level)
            verbosity <- if (consumed) as.integer(level) else 1L
        } else {
            paths <- c(paths, token)
        }
    }
    list(verbosity = verbosity, paths = paths)
}

.read_attribute_lines <- function(info) {
    stopifnot(is.character(info), length(info) == 1L, file.exists(info))
    lines <- readr::read_lines(info, progress = FALSE)
    .log("DEBUG", sprintf("read %d attribute records from %s", length(lines), info))
    unique(lines[nzchar(lines)])
}

.extract_attribute <- function(lines, attribute) {
    stopifnot(is.character(lines), is.character(attribute))
    pattern <- paste0('(?:^|[;\\s])', attribute, ' "([^"]*)"')
    candidates <- stringr::str_detect(lines, stringr::fixed(attribute))
    values <- rep(NA_character_, length(lines))
    values[candidates] <- stringr::str_extract(lines[candidates], pattern, group = 1L)
    .log("DEBUG", sprintf(
        "matched %s on %d of %d records",
        attribute, sum(!is.na(values)), length(values)
    ))
    values
}

.build_attribute_table <- function(lines) {
    stopifnot(is.character(lines), length(lines) > 0L)
    columns <- lapply(attributes, function(a) .extract_attribute(lines, a))
    names(columns) <- attributes
    present <- vapply(columns, function(v) any(!is.na(v)), logical(1L))
    if (!any(present)) {
        .log("FATAL", "no known transcript attributes found in input")
        stop("no known transcript attributes found in input")
    }
    if (!all(present)) {
        .log("WARN", sprintf(
            "absent attributes dropped: %s",
            paste(attributes[!present], collapse = ", ")
        ))
    }
    dplyr::as_tibble(columns[present])
}

quick_parse <- function(info, out) {
    stopifnot(is.character(out), length(out) == 1L, nzchar(out))
    .log("INFO", sprintf("parsing transcript attributes from %s", info))
    lines <- .read_attribute_lines(info)
    if (length(lines) == 0L) {
        .log("FATAL", sprintf("no attribute records in %s", info))
        stop("empty transcript info file")
    }
    .log("INFO", sprintf("parsing %d unique attribute records", length(lines)))
    attrs <- dplyr::distinct(.build_attribute_table(lines))
    .log("INFO", sprintf("writing %d transcript records to %s", nrow(attrs), out))
    readr::write_csv(attrs, out, col_names = TRUE, progress = FALSE)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    options <- .read_options(args)
    base::options(create.verbosity = options$verbosity)
    if (length(options$paths) < 2) {
        .log("FATAL", "missing transcript info or output path")
        stop(
            "Usage: Rscript parse_transcript_info.R [-v 1-4] <transcript_info> <output_file>"
        )
    }

    transcript_info <- options$paths[1]
    output_file <- options$paths[2]

    quick_parse(transcript_info, output_file)
}
main()
