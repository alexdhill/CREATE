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

quick_parse <- function(info, out) {
    `%>%` <- tidyr::`%>%`
    attrs <- info %>%
        readLines() %>%
        stringr::str_split(";") %>%
        lapply(function(s) stringr::str_split(s, ' "')) %>%
        lapply(function(s) {
            s %>%
                lapply(function(a) {
                    r <- a[[2]] %>%
                        stringr::str_replace_all(";", "") %>%
                        stringr::str_replace_all('"', "")
                    names(r) <- a[[1]]
                    r
                }) %>%
                do.call(what = "c")
        }) %>%
        dplyr::bind_rows() %>%
        dplyr::select(dplyr::any_of(
            "transcript_id",
            "gene_id",
            "gene_name",
            "gene_biotype",
            "source_gene",
            "source_transcript"
        )) %>%
        dplyr::distinct() %>%
        readr::write_csv(., out, col_names = TRUE)
}

parse_info <- function(info, out) {
    `%>%` <- tidyr::`%>%`
    data.frame(attrs = readLines(info)) %>%
        dplyr::mutate(uid = dplyr::row_number()) %>%
        tidyr::separate_longer_delim(cols = attrs, delim = "; ") %>%
        dplyr::filter(attrs != "") %>%
        tidyr::separate_wider_regex(
            cols = attrs,
            patterns = c(key = "[^ ]+", " ", value = ".*")
        ) %>%
        dplyr::mutate(value = stringr::str_remove_all(value, '"')) %>%
        tidyr::pivot_wider(names_from = key, values_from = value) %>%
        dplyr::distinct() %>%
        dplyr::select(dplyr::any_of(c(
            "transcript_id",
            "gene_id",
            "gene_biotype",
            "gene_name",
            "source_gene",
            "source_transcript"
        ))) %>%
        readr::write_csv(., out, col_names = TRUE)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop(
            "Usage: Rscript parse_transcript_info.R <transcript_info> <output_file>"
        )
    }

    transcript_info <- args[1]
    output_file <- args[2]

    parse_info(transcript_info, output_file)
}
main()
