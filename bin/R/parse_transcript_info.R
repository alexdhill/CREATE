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

suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(stringr))
suppressMessages(library(magrittr))

filter_info <- function(i) {
    grepl("transcript_id|gene_id|gene_biotype|gene_name", i) %>%
        return()
}

parse_info <- function(info, out) {
    readLines(info) %>%
        lapply(., function(i) {
            unlist(strsplit(i, ';')) %>%
            .[grepl(pattern="transcript_id|gene_id|gene_biotype|gene_name", x=.)] %>%
            gsub("^\\s+|\"", "", .) %>%
            set_names(lapply(., function(x) {strsplit(x, "\\s+")[[1]][1]})) %>%
            lapply(function(x) strsplit(x, "\\s+")[[1]][2]) %>%
            return()
        }) %>%
        bind_rows() %>%
        mutate(gene_name=case_when(is.na(gene_name)~gene_id, .default=gene_name)) %>%
        distinct() %>%
        write_csv(., out, col_names = TRUE)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 2) {
        stop("Usage: Rscript parse_transcript_info.R <transcript_info> <output_file>")
    }

    transcript_info <- args[1]
    output_file <- args[2]

    parse_info(transcript_info, output_file)
}
main()