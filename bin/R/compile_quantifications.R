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

suppressPackageStartupMessages({
    library(argparse)
    library(tidyverse)
    library(dplyr)
    library(purrr)
    library(stringr)
    library(readr)
    library(tximeta)
    library(AnnotationDbi)
    library(BiocFileCache)
    library(SummarizedExperiment)
    library(HDF5Array)
})

get_ranges <- function(txdb, txome_info) {
    range_path <- paste0("gene", "Rngs-", basename(txome_info$gtf))
    range_counts <- getTximetaBFC() %>%
        list.files(full.names = TRUE) %>%
        str_subset(range_path)
    if (length(range_counts) == 0) {
        ranges <- genes(txdb, single.strand.genes.only = FALSE)
        saveRDS(object = ranges, file = file.path(getTximetaBFC(), paste0(range_path, ".rds")))
    } else {
        ranges <- readRDS(range_counts)
    }

    return(ranges)
}

check_transcripts <- function(assays, ranges) {
    assay_names <- rownames(assays[["counts"]])
    transcripts <- assay_names %in% names(ranges)

    if (!all(transcripts)) {
        if (all(!transcripts)) {
            ## Missing ALL
            stop("All quantified transcripts are missing from gtf")
        } else {
            ## Missing some
            warning(paste0(
                "The annotation is missing some of the transcripts quantified. (",
                sum(!(assay_names %in% names(ranges))), "/", length(assay_names), ")
Dropping missing transcripts...
",
                paste0(
                    "[",
                    head(assay_names[!(assay_names %in% ranges)], 3),
                    ifelse(sum(!(assay_names %in% names(ranges))) > 3, ", ...", ""),
                    "]
",
                    collapse = ", "
                )
            ))
            for (as in names(assays))
            {
                assays[[as]] <- assays[[as]][assay_names %in% names(ranges), , drop = FALSE]
            }
        }
    }
    return(assays)
}

summarize_genes <- function(transcript_quants, txdb) {
    txome_info <- metadata(transcript_quants)$txomeInfo
    tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")
    ranges <- get_ranges(txdb, txome_info)
    seqinfo(ranges) <- seqinfo(transcript_quants)
    tx_tmp <- list(
        abundance = assays(transcript_quants)[["abundance"]],
        counts = assays(transcript_quants)[["counts"]],
        length = assays(transcript_quants)[["length"]],
        countsFromAbundance = "no"
    )
    inferential_indexes <- grep("infRep", assayNames(transcript_quants))
    if (length(inferential_indexes) > 0) {
        replicates <- list(assays(transcript_quants)[inferential_indexes])
        tx_tmp <- c(tx_tmp, infReps = replicates)
    }
    agg_transcript_quants <- summarizeToGene(
        object = tx_tmp,
        tx2gene = tx2gene
    )
    assays <- agg_transcript_quants[c("counts", "abundance", "length")]
    if (length(inferential_indexes) > 0) assays <- c(assays, infReps)
    assays <- check_transcripts(assays, ranges)
    ranges <- ranges[rownames(assays[["counts"]])]
    transcript_ids <- CharacterList(split(tx2gene$TXNAME, tx2gene$GENEID))
    if (all(names(ranges) %in% names(transcript_ids))) {
        transcript_ids <- transcript_ids[names(ranges)]
        mcols(ranges)$tx_ids <- transcript_ids
    }
    if ("hasDuplicate" %in% colnames(mcols(transcript_quants))) {
        stopifnot(all(rownames(transcript_quants) %in% tx2gene[, 1]))
        tx2gene_filter <- tx2gene[match(rownames(transcript_quants), tx2gene[, 1]), ]
        duplicates <- LogicalList(split(mcolds(transcript_quants)$hasDuplicate, tx2gene_filter$GENEID))
        mcols(ranges)$numDupObjects <- sum(duplicates)
    }
    metadata <- metadata(transcript_quants)
    metadata$countsFromAbundance <- "no"
    metadata$level <- "gene"
    gene_quants <- SummarizedExperiment(
        assays = assays,
        rowRanges = ranges,
        metadata = metadata,
        colData = colData(transcript_quants)
    )
    return(gene_quants)
}

read_map <- function(reference_dir, class = TRUE) {
    map <- reference_dir %>%
        list.files(full.names = TRUE) %>%
        str_subset("complete_map.tx2g") %>%
        read_csv(col_names = TRUE, progress = FALSE, show_col_types = FALSE) %>%
        dplyr::select(gene_id, gene_name, gene_biotype) %>%
        distinct() %>%
        mutate(
            gene_biotype = case_when(
                !str_detect(gene_biotype, ",") ~ gene_biotype,
                TRUE ~ unlist(lapply(gene_biotype, function(x) {
                    strsplit(x, ",")[[1]][ifelse(class, 2, 1)]
                }))
            )
        ) %>%
        mutate(
            gene_biotype = case_when(
                startsWith(gene_biotype, "Mt-") | startsWith(gene_name, "Mt-") ~ "Mitochondrial",
                gene_biotype == "protein_coding" ~ "Coding",
                gene_biotype %in% c("lncRNA", "LINE", "SINE", "LTR", "DNA") ~ gene_biotype,
                gene_biotype == "Simple_repeat" ~ "Microsatellite",
                TRUE ~ "Other"
            )
        ) %>%
        return()
}

compile_quants <- function(quants, tx2g, reference, metadata) {
    message("Matching quantifications to the metadata sheet...")
    conditions = read_csv(metadata, col_names=c("prefix", "condition"), progress=F, show_col_types=F)
    samples <- quants %>%
        list.files(full.names = TRUE) %>%
        lapply(function(path) {
            quant <- file.path(path, "quant.sf")
            sample <- basename(path)
            name_list <- lst("files" = quant, "names" = sample)
        }) %>%
        bind_rows() %>%
        as.data.frame() %>%
        apply(X=conditions, MARGIN=1, FUN=function(row, quant) {
            quant %>%
                mutate(prefix=unlist(lapply(names, function(s){substr(s, 1, nchar(row[["prefix"]]))}))) %>%
                inner_join(as.data.frame(t(row)), by="prefix", relationship="one-to-one") %>%
                select(-c("prefix"))
        }, .) %>%
        do.call(rbind, .)

    cache <- paste0(reference, "/.bfccache")
    bfc <- BiocFileCache(cache)
    setTximetaBFC(cache)

    message("Compiling quantifications...")
    transcript_quants <- tximeta::tximeta(coldata = samples, type = "salmon", skipMeta = FALSE)

    message("Loading TxDb...")
    txdb <- bfcinfo(bfc) %>%
        as.data.frame() %>%
        filter(endsWith(rname, "complete_annotation.gtf.gz")) %>%
        slice_head(n = 1) %>%
        dplyr::select(rpath) %>%
        pull() %>%
        loadDb()

    message("Summarizing genes...")
    gene_quants <- summarize_genes(transcript_quants, txdb)
    rowData(gene_quants) <- rowData(gene_quants) %>%
        as.data.frame() %>%
        left_join(tx2g, multiple = "any", by = "gene_id")
    
    message("Adding normalized counts")
    assays(gene_quants)$normalized = DESeqDataSet(gene_quants) %>%
        estimateSizeFactors() %>%
        counts(noramlize=TRUE)

    message("Saving gene quantifications...")
    saveHDF5SummarizedExperiment(gene_quants, dir = "counts", replace = TRUE)
}

compile_af <- function(quants, tx2g) {
    suppressPackageStartupMessages({
        library(fishpond)
        library(SingleCellExperiment)
    })

    message("Importing quantifications...")
    gene_quants <- quants %>%
        list.files(full.names = TRUE) %>%
        lapply(function(quant) {
            message(paste0("Loading ", basename(quant), "..."))
            raw <- loadFry(quant, outputFormat = "scRNA")
            counts <- SingleCellExperiment(
                assays(raw),
                colData = colData(raw) %>%
                    as.data.frame() %>%
                    mutate(sample = quant),
                rowData = rowData(raw) %>%
                    as.data.frame() %>%
                    left_join(tx2g, by = c("gene_ids" = "gene_id"), multiple = "any"),
                metadata = metadata(raw),
                checkDimnames = TRUE
            )
            return(counts)
        }) %>%
        do.call(args = ., what = "cbind")
    message("Chunking and storing data...")
    saveHDF5SummarizedExperiment(
        gene_quants,
        dir = "counts",
        replace = TRUE, chunkdim = c(1e4, 1e4)
    )
}

main <- function() {
    parser <- ArgumentParser()
    parser$add_argument(
        "-q", "--quants",
        action = "store",
        help = "The quantification directory"
    )
    parser$add_argument(
        "-r", "--reference",
        action = "store",
        help = "The CREATE reference"
    )
    parser$add_argument(
        "-m", "--metadata",
        action = "store",
        help = "The metadata samplesheet with prefixes and conditions"
    )
    parser$add_argument(
        "-s", "--splintr",
        action = "store_true", default = FALSE, required = FALSE,
        help = "The split gene experiment for splintr"
    )
    args <- parser$parse_args()

    tx2g <- read_map(args$reference)

    if (args$splintr) {
        compile_af(args$quants, tx2g, args$metadata)
    } else {
        compile_quants(args$quants, tx2g, args$reference, args$metadata)
    }
}
main()
