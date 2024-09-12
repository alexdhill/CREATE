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
    library(tidyr)
    library(argparse)
    library(tximeta)
    library(rjson)
    library(txdbmaker)
    library(GenomicFeatures)
    library(BiocFileCache)
})

link_transcriptome <- function(annotation, transcripts, index)
{
    build = strsplit(annotation, '_')[[1]][1]
    out_json = paste0(build, "_complete_txome.json")
    message("Making linked txome...")
    makeLinkedTxome(
        indexDir = index,
        fasta = transcripts,
        gtf = annotation,
        source = "CREATE",
        organism = "Homo sapiens",
        genome = ifelse(build=="T2Tv2", "CHM13", "GRCh38"),
        release = build,
        jsonFile = out_json,
        write = TRUE
    )
}

build_txdb <- function(bfc, annotation)
{
    txdb = txdbmaker::makeTxDbFromGFF(
        annotation,
        dataSource="CREATE",
        organism="Homo sapiens"
    )
    saveDb(txdb, bfcnew(bfc, rname=annotation, ext=".sqlite"))
}

main <- function()
{
    parser = ArgumentParser()
    parser$add_argument(
        "-a", "--annotation", action="store",
        help="The reference annotation"
    )
    parser$add_argument(
        "-t", "--transcripts", action="store",
        help="The reference transcriptome"
    )
    parser$add_argument(
        "-i", "--index", action="store",
        help="The salmon index/digest"
    )
    # parser$add_argument(
    #     "-p", "--prefix", action="store",
    #     help="The type of txome being linked"
    # )
    args = parser$parse_args()
    
    message("Making file cache...")
    dir.create(".bfccache")
    bfc = BiocFileCache(".bfccache")
    setTximetaBFC(".bfccache")

    link_transcriptome(args$annotation, args$transcripts, args$index)

    message("Creating transcript database...")
    build_txdb(bfc, args$annotation)
}
main()