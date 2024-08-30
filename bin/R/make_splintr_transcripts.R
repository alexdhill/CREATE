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

make_splintr_transcripts <- function(genome, genes, version)
{
    suppressPackageStartupMessages({
        library(tidyr)
        library(stringr)
        library(Biostrings)
        library(GenomicFeatures)
        library(BSgenome)
        library(eisaR)
    })

    message("Reading in annotation...")
    grl = suppressMessages(eisaR::getFeatureRanges(
        gtf=file.path(genes),
        featureType=c("intron", "spliced"),
        intronType="separate",
        verbose=FALSE
    ))
    
    make_exons = function(x) {
        x$type = "exon"
        x$exon_rank = 1L
        x$gene_id = names(x)
        x$transcript_id = make.unique(x$gene_id, sep="")
        x$exon_id = x$transcript_id
        names(x) = NULL
        return(x)
    }

    message("Reformatting introns...")
    spliced_grl = grl[!grepl('-I[0-9]*', names(grl))]
    intron_gr = grl[grepl('-I[0-9]*', names(grl))] %>%
        BiocGenerics::unlist() %>%
        split(., .$gene_id) %>%
        reduce() %>%
        BiocGenerics::unlist() %>%
        make_exons()
    
    message("Reintroducing introns...")
    mcols(intron_gr) = S4Vectors::mcols(intron_gr)[,c("exon_id", "exon_rank", "transcript_id", "gene_id", "type")]
    intron_grl = BiocGenerics::relist(intron_gr, lapply(
        structure(seq_along(intron_gr), names=intron_gr$transcript_id),
        function(i) i
    ))
    splintr = c(spliced_grl, intron_grl)

    message("Reading genome")
    genome_seqs = Biostrings::readDNAStringSet(file.path(genome))
    names(genome_seqs) = word(names(genome_seqs), 1)

    message("Making splintr annotation")
    rm(intron_gr, spliced_grl, intron_grl, grl)
    seqlevels(splintr) = seqlevels(genome_seqs)
    seqlengths(splintr) = suppressWarnings(seqlengths(genome_seqs))
    splintr = trim(splintr)

    message("Making splintr transcripts")
    splintr_transcripts = GenomicFeatures::extractTranscriptSeqs(
        x=genome_seqs, transcripts=splintr
    )
    rm(genome_seqs)

    message("Writing transcripts")
    writeXStringSet(
        splintr_transcripts,
        file.path(paste0(version, "_splintr_transcripts.fa.gz")),
        format="fasta", compress=TRUE
    )

    message("Making transcript maps")
    splintr_tx2g = getTx2Gene(splintr)
    write.table(
        splintr_tx2g, file.path(paste0(version, "_splintr_map.tx2g")),
        sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
    )
    splintr_tx3g = splintr_tx2g %>%
        dplyr::mutate(
            gene_id=word(gene_id, 1, sep='-'),
            status=ifelse(grepl('-I[0-9]*$', transcript_id), 'U', 'S')
        )
    write.table(
        splintr_tx3g, file.path(paste0(version, "_splintr_map.tx3g")),
        sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE
    )
}

main <- function()
{
    suppressPackageStartupMessages({
        library(argparse)
    })

    parser = ArgumentParser()
    parser$add_argument(
        "-a", "--annotation", action="store",
        help="The Genode GTF annotation"
    )
    parser$add_argument(
        "-g", "--genome", action="store",
        help="The reference genome"
    )
    parser$add_argument(
        "-v", "--version", action="store",
        help="The gencode annotation version"
    )
    args = parser$parse_args()

    make_splintr_transcripts(args$genome, args$annotation, args$version)
}
main()