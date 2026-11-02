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

## Repeat loci in the complete annotation are named <repName>_range=<chr>:<start>-<end>_strand=<+/->
repeat_pattern = "_range=[^:]+:[0-9]+-[0-9]+_strand=[+-]$"

mask_introns = function(intron_gr, repeat_gr, stranded, min_width)
{
    ## Remove repeat sequence from intron targets so intronic repeat reads are
    ## not split between the repeat locus and its host gene's intron. With a
    ## stranded library only same-strand repeats compete with the intron.
    n_before = sum(BiocGenerics::width(intron_gr))
    gene_ids = names(intron_gr)
    masked = GenomicRanges::subtract(intron_gr, repeat_gr, ignore.strand=!stranded)
    names(masked) = gene_ids
    masked = BiocGenerics::unlist(masked)
    masked = masked[BiocGenerics::width(masked) >= min_width]
    message(sprintf(
        "Kept %.1f%% of intron sequence after masking (%d segments)",
        100*sum(BiocGenerics::width(masked))/n_before, length(masked)
    ))
    return(masked)
}

tag_hosts = function(repeat_gr, gene_exons)
{
    ## A repeat locus wholly inside a gene's span, on the gene's strand, and
    ## clear of that gene's exons is labelled <repName>@<gene_id>. Its UMIs are
    ## still counted once; the label lets the count be summed to the repeat or
    ## returned to the host gene's unspliced count. Nested genes resolve to the
    ## shortest containing host.
    exon_grl = GenomicRanges::reduce(split(gene_exons, gene_exons$gene_id))
    span_grl = range(exon_grl)
    span_gr = BiocGenerics::unlist(span_grl[lengths(span_grl) == 1])
    within = as.data.frame(GenomicRanges::findOverlaps(
        repeat_gr, span_gr, type="within", ignore.strand=FALSE
    ))
    within$host = names(span_gr)[within$subjectHits]
    within$width = BiocGenerics::width(span_gr)[within$subjectHits]
    exonic = as.data.frame(GenomicRanges::findOverlaps(
        repeat_gr, gene_exons, ignore.strand=FALSE
    ))
    exonic_keys = unique(paste(exonic$queryHits, gene_exons$gene_id[exonic$subjectHits]))
    hosts = within %>%
        dplyr::filter(!(paste(queryHits, host) %in% exonic_keys)) %>%
        dplyr::arrange(queryHits, width) %>%
        dplyr::distinct(queryHits, .keep_all=TRUE)
    message(sprintf("Host-tagged %d of %d repeat loci", nrow(hosts), length(repeat_gr)))
    tags = stats::setNames(rep(NA_character_, length(repeat_gr)), names(repeat_gr))
    tags[hosts$queryHits] = hosts$host
    return(tags)
}

make_splintr_transcripts <- function(genome, genes, version, mask="none", min_width=32L, repeats_only=FALSE)
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
    ## eisaR cannot build introns from an annotation without multi-exon
    ## transcripts, so a repeat-only build asks for spliced features only
    grl = tryCatch(
        suppressMessages(eisaR::getFeatureRanges(
            gtf=file.path(genes),
            featureType=if (repeats_only) "spliced" else c("intron", "spliced"),
            intronType="separate",
            joinOverlappingIntrons = TRUE,
            verbose=TRUE
        )),
        error=function(e) {
            if (grepl("in value to replace 0 elements", conditionMessage(e), fixed=TRUE)) {
                stop(paste0(
                    "No introns found in ", genes, ": the annotation has no multi-exon transcripts. ",
                    "For a repeat-aware build the gene annotation is missing or empty ",
                    "(check the *_gencode_annotation.gtf.gz input); ",
                    "for a repeat-only build pass --repeats_only."
                ), call.=FALSE)
            }
            stop(e)
        }
    )
    
    make_exons = function(x) {
        x$type = "exon"
        x$exon_rank = 1L
        x$gene_id = names(x)
        x$transcript_id = make.unique(x$gene_id, sep="")
        x$exon_id = x$transcript_id
        names(x) = NULL
        return(x)
    }

    ## anchored: RepeatMasker names such as MamGypsy2-I also contain '-I'
    is_intron = grepl('-I[0-9]*$', names(grl))
    spliced_grl = grl[!is_intron]
    splintr = spliced_grl
    host_tags = NULL
    if (any(is_intron))
    {
        message("Reformatting introns...")
        intron_gr = grl[is_intron] %>%
            BiocGenerics::unlist() %>%
            split(., .$gene_id) %>%
            reduce() %>%
            BiocGenerics::unlist()

        repeat_gr = BiocGenerics::unlist(spliced_grl[grepl(repeat_pattern, names(spliced_grl))])
        if (mask != "none" && length(repeat_gr) > 0)
        {
            message(paste0("Masking ", ifelse(mask=="stranded", "same-strand ", ""), "repeats out of introns..."))
            intron_gr = mask_introns(intron_gr, repeat_gr, mask=="stranded", min_width)
            message("Tagging intronic repeat loci with their host genes...")
            gene_exons = BiocGenerics::unlist(spliced_grl[!grepl(repeat_pattern, names(spliced_grl))])
            host_tags = tag_hosts(repeat_gr, gene_exons)
            rm(gene_exons)
        }
        rm(repeat_gr)
        intron_gr = make_exons(intron_gr)

        message("Reintroducing introns...")
        mcols(intron_gr) = S4Vectors::mcols(intron_gr)[,c("exon_id", "exon_rank", "transcript_id", "gene_id", "type")]
        intron_grl = BiocGenerics::relist(intron_gr, lapply(
            structure(seq_along(intron_gr), names=intron_gr$transcript_id),
            function(i) i
        ))
        splintr = c(spliced_grl, intron_grl)
        rm(intron_gr, intron_grl)
    } else
    {
        message("No intron features (repeat-only annotation): building spliced targets only")
    }

    message("Reading genome")
    genome_seqs = Biostrings::readDNAStringSet(file.path(genome))
    names(genome_seqs) = word(names(genome_seqs), 1)

    message("Making splintr annotation")
    rm(spliced_grl, grl)
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
            status=ifelse(grepl('-I[0-9]*$', transcript_id), 'U', 'S'),
            gene_id=ifelse(status=='U', sub('-I$', '', gene_id), gene_id)
        )
    if (!is.null(host_tags))
    {
        hosted = match(splintr_tx3g$transcript_id, names(host_tags))
        hosted = ifelse(is.na(hosted), NA_character_, host_tags[hosted])
        splintr_tx3g = splintr_tx3g %>%
            dplyr::mutate(gene_id=ifelse(is.na(hosted), gene_id, paste0(gene_id, "@", hosted)))
    }
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
    parser$add_argument(
        "-m", "--mask", action="store", default="none",
        choices=c("stranded", "unstranded", "none"),
        help="Remove repeat sequence from intron targets and host-tag intronic repeats"
    )
    parser$add_argument(
        "--repeats_only", action="store_true", default=FALSE,
        help="The annotation contains only repeats (no introns to extract)"
    )
    parser$add_argument(
        "--min_width", action="store", type="integer", default=32L,
        help="Minimum intron segment kept after masking (salmon drops targets <= k)"
    )
    args = parser$parse_args()

    make_splintr_transcripts(args$genome, args$annotation, args$version, args$mask, args$min_width, args$repeats_only)
}
main()