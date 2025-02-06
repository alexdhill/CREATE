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

link_transcriptome <- function(annotation, transcripts, index) {
    build <- strsplit(annotation, "_")[[1]][1]
    out_json <- paste0(build, "_complete_txome.json")
    message("Making linked txome...")
    makeLinkedTxome(
        indexDir = index,
        fasta = transcripts,
        gtf = annotation,
        source = "CREATE",
        organism = "Homo sapiens",
        genome = ifelse(build == "T2Tv2", "CHM13", "GRCh38"),
        release = build,
        jsonFile = out_json,
        write = TRUE
    )
}

build_txdb <- function(bfc, annotation) {
    genome=strsplit(annotation, "_")[[1]][1]
    if (genome=="T2Tv2")
    {
        message("Parsing T2T GFF3 attributes...")
        .GENE_TYPES <- c(
            "gene",
            "protein_coding_gene",
                "gene_with_edited_transcript",
                "gene_with_mRNA_with_frameshift",
                "gene_with_polyadenylated_mRNA",
                "gene_with_recoded_mRNA",
                    "gene_with_mRNA_recoded_by_translational_bypass",
                    "gene_with_transcript_with_translational_frameshift",
                    "gene_with_stop_codon_read_through",
                        "gene_with_stop_codon_redefined_as_pyrrolysine",
                        "gene_with_stop_codon_redefined_as_selenocysteine",
            "ncRNA_gene",
                "snoRNA_gene",
                "gRNA_gene",
                "RNase_P_RNA_gene",
                "tmRNA_gene",
                "SRP_RNA_gene",
                "tRNA_gene",
                "rRNA_gene",
                    "rRNA_25S_gene",
                    "rRNA_16S_gene",
                    "rRNA_21S_gene",
                    "rRNA_5S_gene",
                    "rRNA_28S_gene",
                    "rRNA_23S_gene",
                    "rRNA_5_8S_gene",
                    "rRNA_18S_gene",
                "lncRNA_gene",
                    "bidirectional_promoter_lncRNA",
                    "sense_overlap_ncRNA_gene",
                    "lincRNA_gene",
                    "antisense_lncRNA_gene",
                    "sense_intronic_ncRNA_gene",
                "enzymatic_RNA_gene",
                    "ribozyme_gene",
                "RNase_MRP_RNA_gene",
                "snRNA_gene",
                "telomerase_RNA_gene",
                "miRNA_gene",
                "piRNA_gene",
                "scRNA_gene",
            "pseudogene",
                "non_processed_pseudogene",
                    "duplicated_pseudogene",
                    "cassette_pseudogene",
                    "nuclear_mt_pseudogene",
                    "translated_unprocessed_pseudogene",
                    "pseudogene_by_unequal_crossing_over",
                    "transcribed_unprocessed_pseudogene",
                    "unitary_pseudogene",
                        "transcribed_unitary_pseudogene",
                        "allelic_pseudogene",
                "processed_pseudogene",
                    "transcribed_processed_pseudogene",
                    "translated_processed_pseudogene",
                "polymorphic_pseudogene",
                    "polymorphic_pseudogene_with_retained_intron",
                "vertebrate_immune_system_pseudogene",
                    "T_cell_receptor_pseudogene",
                        "TR_J_pseudogene",
                        "TR_V_pseudogene",
                    "immunoglobulin_pseudogene",
                        "IG_J_pseudogene",
                        "IG_V_pseudogene",
                        "IG_C_pseudogene",
                "transposable_element_pseudogene",
            "transposable_element_gene",
                "engineered_foreign_transposable_element_gene"
        )
        .TX_TYPES <- c(
            "transcript",
            "primary_transcript",
            "mRNA",
            "rRNA",
            "snoRNA",
            "snRNA",
            "tRNA",
            "tmRNA",
            "miRNA",
            "miRNA_primary_transcript",
            "RNase_P_RNA",
            "RNase_MRP_RNA",
            "SRP_RNA",
            "misc_RNA",
            "antisense_RNA",
            "antisense",
            "lnc_RNA",
            "antisense_lncRNA",
            "transcript_region",
            "scRNA",
            "guide_RNA",
            "telomerase_RNA",
            "vault_RNA",
            "Y_RNA",
            "ncRNA",
                "antisense_RNA",
                "three_prime_overlapping_ncrna",
                "vault_RNA",
            "predicted_transcript",
                "unconfirmed_transcript",
            "pseudogenic_transcript",
                "pseudogenic_rRNA",
                    "unitary_pseudogenic_rRNA",
                    "allelic_pseudogenic_rRNA",
                    "unprocessed_pseudogenic_rRNA",
                    "processed_pseudogenic_rRNA",
                "pseudogenic_tRNA",
                    "unitary_pseudogenic_tRNA",
                    "allelic_pseudogenic_tRNA",
                    "processed_pseudogenic_tRNA",
                    "unprocessed_pseudogenic_tRNA",
            "vertebrate_immunoglobulin_T_cell_receptor_segment",
                "V_gene_segment",
                "D_gene_segment",
                "C_gene_segment",
                "J_gene_segment",
            "gene_segment",
                "pseudogenic_gene_segment"
        )
        .EXON_TYPES <- c(
            "exon",
            "pseudogenic_exon",
            "coding_exon",
            "five_prime_coding_exon",
            "interior_coding_exon",
            "three_prime_coding_exon",
            "exon_of_single_exon_gene",
            "interior_exon",
            "noncoding_exon",
            "five_prime_noncoding_exon",
            "three_prime_noncoding_exon"
        )
        .CDS_TYPES <- c(
            "CDS",
            "transposable_element_CDS",
            "CDS_predicted",
            "edited_CDS"
        )
        .STOP_CODON_TYPES <- "stop_codon"
        GFF_FEATURE_TYPES <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES, .CDS_TYPES, .STOP_CODON_TYPES)

        granges = rtracklayer::import(
            annotation, format="gtf",
            feature.type=GFF_FEATURE_TYPES,
            colnames=c(
                "transcript_id", "gene_id", "gene_name", "transcript_biotype", "gene_biotype",
                "Name", "Parent", "Dbxref", "new_transcript_id", "new_gene_id", "type", "ID"
            )
        )
        txdb <- txdbmaker::makeTxDbFromGRanges(granges)
    }
    else
    {
        message("Parsing HG38 GTF attributes...")
        txdb <- txdbmaker::makeTxDbFromGFF(annotation, format="gtf", dataSource="CREATE", organism="Homo sapiens")
    }
    saveDb(txdb, bfcnew(bfc, rname = annotation, ext = ".sqlite"))
} 

main <- function() {
    parser <- ArgumentParser()
    parser$add_argument(
        "-a", "--annotation",
        action = "store",
        help = "The reference annotation"
    )
    parser$add_argument(
        "-t", "--transcripts",
        action = "store",
        help = "The reference transcriptome"
    )
    parser$add_argument(
        "-i", "--index",
        action = "store",
        help = "The salmon index/digest"
    )

    args <- parser$parse_args()

    message("Making file cache...")
    dir.create(".bfccache")
    bfc <- BiocFileCache(".bfccache")
    setTximetaBFC(".bfccache")

    link_transcriptome(args$annotation, args$transcripts, args$index)

    message("Creating transcript database...")
    build_txdb(bfc, args$annotation)
}
main()
