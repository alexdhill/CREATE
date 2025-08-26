main = function()
{
    suppressMessages({
        library(tidyr)
        library(dplyr)
        library(purrr)
        library(magrittr)
        library(stringr)
        library(readr)
        library(argparse)
    })

    parser = ArgumentParser()
    parser$add_argument("--gtf", help = "The repeat annotation", type = "character")
    parser$add_argument("--info", help = "The RepeatMasker info file", type = "character")
    parser$add_argument("--list", help = "The list of repeats in the ann file", type = "character")
    parser$add_argument("--out", help = "The output file", type = "character")
    args = parser$parse_args()

    message("Importing repeat annotation...")
    ann = read_tsv(args$gtf, show_col_types = FALSE, col_names = FALSE) %>%
        select(1:9) %>%
        set_colnames(c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")) %>%
        mutate(feature="gene")

    message("Importing repeat list...")
    list = read_tsv(args$list, show_col_types = FALSE, col_names = FALSE) %>%
        set_colnames(c("family"))

    message("Importing RepeatMasker info...")
    info = read_tsv(args$info, show_col_types = FALSE, col_names = FALSE) %>%
        set_colnames(c("family", "subclass", "super")) %>%
        distinct() %>%
        select(family, subclass)

    message("Adding subclass biotypes...")
    biotypes = list %>%
        left_join(info, by="family", multiple="any")
    ann$attribute = paste0(ann$attribute, " gene_type \"", biotypes$subclass, "\";")
    
    message("Structuring and writing annotation...")
    gtf = ann %>%
        mutate(count=3) %>%
        uncount(count) %>%
        mutate(
            feature=rep(c("gene", "transcript", "exon"), n()/3),
            attribute=case_when(
                feature=="exon" ~ paste0(attribute, " exon_number 1;"),
                TRUE ~ attribute
            )
        ) %>%
        as.data.frame()
    head(gtf)
    write_tsv(gtf, args$out, col_names=FALSE)
}
main()