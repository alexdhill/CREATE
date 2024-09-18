suppressMessages({
    library(tidyr)
    library(dplyr)
    library(magrittr)
    library(stringr)
    library(readr)
})

loadQuantsH5 <- function(quants)
{
    suppressMessages({
        library(HDF5Array)
        library(SummarizedExperiment)
        library(SingleCellExperiment)
    })

    loadHDF5SummarizedExperiment(quants) %>%
        return()
}

readTranscriptMap <- function(reference_dir, class=TRUE)
{
    map = reference_dir %>% 
        list.files(full.names=TRUE) %>%
        str_subset("complete_map.tx2g") %>%
        read_csv(col_names=TRUE, progress=FALSE, show_col_types=FALSE) %>%
        select(gene_id, gene_name, gene_biotype) %>%
        distinct() %>%
        mutate(
            gene_biotype=case_when(
                !str_detect(gene_biotype, ",") ~ gene_biotype,
                TRUE ~ unlist(lapply(gene_biotype, function(x){strsplit(x, ",")[[1]][ifelse(class, 2, 1)]}))
            )
        )%>%
        mutate(
            gene_biotype=case_when(
                startsWith(gene_biotype, "Mt-") | startsWith(gene_name, "Mt-") ~ "Mitochondrial",
                gene_biotype == "protein_coding" ~ "Coding",
                gene_biotype %in% c("lncRNA", "LINE", "SINE", "LTR", "DNA") ~ gene_biotype,
                gene_biotype == "Simple_repeat" ~ "Microsatellite",
                TRUE ~ "Other"
            )
        ) %>%
        return()
}

runDESeq2 <- function(quants, metadata)
{
    suppressMessages({
        library(DESeq2)
    })

    valid_conditions = metadata$condition %>%
        unique() %>%
        lapply(function(c){if (sum(metadata$condition==c)>2) return(c)}) %>%
        unlist()

    colData = colData(quants) %>%
        as.data.frame() %>%
        mutate(sample=names, names=NULL) %>%
        left_join(metadata, by="sample") %>%
        filter(condition %in% valid_conditions)

    deseq_data = quants %>%
        assay() %>%
        as.data.frame() %>%
        add(0.5) %>%
        floor() %>%
        select(all_of(sample))
    
    if (ncol(deseq_data) == 0) return(NA)

    deseq = DESeq2::DESeqDataSetFromMatrix(
            countData=deseq_data,
            colData=colData,
            design=~condition
        ) %>%
        DESeq2::DESeq()
    
    return(deseq)
}

createTheme <- function()
{
    suppressMessages({
        library(ggplot2)
        library(ggthemes)
        library(ggsci)
    })
    base_size = 10
    ggplot2::theme_set(
        ggthemes::theme_foundation(
            base_size = base_size,
            base_family = 'Helvetica'
        )+
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = base_size,
                face = 'bold'
            ),
            plot.background = ggplot2::element_rect(colour = NA),
            plot.tag = ggplot2::element_text(
                size = base_size,
                face = "bold"
            ),
            plot.margin = ggplot2::margin(0.04, 0.04, 0.04, 0.04, unit = "in"),
            panel.background = ggplot2::element_rect(colour = NA),
            panel.border = ggplot2::element_rect(colour = NA),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(),
            axis.line.x = NULL,
            axis.line.y = NULL,
            axis.text = ggplot2::element_text(size = rel(0.95)),
            axis.text.x = ggplot2::element_text(
                margin = ggplot2::margin(t = 0.8 * base_size/4),
                vjust = 1
            ),
            axis.text.x.top = ggplot2::element_text(
                margin = ggplot2::margin(b = 0.8 * base_size/4),
                vjust = 0
            ),
            axis.text.y = ggplot2::element_text(
                margin = ggplot2::margin(r = 0.5 * base_size/4),
                hjust = 1
            ),
            axis.text.y.right = ggplot2::element_text(
                margin = ggplot2::margin(l = 0.5 * base_size/4),
                hjust = 0
            ),
            axis.ticks = ggplot2::element_line(), 
            axis.ticks.length = ggplot2::unit(base_size/2.5, "pt"),
            axis.ticks.length.x = NULL,
            axis.ticks.length.x.top = NULL,
            axis.ticks.length.x.bottom = NULL,
            axis.ticks.length.y = NULL,
            axis.ticks.length.y.left = NULL,
            axis.ticks.length.y.right = NULL,
            strip.text = ggplot2::element_text(size = rel(0.8), face = "bold"),
            strip.background = ggplot2::element_blank(),
            strip.placement = "outside",
            legend.key.size= ggplot2::unit(0.25, "in"),
            legend.spacing = ggplot2::unit(0, "in"),
            legend.key = ggplot2::element_rect(colour = NA),
            legend.title = ggplot2::element_text(face="italic"),
            legend.text = ggplot2::element_text(face = "bold"),
            legend.justification = c(0, 0.75),
            legend.box.just = "right",
            legend.margin = ggplot2::margin(6, 6, 6, 6),
            legend.box.spacing = ggplot2::unit(-0.02, "in")
        )
    )
}

biotype_colors = list(
    "Coding" = "#404040",
    "lncRNA" = "#00DAE6",
    "SINE" = "#BD4F00",
    "LTR" = "#88D100",
    "LINE" = "#6A0094",
    "DNA" = "#010A80",
    "Microsatellite" = "#A80202",
    "Mitochondrial" = "#006B1E",
    "Other" = "#574400"
)