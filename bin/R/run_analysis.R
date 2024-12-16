main = function()
{
    suppressMessages({
        library(tidyr)
        library(dplyr)
        library(magrittr)
        library(stringr)
        library(readr)
        library(argparse)
        library(ggplot2)
        library(patchwork)
        library(HDF5Array)
        library(SummarizedExperiment)
    })

    parser = ArgumentParser()
    parser$add_argument("--quants", help = "The `CREATE quant` counts dir", type = "character")
    parser$add_argument("--helpers", help = "The Kim Lab R helpers file", type = "character")
    parser$add_argument("--output", help = "The output directory", type = "character")
    args = parser$parse_args()

    source(args$helpers)
    dir.create("data")
    dir.create("plots")
    createTheme()

    message("Importing quantifications...")
    quants = loadHDF5SummarizedExperiment(args$quants)

    message("Saving normalized counts...")
        quants %>%
            assays() %>%
            extract2("normalized") %>%
            as.data.frame() %>%
            mutate(gene_id=rownames(.)) %>%
            write_csv("data/normalized_counts.csv", col_names=T, progress=F)
    
    master = quants %>%
        assays() %>%
        extract2("normalized") %>%
        as.data.frame() %>%
        mutate(gene_id = rownames(.)) %>%
        left_join(as.data.frame(rowData(quants)), by="gene_id") %>%
        pivot_longer(
            cols = -c("gene_id", "gene_name", "gene_biotype", "tx_ids"),
            names_to = "sample",
            values_to = "count"
        ) %>%
        select(-c(tx_ids)) %>%
        left_join(as.data.frame(colData(quants)), by=c("sample"="names"))

    message("Plotting biotype detection...")
    detection_plot = master %>%
        filter(count >= 10) %>%
        group_by(gene_biotype, condition, sample) %>%
        summarise(ndetected=n()) %>%
        ggplot(aes(x=factor(gene_biotype, levels=names(biotype_colors)), y=ndetected, fill=factor(gene_biotype, levels=names(biotype_colors))))+
            facet_wrap(~condition, scales="free")+
            geom_boxplot(alpha=0.5)+
            labs(x="Condition", y="# Genes (>=10 transcripts)")+
            scale_y_continuous(trans="log10")+
            scale_fill_manual(values=biotype_colors, guide="none")
    ggsave(
        plot=detection_plot, file="plots/gene_detection.pdf",
        width=6, height=3, units="in", dpi=300
    )

    message("Plotting biotype composition...")
    composition_plot = master %>%
        group_by(sample, condition, gene_biotype) %>%
        summarise(count=sum(count)) %>%
        ggplot(aes(x=sample, y=count, fill=factor(gene_biotype, levels=names(biotype_colors))))+
            facet_grid(~condition, scales="free_x", space="free")+
            geom_bar(stat="identity")+
            labs(x="Sample ID", y="Normalized Counts", fill="Gene Biotype")+
            scale_fill_manual(values=biotype_colors)+
            theme(axis.text.x=element_text(angle=90, hjust=1))
    ggsave(
        plot=composition_plot, file="plots/biotype_composition.pdf",
        width=6, height=4, units="in", dpi=300
    )

    message("Plotting PCA of normalized counts...")
    pca_vectors = quants %>%
        assays() %>%
        extract2("normalized") %>%
        as.data.frame() %>%
        prcomp()
    pca_variance = summary(pca_vectors)$importance[2,]
    pca_plot = pca_vectors$rotation %>%
        as.data.frame() %>%
        mutate(sample=rownames(.)) %>%
        left_join(as.data.frame(colData(quants)), by=c("sample"="names")) %>%
        ggplot(aes(x=PC1, y=PC2, color=condition))+
            geom_point()+
            labs(
                x=paste0("PC1 (", round(pca_variance[1]*100, digits=2), "%)"),
                y=paste0("PC2 (", round(pca_variance[2]*100, digits=2), "%)"),
                color="Condition"
            )
    variance_plot = data.frame(importance=unname(pca_variance), comps=factor(names(pca_variance))) %>%
        slice_head(n=5) %>%
        ggplot(aes(x=comps, y=importance*100))+
            geom_bar(stat="identity")+
            labs(x="principal Components", y="Variance Explained (%)")
    ggsave(
        plot=pca_plot+variance_plot+plot_layout(guides="collect"),
        file="plots/pca.pdf",
        width=8, height=4, units="in", dpi=300
    )

    message("Running DESeq2...")
    deseq = runDESeq2(quants)
    if (length(is.na(deseq))>1 || !is.na(deseq))
    {
        conditions = unlist(lapply(unique(colData(quants)$condition), function(c) if (sum(colData(quants)$condition==c)>2) return(c))) # alpha
        if (length(conditions) >= 2)
        {
            message("Printing differential expression results...")
            comps = apply(combn(conditions, 2), 2, paste, collapse="_vs_")
            lapply(comps, function(comp)
            {
                message(paste("Running", comp, "comparison..."))
                deseq_res = results(deseq, contrast=c("condition", strsplit(comp, "_vs_")[[1]][1], strsplit(comp, "_vs_")[[1]][2])) %>%
                    as.data.frame() %>%
                    mutate(gene_id=rownames(.)) %>%
                    left_join(as.data.frame(rowData(quants)), by="gene_id") %>%
                    select(-c(tx_ids))
                write_csv(deseq_res, paste0("data/", comp, "_deseq_results.csv"), col_names=T, progress=F)

                volcano = deseq_res %>%
                    filter(!unlist(apply(., 1, anyNA))) %>%
                    mutate(status=case_when(
                        (padj < 0.05) & (abs(log2FoldChange) > 1) ~ "Sig",
                        TRUE ~ "Not"
                    )) %>%
                    mutate(gene_biotype=case_when(status=="Sig"~gene_biotype, TRUE~"Not")) %>%
                    ggplot(aes(x=log2FoldChange, y=-log10(padj), color=factor(gene_biotype, levels=c(names(biotype_colors), "Not"))))+
                        geom_point(aes(alpha=status), show.legend=TRUE)+
                        geom_text(
                            data=data.frame(
                                gene_id="", gene_biotype="Not",
                                log2FoldChange=min(deseq_res[!is.na(deseq_res$padj),]$log2FoldChange),
                                padj=min(deseq_res[!is.na(deseq_res$padj),]$padj),
                                lab=strsplit(comp, "_vs_")[[1]][2]
                            ),
                            aes(label=lab), vjust=-1, hjust=0
                        )+
                        geom_text(
                            data=data.frame(
                                gene_id="", gene_biotype="Not",
                                log2FoldChange=max(deseq_res[!is.na(deseq_res$padj),]$log2FoldChange),
                                padj=min(deseq_res[!is.na(deseq_res$padj),]$padj),
                                lab=strsplit(comp, "_vs_")[[1]][1]
                            ),
                            aes(label=lab), vjust=-1, hjust=1
                        )+
                        scale_color_manual(
                            breaks=names(biotype_colors),
                            limits=c(names(biotype_colors), "Not"),
                            values=c(biotype_colors, c("Not"='black'))
                        )+
                        scale_alpha_manual(guide="none", values=c("Sig"=1, "Not"=0.15))+
                        geom_hline(yintercept=-log10(0.05), linetype="dashed")+
                        geom_vline(xintercept=c(-1, 1), linetype="dashed")+
                        labs(x="log2(Fold Change)", y="-log10(adjusted p-value)", color="Gene Biotype")
                ggsave(
                    plot=volcano, file=paste0("plots/", comp, "_diff_expression.pdf"),
                    width=7, height=6, units="in", dpi=300
                )
            })
        }
        else
        {
            message("Less than 2 conditions detected, skipping differential expression...")
        }
    }
}
main()
