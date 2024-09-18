main = function()
{
    suppressMessages({
        library(argparse)
        library(ggplot2)
        library(patchwork)
    })

    parser = ArgumentParser()
    parser$add_argument("--quants", help = "The `CREATE quant` counts dir", type = "character")
    parser$add_argument("--helpers", help = "The Kim Lab R helpers file", type = "character")
    parser$add_argument("--samplesheet", help = "The sample,condition csv", type = "character")
    parser$add_argument("--output", help = "The output directory", type = "character")
    args = parser$parse_args()

    source(args$helpers)

    message("Importing quantifications...")
    quants = loadQuantsH5(args$quants)

    message("Reading sample sheet...")
    metadata = read_csv(args$samplesheet, col_names=F, progress=F, show_col_types=F) %>%
        set_colnames(c("sample", "condition"))

    message("Merging biotypes with quantifications...")
    raw_master = quants %>%
        assay() %>%
        as.data.frame() %>%
        mutate(gene_id=rownames(.)) %>%
        left_join(as.data.frame(rowData(quants)), by="gene_id") %>%
        pivot_longer(
            cols=-c("gene_id", "gene_name", "gene_biotype", "tx_ids"),
            names_to="sample",
            values_to="count"
        ) %>%
        select(-c(tx_ids)) %>%
        left_join(metadata, by="sample")

    dir.create("data")
    message("Saving raw counts...")
    quants %>%
        assay() %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        mutate(gene_id=rownames(.)) %>%
        write_csv("data/raw_counts.csv", col_names=T, progress=F)

    createTheme()

    message("Running DESeq2...")
    deseq = runDESeq2(quants, metadata)
    if (is.na(deseq))
    {
        message("Insufficient replicates for normalized counts...")
        message("Plotting with raw counts...")
        norm_master = raw_master
    }
    else
    {
        message("Saving normalized counts...")
        counts(deseq, normalized=T) %>%
            as.data.frame() %>%
            mutate(gene_id=rownames(.)) %>%
            write_csv("data/normalized_counts.csv", col_names=T, progress=F)
        
        norm_master = deseq %>%
            counts(normalized=T) %>%
            as.data.frame() %>%
            mutate(gene_id=rownames(.)) %>%
            left_join(as.data.frame(rowData(quants)), by="gene_id") %>%
            pivot_longer(
                cols=-c("gene_id", "gene_name", "gene_biotype", "tx_ids"),
                names_to="sample",
                values_to="count"
            ) %>%
            select(-c(tx_ids)) %>%
            left_join(metadata, by="sample")
    }

    dir.create("plots")

    message("Plotting gene detection by biotype...")
    detection = norm_master %>%
        filter(count >= 1) %>%
        group_by(gene_biotype, sample, condition) %>%
        summarise(ndetected=n()) %>%
        ggplot(aes(x=gene_biotype, y=ndetected, fill=factor(gene_biotype, levels=names(biotype_colors)), linetype=condition))+
            geom_boxplot(alpha=0.5)+
            labs(x="Gene Biotype", y="# Detected Genes", linetype="Sample Condition")+
            scale_y_continuous(trans="log10")+
            scale_fill_manual(values=biotype_colors, guide="none")
    ggsave(
        plot=detection, file="plots/gene_detection.pdf",
        width=8, height=4, units="in", dpi=300
    )

    message("Plotting biotype expression and composition by sample...")
    expression = norm_master %>%
        group_by(sample, condition, gene_biotype) %>%
        summarise(count=sum(count)) %>%
        ggplot(aes(x=sample, y=count, fill=factor(gene_biotype, levels=names(biotype_colors))))+
            facet_grid(~condition, scales="free_x", space="free")+
            geom_bar(stat="identity")+
            scale_fill_manual(values=biotype_colors)+
            labs(x="Sample ID", y="Normalized Counts", fill="Gene Biotype")+
            theme(axis.text.x=element_text(angle=90, hjust=1))
    ggsave(
        plot=expression, file="plots/biotype_expression.pdf",
        width=9, height=6, units="in", dpi=300
    )

    composition = norm_master %>%
        group_by(sample, condition, gene_biotype) %>%
        summarise(count=sum(count)) %>%
        group_by(sample, condition) %>%
        mutate(count=count/sum(count)) %>%
        ggplot(aes(x=sample, y=count, fill=factor(gene_biotype, levels=names(biotype_colors))))+
            geom_bar(stat="identity", position="stack")+
            facet_grid(~condition, scales='free_x', space="free")+
            scale_fill_manual(values=biotype_colors)+
            scale_y_continuous(lim=c(0,1.1), breaks=c(0,0.25,0.5,0.75,1))+
            geom_text(
                data = norm_master %>%
                    group_by(sample, condition) %>%
                    summarise(nreads=sum(count)) %>%
                    mutate(
                        label=paste0(round(nreads/1e6, 2), "M"),
                        count=1.01,
                        gene_biotype="Coding"
                    ),
                aes(label=label), color="black",
                angle=90, hjust=0, vjust=0.5, size=3
            )+
            geom_segment(aes(x=-Inf, y=-Inf, xend=-Inf, yend=1), linewidth=0.75)+
            labs(x="Sample", y="Fraction of Reads", fill="Gene Biotype")+
            theme(
                axis.text.x=element_text(angle=90, hjust=1),
                axis.line.y=element_blank(),
                legend.justification=0.75,
            )
    ggsave(
        plot=composition, file="plots/biotype_composition.pdf",
        width=12, height=6, units="in", dpi=300
    )

    message("Plotting PCA...")
    pca_res = deseq %>%
        counts(normalized=T) %>%
        as.data.frame() %>%
        prcomp()
    pca_vars = summary(pca_res)$importance[2,]
    pca_plot = pca_res$rotation %>%
        as.data.frame() %>%
        mutate(sample=rownames(.)) %>%
        left_join(metadata, by="sample") %>%
        ggplot(aes(x=PC1, y=PC2, color=condition))+
            geom_point()+
            labs(
                x=paste0("PC1 (", round(pca_vars[1]*100, digits=2), "%)"),
                y=paste0("PC2 (", round(pca_vars[2]*100, digits=2), "%)"),
                color="Condition"
            )
    var_plot = data.frame(importance=unname(pca_vars), comps=factor(names(pca_vars))) %>%
        slice_head(n=5) %>%
        ggplot(aes(x=comps, y=importance*100))+
            geom_bar(stat="identity")+
            labs(x="Principal Components", y="Variance Explained (%)")
    ggsave(
        plot=pca_plot+var_plot+plot_layout(guides="collect"),
        file="tmp_pc_analysis.pdf",
        width=12, height=6, units="in", dpi=300
    )

    message("Printing differential expression results...")
    conditions = unlist(lapply(unique(metadata$condition), function(c) if (sum(metadata$condition==c)>2) return(c))) # alpha
    if (length(conditions) > 2)
    {
        message("More than two conditions detected, running pairwise comparisons...")
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
    else if (length(conditions) == 2)
    {
        deseq_res = results(deseq, contrast=c("condition", conditions[1], conditions[2])) %>%
            as.data.frame() %>%
            mutate(gene_id=rownames(.)) %>%
            left_join(as.data.frame(rowData(quants)), by="gene_id") %>%
            select(-c(tx_ids))
        write_csv(deseq_res, "data/deseq_results.csv", col_names=T, progress=F)

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
                        lab=conditions[2]
                    ),
                    aes(label=lab), vjust=-1, hjust=0
                )+
                geom_text(
                    data=data.frame(
                        gene_id="", gene_biotype="Not",
                        log2FoldChange=max(deseq_res[!is.na(deseq_res$padj),]$log2FoldChange),
                        padj=min(deseq_res[!is.na(deseq_res$padj),]$padj),
                        lab=conditions[1]
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
            plot=volcano, file="plots/diff_expression.pdf",
            width=7, height=6, units="in", dpi=300
        )
    }
}
main()