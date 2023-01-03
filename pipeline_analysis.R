rm(list = ls())

# import library
{
    library("patchwork")
    library("reshape2")
    library("ggplot2")
    library("dendextend")
    library("ggrepel")
    library("ggplot2")
    library("dplyr")
    library("clusterProfiler")
    library("org.Mm.eg.db")
    library("cowplot")
    library("enrichplot")
}

# load data
{
    df_TPM <- read.table("./02_result/Expression_TPM.tsv", header = T, sep = "\t", row.names = 1)
    df_phe <- read.table("./00_metadata/group.csv", header = T, sep = ",")

    # df_phe$Sample 转化为 factor
    df_phe$Sample <- factor(df_phe$Sample, levels = df_phe$Sample)
    df_phe$SRR_ID <- as.factor(df_phe$SRR_ID)
    df_phe$Group <- as.factor(df_phe$Group)
}

# PCA
{
    pca.input <- t(df_TPM) # transpose for PCA
    pca.input <- pca.input[, which(apply(pca.input, 2, var) != 0)] # filter out genes with zero variance
    pca <- prcomp( # PCA
        pca.input,
        center = TRUE,
        scale = TRUE
    )
    pca.data <- data.frame(pca$x) # PCA data
    pca.variance <- pca$sdev^2 / sum(pca$sdev^2) # PCA variance

    pcaPlotDat <- cbind(df_phe, pca$x[, 1:3]) # PCA plot data
    pcaPlotDat

    # plot PCA result
    p_PCA <- ggplot(pcaPlotDat) +
        theme_classic() +
        geom_point(aes(x = PC1, y = PC2, shape = Group, color = Group), size = 3) +
        xlab(paste("PC1 (", round(pca.variance[1] * 100, 1), "%)", sep = "")) +
        ylab(paste("PC2 (", round(pca.variance[2] * 100, 1), "%)", sep = "")) +
        # theme(legend.title = element_blank()) +
        scale_color_manual(values = c("#44af8f", "#d95f02"))
    # stat_ellipse(aes(x = PC1, y = PC2, color = Group),level = 0.95, show.legend = F) +

    p_PCA
    ggsave(p_PCA, filename = "./03_analysis/PCA.pdf", width = 4, height = 3)
}

# hierarchical clustering analysis.
{
    # plot Cluster Dendrogram
    tree_corD <- hclust(as.dist(1 - cor(t(pca$x[, 1:6]))))
    dend <- as.dendrogram(tree_corD)
    pdf("./03_analysis/hclust.pdf", 7.5, 6)
    par(mar = c(6, 6, 6, 6))
    dend %>%
        set("labels_col", value = c("#44af8f", "#d95f02"), k = 2) %>%
        plot(main = "")

    dev.off()
}

# Differential expression
{
    gene_exp <- read.table("./02_result/gene_exp.diff", header = T)
    gene_exp <- gene_exp[, c(2, 3, 10, 12, 13)]
    gene_exp <- gene_exp[is.finite(gene_exp$log2.fold_change.), ]

    # add sig column: Up, Down, Nosig
    gene_exp$Sig <- "Nosig"
    gene_exp$Sig[gene_exp$p_value < 0.05 & gene_exp$log2.fold_change. > 1] <- "Up"
    gene_exp$Sig[gene_exp$p_value < 0.05 & gene_exp$log2.fold_change. < -1] <- "Down"
    gene_exp$Sig <- factor(gene_exp$Sig, levels = c("Up", "Nosig", "Down"))

    table(gene_exp$Sig)


    # select DEG
    DEG <- gene_exp[gene_exp$Sig != "Nosig", ]
    DEG <- DEG[order(abs(DEG$log2.fold_change.), decreasing = T), ] # sort by log2FC
    write.table(DEG, file = "./03_analysis/DiffGenes_FC2.csv ", sep = ",", quote = F, row.names = F)

}


# plot barplot

{
    tab <- as.data.frame(table(DEG$Sig))
    tab <- tab[tab$Var1 != "Nosig", ]


    tab$Var1 <- factor(tab$Var1, levels = c("Up", "Down"))

    p <- ggplot(tab, aes(x = Var1, y = Freq, label = Freq, fill = Var1)) +
        geom_bar(stat = "identity", color = "black", size = 0.7) +
        geom_text(position = position_dodge(0.9), vjust = -0.2, size = 2.5) +
        ylim(0, max(tab$Freq) * 1.1) +
        theme_classic() +
        xlab("Differential expression") +
        ylab("Number of genes") +
        ggtitle("") +
        theme(legend.position = "none") +
        scale_fill_manual(values = c("#b2182b99", "#2166ac99"))

    p
    ggsave(p, filename = "./03_analysis/diffgene_number_barplot.pdf", width = 4, height = 4)
}

# Volcano plot
{
    color <- c(Up = "#b2182b", Nosig = "grey", Down = "#2166ac")

    p <- ggplot(gene_exp, aes(log2.fold_change., -log10(p_value), color = Sig)) +
        geom_point(alpha = 0.4) +
        theme_bw() +
        scale_color_manual(values = color) +
        # auxiliary line
        geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
        geom_vline(xintercept = c(-1, 1), lty = 4, col = "grey", lwd = 0.6) +
        xlim(-6, 6) +
        ylab(expression(-log[10](p ~ value))) +
        xlab(expression(log[2](Fold ~ Change))) +
        theme(
            legend.position = "right",
            panel.grid = element_blank(),
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 14)
        ) +
        theme(legend.title = element_blank())
    p

    # selected_genes
    selected_genes <- c(
        "Stat1", "Irf9", "Oas3", "Usp18", "Ddx58", "Isg15", "Ifi203", "Ifit2", "Ifi44", "Cd300If",
        "Nefh", "Rorb", "Nefm", "Cplx1", "Syt2", "Chrm2", "Kcnc1", "Kcnc3", "Kcnip2", "Shroom2", "Atg9a", "Chtf18", "Etnppl"
    )
    selected_genes <- gene_exp[gene_exp$gene %in% selected_genes, ]

    # plot volcano label
    p_topn <- p + geom_label_repel(
        data = selected_genes,
        aes(log2.fold_change., -log10(p_value), label = gene, fill = Sig),
        color = "white",
        size = 4, # font size
        box.padding = 0.8, # font to point distance
        segment.colour = "black",
        segment.size = 0.7, # line width
        min.segment.length = 0, # min line length
        show.legend = FALSE,
    ) + scale_fill_manual(values = c("#b2182b", "#2166ac"))
    p_topn
    ggsave(p_topn, filename = "./03_analysis/diffgene_volcano.pdf", width = 6, height = 6)
}

# Heatmap

{
    # selected genes
    selected_up <- DEG[DEG$Sig == "Up", ][1:40, ]
    selected_down <- DEG[DEG$Sig == "Down", ][1:20, ]
    selected_DEG <- rbind(selected_up, selected_down) -> selected_genes
    sigGeneDat <- df_TPM[selected_DEG$gene_id, ]
    rownames(sigGeneDat) <- selected_DEG$gene


    periodCols <- c("#44af8f", "#d95f02")
    names(periodCols) <- levels(df_phe$Group)
    periodColors <- periodCols[as.character(df_phe$Group)]

    rowHClust <- hclust(as.dist(1 - cor(t(sigGeneDat))))
    rowTree <- as.dendrogram(rowHClust)

    # cutree from hierarchical clustering
    gClust <- cutree(rowHClust, k = 2)


    pdf("./03_analysis/diffgene_heatmap.pdf", 6, 12)

    heatmap.2(as.matrix(sigGeneDat),
        Rowv = rowTree, dendrogram = "both", scale = "row", trace = "none", col = redgreen(75), ColSideColors = periodColors, margins = c(8, 10), lhei = c(1, 9),
        cexRow = 0.8, cexCol = 1.0 # fontszie
    )
    dev.off()
}


# GO&KEGG
{
    rownames(df_TPM) <- gsub("\\..*", "", rownames(df_TPM))
    rownames(DEG) <- gsub("\\..*", "", DEG$gene_id)
    DEG_up <- DEG[DEG$Sig == "Up", ]
    DEG_down <- DEG[DEG$Sig == "Down", ]

    # gene id translaste to ENTREZID
    DEG_up_gene <- bitr(rownames(DEG_up),
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Mm.eg.db
    )

    DEG_down_gene <- bitr(rownames(DEG_down),
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Mm.eg.db
    )

    geneList <- bitr(rownames(df_TPM),
        fromType = "ENSEMBL",
        toType = "ENTREZID",
        OrgDb = org.Mm.eg.db
    )


    # GO term  enrichment analysis
    DEG_up_go <- enrichGO(
        gene = DEG_up_gene$ENTREZID,
        universe = names(geneList$ENTREZID),
        OrgDb = org.Mm.eg.db,
        ont = "BP", # BP: biological process, MF: molecular function, CC: cellular component
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
    )

    DEG_up_go_read <- DOSE::setReadable(DEG_up_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENTREZID"
    ) # ENTREZID to gene Symbol
    write.csv(DEG_up_go_read@result, "./03_analysis/GO_DEG_up_enrichresults.csv")



    DEG_down_go <- enrichGO(
        gene = DEG_down_gene$ENTREZID,
        universe = names(geneList$ENTREZID),
        OrgDb = org.Mm.eg.db,
        ont = "BP", # BP: biological process, MF: molecular function, CC: cellular component
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
    )


    DEG_down_go_read <- DOSE::setReadable(DEG_down_go,
        OrgDb = "org.Mm.eg.db",
        keyType = "ENTREZID"
    ) # ENTREZID to gene Symbol
    write.csv(DEG_down_go_read@result, "./03_analysis/GO_DEG_down_enrichresults.csv")

    # dotplot
    p1 <- dotplot(DEG_up_go, showCategory = 10) + ggtitle("Downregulated")
    p2 <- dotplot(DEG_down_go, showCategory = 10) + ggtitle("Upregulated")
    pp <- plot_grid(p1, p2, ncol = 2, labels = LETTERS[1:2])
    ggsave(pp, filename = "./03_analysis/GO_dotplot.pdf", width = 12, height = 10)

    # cnetplot
    pp_go_cnetplot <- cnetplot(DEG_up_go,
        showCategory = 6, colorEdge = T,
        node_label = "all",
        color_category = "steelblue"
    )
    pp_down_cnetplot <- cnetplot(DEG_down_go,
        showCategory = 6, colorEdge = T,
        node_label = "all",
        color_category = "steelblue"
    )

    pp_cnetplot <- plot_grid(pp_go_cnetplot, pp_down_cnetplot, ncol = 2, labels = LETTERS[1:4])

    ggsave(pp_cnetplot, filename = "./03_analysis/GO_cnetplot.pdf", width = 24, height = 10)


    # emapplot
    DEG_up_go_pair <- pairwise_termsim(DEG_up_go)
    pp_go_emapplot <- emapplot(DEG_up_go_pair)
    DEG_down_go_pair <- pairwise_termsim(DEG_down_go)
    pp_down_emapplot <- emapplot(DEG_down_go_pair)

    pp_emapplot <- plot_grid(pp_go_emapplot, pp_down_emapplot, ncol = 2, labels = LETTERS[1:2])

    ggsave(pp_emapplot, filename = "./03_analysis/GO_emapplot.pdf", width = 24, height = 10)

    # goplot
    pp_go_goplot <- goplot(DEG_up_go)
    pp_down_goplot <- goplot(DEG_down_go)

    pp_goplot <- plot_grid(pp_go_goplot, pp_down_goplot, ncol = 2, labels = LETTERS[1:2])

    ggsave(pp_goplot, filename = "./03_analysis/GO_goplot.pdf", width = 24, height = 10)


    # KEGG pathway enrichment analysis
    DEG_up_kk <- enrichKEGG(
        gene = DEG_up_gene$ENTREZID,
        organism = "mmu",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    DEG_down_kk <- enrichKEGG(
        gene = DEG_down_gene$ENTREZID,
        organism = "mmu",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
    )

    p_dotplot_up_kk <- dotplot(DEG_up_kk, showCategory = 10) + ggtitle("Upregulated")
    p_dotplot_down_kk <- dotplot(DEG_down_kk, showCategory = 10) + ggtitle("Downregulated")
    pp2 <- plot_grid(p_dotplot_up_kk, p_dotplot_down_kk, ncol = 2, labels = LETTERS[1:2])
    ggsave(pp2, filename = "./03_analysis/KEGG_dotplot.pdf", width = 11, height = 7)
}