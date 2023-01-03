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
}