#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

gene_expression <- as.matrix(read.csv(args[1], header=FALSE, sep=",")[-1,-1])
rownames(gene_expression) <- read.csv(args[1], header=FALSE, sep=",")[-1,1]
colnames(gene_expression) <- read.csv(args[1], header=FALSE, sep=",")[1,-1]

png(args[2],
    width = 1500,
    height = 1500,
    res = 300,
    pointsize = 8)

heatmap(gene_expression)
dev.off()