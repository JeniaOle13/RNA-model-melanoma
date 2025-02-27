setwd("/home/acari/Documents/github/probiotics-melanoma/")

# Import libraries
library(vegan)
library(fgsea)
library(limma)
library(DESeq2)
library(lmerTest)
library(permutes)
library(tidyverse)
library(immunedeconv)
library(clusterProfiler)

library(ggpubr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnhancedVolcano)

library(msigdbr)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)

mypal <- brewer.pal(8, "Set1")

set.seed(10)

# Import data
meta_data <- read.csv("data/metadata.tsv", sep = "\t")
meta_data$group <- as.factor(meta_data$group)
meta_data$batch <- as.factor(meta_data$batch)
meta_data$File <- meta_data$sampleid
meta_data <- meta_data[c(1,4,2,3)]
meta_data$File <- paste0(meta_data$File, ".counts")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data,
                                       directory = "data/htseq/",
                                       design= ~ 0 + group + batch)
rownames(ddsHTSeq) <- sapply(str_split(rownames(ddsHTSeq), "\\."), function(x) x[1])

# Filter gene table by counts
keep <- rowSums(counts(ddsHTSeq)>= 10) > ncol(ddsHTSeq)*0.3
ddsHTSeq <- ddsHTSeq[keep,]

# Diff-exp analysis
ddsHTSeq$group <- relevel(ddsHTSeq$group, ref = "M")
ddsHTSeq <- DESeq(ddsHTSeq)

COU <- counts(ddsHTSeq, normalized=TRUE)
# write.table(COU, "results/normalized.counts.tsv", sep = "\t", quote = F, row.names = T)

contr_1 <- makeContrasts(groupM - groupM_BIF, levels = resultsNames(ddsHTSeq))
contr_2 <- makeContrasts(groupM - groupM_LAC, levels = resultsNames(ddsHTSeq))
contr_3 <- makeContrasts(groupM_BIF - groupM_LAC, levels = resultsNames(ddsHTSeq))

res_1 <- results(ddsHTSeq, contrast=contr_1)
res_2 <- results(ddsHTSeq, contrast=contr_2)
res_3 <- results(ddsHTSeq, contrast=contr_3)

# Make MDS visualization
rld <- as.data.frame(t(assay(rlog(ddsHTSeq))))

mds <- metaMDS(rld, distance = "euclidean")
mds.points <- as.data.frame(mds$points)
mds.points <- merge(meta_data, cbind(rownames(mds.points), mds.points), by = 1)[-2]

nmds.plot <- ggplot(mds.points, aes(MDS1, MDS2, col = group, shape = batch))+
    geom_point(size = 2.75)+
    # geom_text(size = 3)+
    theme_classic()+
    theme(legend.position = "bottom")+
    scale_color_brewer(palette = "Set1")

# ggsave(filename = "figures/nmds.plot.pdf", plot = nmds.plot, device = "pdf", width = 5, height = 5)

# Make DEG tables
get_DEG <- function(res){
    DEG <- cbind(as.data.frame(res), COU)
    DEG <- DEG[!is.na(DEG$padj),]
    DEG <- DEG[order(DEG$log2FoldChange, decreasing = T),]
    return(DEG)
}

DEG_BIF <- get_DEG(res_1)
DEG_LAC <- get_DEG(res_2)
DEG_BIF_LAC <- get_DEG(res_3)

# write.table(cbind(ensembl = rownames(DEG_BIF), DEG_BIF), "results/DEG_BIF.tsv", sep = "\t", quote = F, row.names = F)
# write.table(cbind(ensembl = rownames(DEG_LAC), DEG_LAC), "results/DEG_LAC.tsv", sep = "\t", quote = F, row.names = F)
# write.table(cbind(ensembl = rownames(DEG_BIF_LAC), DEG_BIF_LAC), "results/DEG_BIF_LAC.tsv", sep = "\t", quote = F, row.names = F)

# Get volcano plots
VolcanoPlot1 <- EnhancedVolcano(DEG_BIF,
                lab = rownames(DEG_BIF),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M vs M_BIF", labSize = 0)

VolcanoPlot2 <- EnhancedVolcano(DEG_LAC,
                lab = rownames(DEG_LAC),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M vs M_LAC", 
                labSize = 0)

VolcanoPlot3 <- EnhancedVolcano(DEG_BIF_LAC,
                lab = rownames(DEG_BIF_LAC),
                x = 'log2FoldChange',
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M_BIF vs M_LAC", 
                labSize = 0)

volcano.all <- ggarrange(VolcanoPlot1, VolcanoPlot2, VolcanoPlot3, nrow = 1, common.legend = T)
# ggsave(filename = "figures/volcano.pdf", plot = volcano.all, device = "pdf", width = 16, height = 6)

# MSigDB GSEA
## Hallmark
pathwaysDF <- msigdbr("mouse", category="H")
pathways <- split(as.character(pathwaysDF$ensembl_gene), pathwaysDF$gs_name)

get_ranks <- function(DEG){
    ranks.df <- DEG$log2FoldChange
    names(ranks.df) <- rownames(DEG)
    names(ranks.df) <- sapply(str_split(names(ranks.df), "\\."), function(x) x[1])    
    ranks.df <- sort(ranks.df, decreasing = T)
    return(ranks.df)
}

ranks.bif <- get_ranks(DEG_BIF)
ranks.lac <- get_ranks(DEG_LAC)
ranks.bif.lac <- get_ranks(DEG_BIF_LAC)

fgsea.bif <- fgsea(pathways = pathways, 
                  stats    = ranks.bif,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

fgsea.lac <- fgsea(pathways = pathways, 
                     stats    = ranks.lac,
                     eps      = 0.0,
                     minSize  = 15,
                     maxSize  = 500)

fgsea.bif.lac <- fgsea(pathways = pathways, 
                   stats    = ranks.bif.lac,
                   eps      = 0.0,
                   minSize  = 15,
                   maxSize  = 500)

filter_pval <- function(fgsea){
    fgsea <- as.data.frame(fgsea)
    fgsea <- fgsea[fgsea$padj < 0.05,]
    fgsea <- fgsea[order(fgsea$NES, decreasing = T),]
    return(fgsea)
}

fgsea.bif.filter <- filter_pval(fgsea.bif)
fgsea.lac.filter <- filter_pval(fgsea.lac)
fgsea.bif.lac.filter <- filter_pval(fgsea.bif.lac)

# write.table(fgsea.bif.filter[-8], "results/fgsea.bif.tsv", sep = "\t", quote = F, row.names = F)
# write.table(fgsea.lac.filter[-8], "results/fgsea.lac.tsv", sep = "\t", quote = F, row.names = F)
# write.table(fgsea.bif.lac.filter[-8], "results/fgsea.bif.lac.tsv", sep = "\t", quote = F, row.names = F)

get_subset <- function(df, group){
    df.sbs <- df[c(1,6)]
    df.sbs$group <- group
    return(df.sbs)
}

fgsea.bif.sbs <- get_subset(fgsea.bif.filter, "M vs M_BIF")
fgsea.lac.sbs <- get_subset(fgsea.lac.filter, "M vs M_LAC")
fgsea.bif.lac.sbs <- get_subset(fgsea.bif.lac.filter, "M_BIF vs M_LAC")

fgsea.all <- rbind(fgsea.bif.sbs, fgsea.lac.sbs, fgsea.bif.lac.sbs)
fgsea.all <- spread(fgsea.all, group, NES, fill = 0)
rownames(fgsea.all) <- fgsea.all$pathway
fgsea.all <- fgsea.all[-1]

heatmap.hallmark <- pheatmap::pheatmap(as.matrix(fgsea.all), cutree_rows = 7, cutree_cols = 2, display_numbers = T)

# pdf(file = "figures/heatmap.hallmark.pdf", width = 6.5, height = 8.5)
# heatmap.hallmark
# dev.off()

# Immune cell types GSEA (Cell Marker)
cell_marker_data <- read.csv("data/cell_marker_mouse.csv",sep = "\t")
cell_names <- unique(cell_marker_data$cell_name)

cells <- cell_marker_data %>%
    dplyr::select(cell_name, GeneID)

t_cells <- cell_names[!is.na(str_extract(cell_names, "T cell"))]
b_cells <- cell_names[!is.na(str_extract(cell_names, "B cell"))]
neutrophil <- cell_names[!is.na(str_extract(cell_names, "Neutrophil|neutrophil"))]
dendritic_cell <- cell_names[!is.na(str_extract(cell_names, "Dendritic|dendritic"))]
basophil <- cell_names[!is.na(str_extract(cell_names, "Basophil|basophil"))]
mast_cell <- cell_names[!is.na(str_extract(cell_names, "Mast|mast"))]
macrophage <- cell_names[!is.na(str_extract(cell_names, "Macrophage|macrophage"))]
# histiocyte <- cell_names[!is.na(str_extract(cell_names, "Histiocyte|histiocyte"))]
kupffer_cell <- cell_names[!is.na(str_extract(cell_names, "Kupffer|kupffer"))]
plasma_cell <- cell_names[!is.na(str_extract(cell_names, "plasma|Plasma"))]
innate_lymphoid <- cell_names[!is.na(str_extract(cell_names, "Innate|Lymphoid|innate|lymphoid"))]
NK <- cell_names[!is.na(str_extract(cell_names, "NK|natural killer"))]

lymphocyte <- cell_names[!is.na(str_extract(cell_names, "Lymphocyte|lymphocyte"))]
monocyte <- cell_names[!is.na(str_extract(cell_names, "Monocyte|monocyte"))]
granulocyte <- cell_names[!is.na(str_extract(cell_names, "Granulocyte|granulocyte"))]

immune.cells <- unique(c(t_cells, b_cells, neutrophil, dendritic_cell, basophil, mast_cell, 
         macrophage, kupffer_cell, plasma_cell, innate_lymphoid, NK, 
         lymphocyte, monocyte, granulocyte))

ens2ent <- mapIdscell_nameens2ent <- mapIds(org.Mm.eg.db, keys = rownames(COU), keytype="ENSEMBL", column = "ENTREZID")
ens2ent <- as.data.frame(ens2ent)
ens2ent <- cbind(ENSEMBL = rownames(ens2ent), ENTREZID = ens2ent$ens2ent)
row.names(ens2ent) <- 1:nrow(ens2ent)

get_ENTREZID_rank <- function(rank.df, ens2ent){
    rank.df <- merge(ens2ent, cbind(names(rank.df), rank.df), by = 1)
    rank.df <- rank.df[!is.na(rank.df$ENTREZID),]
    rank.df <- rank.df[!duplicated(rank.df$ENTREZID),]
    
    rank.ent <- as.numeric(rank.df$rank.df)
    names(rank.ent) <- rank.df$ENTREZID
    rank.ent <- sort(rank.ent, decreasing = T)
    return(rank.ent)
}

ranks.bif.ent <- get_ENTREZID_rank(ranks.bif, ens2ent)
ranks.lac.ent <- get_ENTREZID_rank(ranks.lac, ens2ent)
ranks.bif.lac.ent <- get_ENTREZID_rank(ranks.bif.lac, ens2ent)

get_GSEA <- function(ranks.df, TERM2GENE, 
                     eps, pvalueCutoff, minGSSize, maxGSSize){
    gsea.df <- GSEA(ranks.df, TERM2GENE = TERM2GENE, eps = eps, 
                      pvalueCutoff = pvalueCutoff, minGSSize = minGSSize, 
                      maxGSSize = maxGSSize)
    gsea.df <- as.data.frame(gsea.df)
    gsea.df <- gsea.df[abs(gsea.df$NES) > 1,]
    return(gsea.df)
}

set.seed(100)
cells.bif <- get_GSEA(ranks.bif.ent, TERM2GENE = cells[cells$cell_name %in% immune.cells,],
         eps = 0, pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 500)
cells.lac <- get_GSEA(ranks.lac.ent, TERM2GENE = cells[cells$cell_name %in% immune.cells,],
                      eps = 0, pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 500)
cells.bif.lac <- get_GSEA(ranks.bif.lac.ent, TERM2GENE = cells[cells$cell_name %in% immune.cells,],
                      eps = 0, pvalueCutoff = 0.05, minGSSize = 20, maxGSSize = 500)

get_subset2 <- function(df, group){
    df.sbs <- df[c(1,5)]
    df.sbs$group <- group
    return(df.sbs)
}

cells.bif.sbs <- get_subset2(cells.bif, "M vs M_BIF")
cells.lac.sbs <- get_subset2(cells.lac, "M vs M_LAC")
cells.bif.lac.sbs <- get_subset2(cells.bif.lac, "M_BIF vs M_LAC")

cells.all <- rbind(cells.bif.sbs, cells.lac.sbs, cells.bif.lac.sbs)
cells.all.s <- spread(cells.all, ID, NES, fill =  0)
rownames(cells.all.s) <- cells.all.s$group
cells.all.s <- cells.all.s[-1]

heatmap.immune.gsea <- pheatmap::pheatmap(t(cells.all.s), display_numbers = T, cutree_rows = 4)

# pdf(file = "figures/heatmap.immune.gsea.pdf", width = 4.75, height = 6.5)
# heatmap.immune.gsea
# dev.off()

# Immune cells proportion analysis: mmcp_counter
exp.data <- COU
# rownames(exp.data) <- sapply(str_split(rownames(exp.data), "\\."), function(x) x[1])

immune_cells <- deconvolute_mmcp_counter(exp.data, log2 = F, genome = "GCRm39", gene_id = "ENSEMBL.ID")

# write.table(immune_cells, "results/mmcp_counter.tsv", sep = "\t", row.names = T, quote = F)

column.annot <- meta_data
rownames(column.annot) <- column.annot$sampleid
column.annot <- column.annot[-c(1,2)]
column.annot$batch <- as.factor(column.annot$batch)

annot.colors <- list(group=c(M=mypal[1], M_BIF = mypal[2], M_LAC = mypal[3]), 
                     batch = c(`1` = "darkred", `2` = "darkblue"))

heatmap.immune <- pheatmap::pheatmap(log(immune_cells+0.001), 
                   scale = "column", 
                   annotation_col = column.annot, 
                   cutree_cols = 4, cutree_rows = 6, 
                   annotation_colors = annot.colors, 
                   show_colnames = F)

# pdf(file = "figures/heatmap.immune.pdf", width = 7, height = 3.75)
# heatmap.immune
# dev.off()

immune_cells.df <- as.data.frame(t(immune_cells))
immune_cells.df <- merge(meta_data, cbind(rownames(immune_cells.df), immune_cells.df), by = 1)[-2]

model <- lmer(`T cells` ~ group + (1|batch), data=immune_cells.df)
step_model <- step(model)
step_model$fixed

my_comparison <- list(c("M", "M_BIF"), c("M", "M_LAC"), c("M_BIF", "M_LAC"))

tcells_boxplot <- ggplot(immune_cells.df, aes(`T cells`, group, fill = group))+
    geom_boxplot()+
    theme_bw()+
    theme(legend.position = "none")+
    scale_fill_brewer(palette = "Set1")+
    stat_compare_means(comparisons = my_comparison)+
    xlab("Mmcp counter T cells score")+
    ylab("Group")

# ggsave(filename = "figures/tcells_boxplot.pdf", plot = tcells_boxplot, device = "pdf", width = 4, height = 2.5)

# Co-expression analysis
ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data[-c(4)],
                                       directory = "data/htseq/",
                                       design= ~ 0 + group)

final_exp <- exp_preprocess(
    ddsHTSeq2, min_exp = 10, variance_filter = TRUE, n = 4000
)

rownames(final_exp) <- sapply(str_split(rownames(final_exp), "\\."), function(x) x[1])

sft <- SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power

net <- exp2gcn(
    final_exp, net_type = "signed hybrid", SFTpower = power, 
    cor_method = "pearson"
)

plot_dendro_and_colors(net)

# ggsave(filename = "figures/co-expression/dendrogram_genes_modules.pdf", device = "pdf", width = 10, height = 5)

plot_ngenes_per_module(net)

# ggsave(filename = "figures/co-expression/ngenes_per_module.pdf", device = "pdf", width = 5, height = 4)

ms <- module_stability(final_exp, net, nRuns = 5)

# ggsave(filename = "figures/co-expression/module_stability.pdf", device = "pdf", width = 8, height = 5)

MEtrait <- module_trait_cor(exp = final_exp, MEs = net$MEs)

module_trait_cor_plot <- plot_module_trait_cor(MEtrait)

# pdf(file = "figures/co-expression/module_trait_cor_plot.pdf", width = 5, height = 4.5)
# module_trait_cor_plot
# dev.off()

plot_expression_profile(
    exp = final_exp, 
    net = net, 
    plot_module = TRUE, 
    modulename = "skyblue3"
) + scale_fill_brewer(palette = "Set1")

# ggsave(filename = "figures/co-expression/expression_profile.pdf", device = "pdf", width = 8, height = 5)

genes_and_modules <- net$genes_and_modules
skyblue3_modules <- genes_and_modules[genes_and_modules$Modules == "skyblue3",]

# write.table(genes_and_modules, "results/genes_and_modules.tsv", sep = "\t", quote = F, row.names = F)
# write.table(skyblue3_modules, "results/skyblue3_modules.tsv", sep = "\t", quote = F, row.names = F)

hubs <- get_hubs_gcn(final_exp, net)
hubs_skyblue3 <- hubs[hubs$Module == "skyblue3",]

# write.table(hubs_skyblue3, "results/hubs_skyblue3.tsv", sep = "\t", quote = F, row.names = F)

edges <- get_edge_list(net, module="skyblue3")
edges_filtered <- get_edge_list(net, module = "skyblue3", filter = TRUE, method = "optimalSFT", r_optimal_test = 0.7)

plot_gcn(
    edgelist_gcn = edges_filtered, 
    net = net, 
    color_by = "module", 
    hubs = hubs
)

# ggsave(filename = "figures/co-expression/skyblue3_network.pdf", device = "pdf", width = 5, height = 5)