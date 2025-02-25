setwd("/home/acari/Documents/github/probiotics-melanoma/")

# Import libraries
library(RColorBrewer)
library(tidyverse)
library(DESeq2)
library(limma)
library(vsn)
library(ggpubr)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)
library(ComplexHeatmap)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(vegan)
library(biomaRt)
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(Orthology.eg.db)
library(lmerTest)
library(permutes)

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
DEG_BIF <- cbind(as.data.frame(res_1), COU)
DEG_BIF <- DEG_BIF[!is.na(DEG_BIF$padj),]
DEG_BIF <- DEG_BIF[order(DEG_BIF$log2FoldChange, decreasing = T),]

DEG_LAC <- cbind(as.data.frame(res_2), COU)
DEG_LAC <- DEG_LAC[!is.na(DEG_LAC$padj),]
DEG_LAC <- DEG_LAC[order(DEG_LAC$log2FoldChange, decreasing = T),]

DEG_BIF_LAC <- cbind(as.data.frame(res_3), COU)
DEG_BIF_LAC <- DEG_BIF_LAC[!is.na(DEG_BIF_LAC$padj),]
DEG_BIF_LAC <- DEG_BIF_LAC[order(DEG_BIF_LAC$log2FoldChange, decreasing = T),]

# write.table(cbind(ensembl = rownames(DEG_BIF), DEG_BIF), "results/deg_bif.tsv", sep = "\t", quote = F, row.names = F)
# write.table(cbind(ensembl = rownames(DEG_LAC), DEG_LAC), "results/deg_lac.tsv", sep = "\t", quote = F, row.names = F)
# write.table(cbind(ensembl = rownames(DEG_BIF_LAC), DEG_BIF_LAC), "results/deg_bif_lac.tsv", sep = "\t", quote = F, row.names = F)

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
pathwaysDF <- msigdbr("mouse", category="H")
pathways <- split(as.character(pathwaysDF$ensembl_gene), pathwaysDF$gs_name)

ranks.bif <- DEG_BIF$log2FoldChange
names(ranks.bif) <- rownames(DEG_BIF)
names(ranks.bif) <- sapply(str_split(names(ranks.bif), "\\."), function(x) x[1])

ranks.lac <- DEG_LAC$log2FoldChange
names(ranks.lac) <- rownames(DEG_LAC)
names(ranks.lac) <- sapply(str_split(names(ranks.lac), "\\."), function(x) x[1])

ranks.bif.lac <- DEG_BIF_LAC$log2FoldChange
names(ranks.bif.lac) <- rownames(DEG_BIF_LAC)
names(ranks.bif.lac) <- sapply(str_split(names(ranks.bif.lac), "\\."), function(x) x[1])

set.seed(10)
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

fgsea.bif.sbs <- fgsea.bif[fgsea.bif$padj < 0.05,]
fgsea.lac.sbs <- fgsea.lac[fgsea.lac$padj < 0.05,]
fgsea.bif.lac.sbs <- fgsea.bif.lac[fgsea.bif.lac$padj < 0.05,]

fgsea.bif.sbs <- as.data.frame(fgsea.bif.sbs)[c(1,6)]
fgsea.bif.sbs$group <- "M vs M_BIF"
fgsea.lac.sbs <- as.data.frame(fgsea.lac.sbs)[c(1,6)]
fgsea.lac.sbs$group <- "M vs M_LAC"
fgsea.bif.lac.sbs <- as.data.frame(fgsea.bif.lac.sbs)[c(1,6)]
fgsea.bif.lac.sbs$group <- "M_BIF vs M_LAC"

fgsea.all <- rbind(fgsea.bif.sbs, fgsea.lac.sbs, fgsea.bif.lac.sbs)
fgsea.all <- spread(fgsea.all, group, NES, fill = 0)
rownames(fgsea.all) <- fgsea.all$pathway
fgsea.all <- fgsea.all[-1]

heatmap.hallmark <- pheatmap::pheatmap(as.matrix(fgsea.all), cutree_rows = 7, cutree_cols = 2, display_numbers = T)

# Immune cell types GSEA
cell_marker_data <- read.csv("data/cell_marker_mouse.csv",sep = "\t")
cell_names <- unique(cell_marker_data$cell_name)

# cells <- cell_marker_data %>%
#     dplyr::select(cell_name, GeneID)

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

DEG_BIF.rank <- as.data.frame(DEG_BIF[2])
DEG_BIF.rank <- merge(ens2ent, cbind(rownames(DEG_BIF.rank), DEG_BIF.rank), by = 1)
DEG_BIF.rank <- DEG_BIF.rank[!is.na(DEG_BIF.rank$ENTREZID),]
DEG_BIF.rank <- DEG_BIF.rank[!duplicated(DEG_BIF.rank$ENTREZID),]

ent.bif.rank <- DEG_BIF.rank$log2FoldChange
names(ent.bif.rank) <- DEG_BIF.rank$ENTREZID
ent.bif.rank <- sort(ent.bif.rank, decreasing = T)

cells.bif <- GSEA(ent.bif.rank, TERM2GENE = cells[cells$cell_name %in% immune.cells,], eps = 0, pvalueCutoff = 0.05)
cells.bif <- as.data.frame(cells.bif)
cells.bif <- cells.bif[abs(cells.bif$NES) > 1,]

DEG_LAC.rank <- as.data.frame(DEG_LAC[2])
DEG_LAC.rank <- merge(ens2ent, cbind(rownames(DEG_LAC.rank), DEG_LAC.rank), by = 1)
DEG_LAC.rank <- DEG_LAC.rank[!is.na(DEG_LAC.rank$ENTREZID),]
DEG_LAC.rank <- DEG_LAC.rank[!duplicated(DEG_LAC.rank$ENTREZID),]

ent.lac.rank <- DEG_LAC.rank$log2FoldChange
names(ent.lac.rank) <- DEG_LAC.rank$ENTREZID
ent.lac.rank <- sort(ent.lac.rank, decreasing = T)

cells.lac <- GSEA(ent.lac.rank, TERM2GENE = cells[cells$cell_name %in% immune.cells,], eps = 0, pvalueCutoff = 0.05)
cells.lac <- as.data.frame(cells.lac)
cells.lac <- cells.lac[abs(cells.lac$NES) > 1,]

# Immune cells proportion analysis: mmcp_counter
exp.data <- COU
rownames(exp.data) <- sapply(str_split(rownames(exp.data), "\\."), function(x) x[1])

immune_cells <- deconvolute_mmcp_counter(exp.data, log2 = F, genome = "GCRm39", gene_id = "ENSEMBL.ID")

column.annot <- meta_data
rownames(column.annot) <- column.annot$sampleid
column.annot <- column.annot[-c(1,2,4)]
column.annot$batch <- as.factor(column.annot$batch)

heatmap.immune <- pheatmap::pheatmap(log(immune_cells+0.001), 
                   scale = "column", 
                   annotation_col = column.annot, 
                   cutree_cols = 4, cutree_rows = 6)

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

# Co-expression analysis
library(BioNERO)

ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(sampleTable = meta_data[-c(4,5)],
                                       directory = "data/htseq/",
                                       design= ~ 0 + group)

final_exp <- exp_preprocess(
    ddsHTSeq2, min_exp = 10, variance_filter = TRUE, n = 2000
)

rownames(final_exp) <- sapply(str_split(rownames(final_exp), "\\."), function(x) x[1])

# plot_heatmap(final_exp, type = "samplecor", show_rownames = FALSE)
# plot_PCA(final_exp)+
#     ggtitle("PCA")

sft <- SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power

net <- exp2gcn(
    final_exp, net_type = "signed hybrid", SFTpower = power, 
    cor_method = "pearson"
)

plot_dendro_and_colors(net)
plot_eigengene_network(net)
plot_ngenes_per_module(net)

module_stability(final_exp, net, nRuns = 5)

MEtrait <- module_trait_cor(exp = final_exp, MEs = net$MEs)

plot_module_trait_cor(MEtrait)

plot_expression_profile(
    exp = final_exp, 
    net = net, 
    plot_module = TRUE, 
    modulename = "magenta"
)

genes_and_modules <- net$genes_and_modules
megenta_modules <- genes_and_modules[genes_and_modules$Modules == "magenta",]

# write.table(genes_and_modules, "results/genes_and_modules.tsv", sep = "\t", quote = F, row.names = F)
# write.table(megenta_modules, "results/megenta_modules.tsv", sep = "\t", quote = F, row.names = F)

hubs <- get_hubs_gcn(final_exp, net)
hubs_magenta <- hubs[hubs$Module == "magenta",]

# write.table(hubs_magenta, "results/hubs_magenta.tsv", sep = "\t", quote = F, row.names = F)

edges <- get_edge_list(net, module="magenta")
edges_filtered <- get_edge_list(net, module = "magenta", filter = TRUE, method = "optimalSFT")

plot_gcn(
    edgelist_gcn = edges_filtered, 
    net = net, 
    color_by = "module", 
    hubs = hubs
)