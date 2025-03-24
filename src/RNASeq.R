setwd("/home/acari/Documents/github/probiotics-melanoma/")

# Import libraries
library(vegan)
library(fgsea)
library(limma)
library(stringr)
library(BioNERO)
library(DESeq2)
library(tidyverse)
library(factoextra)
library(clusterProfiler)

library(ggpubr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(EnhancedVolcano)

library(DOSE)
library(msigdbr)
library(msigdbdf)
library(biomaRt)
library(ReactomePA)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(AnnotationHub)
library(GSVA)
library(STRINGdb)
library(immunedeconv)

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

ENTREZID <- mapIds(org.Mm.eg.db, keys = rownames(ddsHTSeq), keytype="ENSEMBL", column = "ENTREZID")
ENTREZID <- as.data.frame(ENTREZID)
ENTREZID <- cbind(rownames(ENTREZID), ENTREZID)
colnames(ENTREZID)[1] <- "ENSEMBL"

SYMBOL <- mapIds(org.Mm.eg.db, keys = rownames(ddsHTSeq), keytype="ENSEMBL", column = "SYMBOL")
SYMBOL <- as.data.frame(SYMBOL)
SYMBOL <- cbind(rownames(SYMBOL), SYMBOL)
colnames(SYMBOL)[1] <- "ENSEMBL"

convertID <- merge(SYMBOL, ENTREZID)

# Filter gene table by counts
keep <- rowSums(counts(ddsHTSeq)>= 10) > ncol(ddsHTSeq)*0.3
ddsHTSeq <- ddsHTSeq[keep,]

# Diff-exp analysis
ddsHTSeq <- DESeq(ddsHTSeq)

resultsNames(ddsHTSeq)

COU <- counts(ddsHTSeq, normalized=TRUE)
# write.table(COU, "results/normalized.counts.tsv", sep = "\t", quote = F, row.names = T)

contr_1 <- makeContrasts(groupM_BIF - groupM, levels = resultsNames(ddsHTSeq))
contr_2 <- makeContrasts(groupM_LAC - groupM, levels = resultsNames(ddsHTSeq))
contr_3 <- makeContrasts(groupM_BIF - groupM_LAC, levels = resultsNames(ddsHTSeq))

res_1 <- results(ddsHTSeq, contrast=contr_1)
res_2 <- results(ddsHTSeq, contrast=contr_2)
res_3 <- results(ddsHTSeq, contrast=contr_3)

# Make PCA visualization
vsd <- as.data.frame(assay(vst(ddsHTSeq)))

PCA <- prcomp(t(vsd))
PCA.points <- as.data.frame(PCA$x)
PCA.points <- merge(meta_data, cbind(rownames(PCA.points), PCA.points), by = 1)[-2]

PCA_plot <- ggplot(PCA.points, aes(PC1, PC2, col = group, shape = batch))+
    geom_point(size = 2.75)+
    # geom_text(size = 3)+
    theme_bw()+
    theme(legend.position = "right")+
    scale_color_brewer(palette = "Set1")+
    xlab("PC1 (48%)")+
    ylab("PC2 (17%)")

prop_plot <- fviz_eig(PCA, col.var="blue") + ggtitle("Proportion of variance")

# ggsave(filename = "figures/nmds.plot.pdf", plot = nmds.plot, device = "pdf", width = 5, height = 5)
# ggsave(filename = "figures/PCA_plot.pdf", plot = PCA_plot, device = "pdf", width = 6, height = 5)

# Make DEG tables
get_DEG <- function(res){
    # DEG <- cbind(as.data.frame(res), COU)
    DEG <- as.data.frame(res)
    DEG <- DEG[!is.na(DEG$padj),]
    DEG <- DEG[order(DEG$log2FoldChange, decreasing = T),]
    return(DEG)
}

DEG_BIF <- get_DEG(res_1)
DEG_LAC <- get_DEG(res_2)
DEG_BIF_LAC <- get_DEG(res_3)

DEG_BIF_add <- merge(convertID[-3], cbind(rownames(DEG_BIF), DEG_BIF), by = 1)
DEG_LAC_add <- merge(convertID[-3], cbind(rownames(DEG_LAC), DEG_LAC), by = 1)
DEG_BIF_LAC_add <- merge(convertID[-3], cbind(rownames(DEG_BIF_LAC), DEG_BIF_LAC), by = 1)

DEG_BIF_add <- DEG_BIF_add[order(DEG_BIF_add$log2FoldChange, decreasing = T),]
rownames(DEG_BIF_add) <- 1:nrow(DEG_BIF_add)

DEG_LAC_add <- DEG_LAC_add[order(DEG_LAC_add$log2FoldChange, decreasing = T),]
rownames(DEG_LAC_add) <- 1:nrow(DEG_LAC_add)

DEG_BIF_LAC_add <- DEG_BIF_LAC_add[order(DEG_BIF_LAC_add$log2FoldChange, decreasing = T),]
rownames(DEG_BIF_LAC_add) <- 1:nrow(DEG_BIF_LAC_add)

# Get volcano plots
VolcanoPlot1 <- EnhancedVolcano(DEG_BIF_add,
                lab = DEG_BIF_add$SYMBOL,
                x = 'log2FoldChange', 
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M + Bif vs M", 
                ###########################
                # labSize = 4,
                drawConnectors = TRUE, 
                # widthConnectors = 0.75, 
                gridlines.major = FALSE, 
                gridlines.minor = FALSE)

VolcanoPlot2 <- EnhancedVolcano(DEG_LAC_add,
                lab = DEG_LAC_add$SYMBOL,
                x = 'log2FoldChange', 
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M + Lac vs M", 
                ###########################
                # labSize = 4,
                drawConnectors = TRUE, 
                # widthConnectors = 0.75, 
                gridlines.major = FALSE, 
                gridlines.minor = FALSE)

VolcanoPlot3 <- EnhancedVolcano(DEG_BIF_LAC_add,
                lab = DEG_BIF_LAC_add$SYMBOL,
                x = 'log2FoldChange', 
                y = 'padj', pCutoff = 0.05, 
                subtitle = "", 
                title = "M + Bif vs M + Lac", 
                ###########################
                # labSize = 4,
                drawConnectors = TRUE, 
                # widthConnectors = 0.75, 
                gridlines.major = FALSE, 
                gridlines.minor = FALSE)

# GSEA
## MSigDB GSEA
### Hallmark
pathwaysH <- msigdbr("mouse", category="H")
pathwaysH.sbs <- split(as.character(pathwaysH$ensembl_gene), pathwaysH$gs_name)

get_ranks <- function(DEG){
    DEG$log2FoldChange[DEG$log2FoldChange > 0] <- 1
    DEG$log2FoldChange[DEG$log2FoldChange < 0] <- -1
    ranks.df <- -log10(DEG$pvalue)*DEG$log2FoldChange
    names(ranks.df) <- rownames(DEG)
    # names(ranks.df) <- sapply(str_split(names(ranks.df), "\\."), function(x) x[1])    
    ranks.df <- sort(ranks.df, decreasing = T)
    return(ranks.df)
}

ranks.bif <- get_ranks(DEG_BIF)
ranks.lac <- get_ranks(DEG_LAC)
ranks.bif.lac <- get_ranks(DEG_BIF_LAC)

fgsea.bif <- fgsea(pathways = pathwaysH.sbs, 
                  stats    = ranks.bif,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

fgsea.lac <- fgsea(pathways = pathwaysH.sbs, 
                     stats    = ranks.lac,
                     eps      = 0.0,
                     minSize  = 15,
                     maxSize  = 500)

fgsea.bif.lac <- fgsea(pathways = pathwaysH.sbs, 
                   stats    = ranks.bif.lac,
                   eps      = 0.0,
                   minSize  = 15,
                   maxSize  = 500)

filter_pval <- function(fgsea, NES, padj){
    fgsea <- as.data.frame(fgsea)
    fgsea <- fgsea[fgsea$padj < padj,]
    fgsea <- fgsea[abs(fgsea$NES) > NES,]
    fgsea <- fgsea[order(fgsea$NES, decreasing = T),]
    return(fgsea)
}

fgsea.bif.filter <- filter_pval(fgsea.bif, 1, 0.01)
fgsea.lac.filter <- filter_pval(fgsea.lac, 1, 0.01)
fgsea.bif.lac.filter <- filter_pval(fgsea.bif.lac, 1, 0.01)

# write.table(fgsea.bif.filter[-8], "results/fgsea.bif.tsv", sep = "\t", quote = F, row.names = F)
# write.table(fgsea.lac.filter[-8], "results/fgsea.lac.tsv", sep = "\t", quote = F, row.names = F)
# write.table(fgsea.bif.lac.filter[-8], "results/fgsea.bif.lac.tsv", sep = "\t", quote = F, row.names = F)

get_subset <- function(df, group){
    df.sbs <- df[c(1,6)]
    df.sbs$group <- group
    return(df.sbs)
}

fgsea.bif.sbs <- get_subset(fgsea.bif.filter, "M+BIF vs M")
fgsea.lac.sbs <- get_subset(fgsea.lac.filter, "M+LAC vs M ")
fgsea.bif.lac.sbs <- get_subset(fgsea.bif.lac.filter, "M+BIF vs M+LAC")

fgsea.all <- rbind(fgsea.bif.sbs, fgsea.lac.sbs, fgsea.bif.lac.sbs)
fgsea.all <- spread(fgsea.all, group, NES, fill = 0)
rownames(fgsea.all) <- fgsea.all$pathway
fgsea.all <- fgsea.all[-1]

heatmap.hallmark <- pheatmap::pheatmap(as.matrix(fgsea.all), cutree_rows = 7, cutree_cols = 2, display_numbers = T)

# Cell Marker GSEA
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

cells.bif.sbs <- get_subset2(cells.bif, "M+BIF vs M")
cells.lac.sbs <- get_subset2(cells.lac, "M+LAC vs M")
cells.bif.lac.sbs <- get_subset2(cells.bif.lac, "M+BIF vs M+LAC")

cells.all <- rbind(cells.bif.sbs, cells.lac.sbs, cells.bif.lac.sbs)
cells.all.s <- spread(cells.all, ID, NES, fill =  0)
rownames(cells.all.s) <- cells.all.s$group
cells.all.s <- cells.all.s[-1]

heatmap.immune.gsea <- pheatmap::pheatmap(t(cells.all.s), display_numbers = T, cutree_rows = 4)

# T cell Mmcp counter
exp.data <- vst(counts(ddsHTSeq))
immune_cells <- deconvolute_mmcp_counter(exp.data, log2 = F, genome = "GCRm39", gene_id = "ENSEMBL.ID")

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

immune_cells.df <- as.data.frame(t(immune_cells))
immune_cells.df <- merge(meta_data, cbind(rownames(immune_cells.df), immune_cells.df), by = 1)[-2]
immune_cells.df$`T cells` <- rank(immune_cells.df$`T cells`)

my_comparison <- list(c("M", "M_BIF"), c("M", "M_LAC"), c("M_BIF", "M_LAC"))

tcells_boxplot <- ggplot(immune_cells.df, aes(group, `T cells`, fill = group))+
    geom_boxplot()+
    theme_bw()+
    theme(legend.position = "none")+
    scale_fill_brewer(palette = "Set1")+
    stat_compare_means(comparisons = my_comparison)+
    xlab("Mmcp counter T cells score")+
    ylab("Group")

summary(aov(`T cells` ~ group + batch, data = immune_cells.df))

# CD8+ ssGSEA
tcell_names <- cell$cell_name[!is.na(str_extract(cell$cell_name, "CD8"))]
tcell_names <- unique(tcell_names)

CD8 <- cell$Symbol[cell$cell_name == "CD8+ T cell"]
CD8 <- CD8[CD8 != ""]

CD8_sets <- list(
    CD8 = CD8
)

gsvapar_2 <- gsvaParam(as.matrix(df.Symbol), CD8_sets)
ssgsea_scores_2 <- gsva(gsvapar_2)
ssgsea_scores_2 <- as.data.frame(t(ssgsea_scores_2))

ssgsea_scores_2$CD8 <- rank(ssgsea_scores_2$CD8)

ssgsea_scores_2 <- merge(meta_data, cbind(rownames(ssgsea_scores_2), ssgsea_scores_2), by = 1)[-2]

ggplot(ssgsea_scores_2, aes(group, CD8))+
    geom_boxplot()+
    stat_compare_means(method = "wilcox", 
                       comparisons = comparisons)

summary(aov(CD8 ~ group + batch, data = ssgsea_scores_2))

# Macrophages ssGSEA
df.Symbol <- as.data.frame(counts(ddsHTSeq))
df.Symbol <- merge(convertID[c(1:2)], cbind(rownames(df.Symbol), df.Symbol), by = 1)
df.Symbol <- df.Symbol[!is.na(df.Symbol$SYMBOL),]
df.Symbol <- df.Symbol[!duplicated(df.Symbol$SYMBOL),]
df.Symbol <- df.Symbol[-1]

rownames(df.Symbol) <- df.Symbol$SYMBOL
df.Symbol <- df.Symbol[-1]

cell <- cell_marker_data[c("cell_name", "Symbol")]

macrophage_names <- cell$cell_name[!is.na(str_extract(cell$cell_name, "macrophage|Macrophage"))]
macrophage_names <- macrophage_names[!is.na(str_extract(macrophage_names, "M1|M2|Macrophage"))]
macrophage_names <- macrophage_names[is.na(str_extract(macrophage_names, "like|inflammatory|Fn1+|Pf4+"))]
macrophage_names <- unique(macrophage_names)
macrophage_names <- macrophage_names[c(1:3)]

cell_macro <- cell[cell$cell_name %in% macrophage_names,]
cell_macro <- cell_macro[cell_macro$Symbol != "",]
cell_macro <- unique(cell_macro)

macro_sets <- list(
    M1 = cell_macro$Symbol[cell_macro$cell_name == "M1 macrophage"],
    M2 = cell_macro$Symbol[cell_macro$cell_name == "M2 macrophage"]
)

gsvapar <- gsvaParam(as.matrix(df.Symbol), macro_sets)
ssgsea_scores <- gsva(gsvapar)

ssgsea_results <- as.data.frame(t(ssgsea_scores))

ssgsea_results$M1_M2_ratio <- ssgsea_results$M1 / ssgsea_results$M2
ssgsea_results <- merge(meta_data, cbind(rownames(ssgsea_results), ssgsea_results), by = 1)[-2]

ssgsea_results$M1 <- rank(ssgsea_results$M1)
ssgsea_results$M2 <- rank(ssgsea_results$M2)

comparisons <- list(c("M", "M_BIF"),  c("M", "M_LAC"), c("M_BIF", "M_LAC"))

ggplot(ssgsea_results, aes(group, M1))+
    geom_boxplot()+
    stat_compare_means(method = "wilcox", 
                       comparisons = comparisons)

ggplot(ssgsea_results, aes(group, M2))+
    geom_boxplot()+
    stat_compare_means(method = "wilcox", 
                       comparisons = comparisons)

ggplot(ssgsea_results[ssgsea_results$M1_M2_ratio < 9,], aes(group, M1_M2_ratio))+
    geom_boxplot()+
    stat_compare_means(method = "wilcox", 
                       comparisons = comparisons)

summary(aov(M1_M2_ratio ~ group + batch, data = ssgsea_results[ssgsea_results$M1_M2_ratio < 9,]))
summary(aov(M1 ~ group + batch, data = ssgsea_results[ssgsea_results$M1_M2_ratio < 9,]))
summary(aov(M2 ~ group + batch, data = ssgsea_results[ssgsea_results$M1_M2_ratio < 9,]))

# Co-expression analysis
## Batch correction
vsd_2 <- vst(ddsHTSeq, blind = F)
mat <- removeBatchEffect(assay(vsd_2), batch = meta_data$batch)

## PCA 
PCA_corr <- prcomp(t(mat))
PCA_corr.points <- as.data.frame(PCA_corr$x)
PCA_corr.points <- merge(meta_data, cbind(rownames(PCA_corr.points), PCA_corr.points), by = 1)[-2]

PCA_corr_plot <- ggplot(PCA_corr.points, aes(PC1, PC2, col = group, shape = batch))+
    geom_point(size = 2.75)+
    # geom_text(size = 3)+
    theme_bw()+
    theme(legend.position = "right")+
    scale_color_brewer(palette = "Set1")+
    xlab("PC1 (33%)")+
    ylab("PC2 (16%)")

prop_corr_plot <- fviz_eig(PCA_corr, col.var="blue") + ggtitle("Proportion of variance")

## Co-expression net
final_exp <- exp_preprocess(
    mat, min_exp = 10, variance_filter = TRUE, n = 4000
)

sft <- SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power

net <- exp2gcn(
    final_exp, net_type = "signed hybrid", SFTpower = power, 
    cor_method = "pearson"
)

plot_dendro <- plot_dendro_and_colors(net)
plot_cor_module <- plot_eigengene_network(net)
plot_ngenes <- plot_ngenes_per_module(net)

dfg <- meta_data[c(1,3)]
rownames(dfg) <- dfg$sampleid
dfg <- dfg[-1]

MEtrait <- module_trait_cor(exp = final_exp, MEs = net$MEs, metadata = dfg)
module_trait_cor_plot <- plot_module_trait_cor(MEtrait)

cor_sets <- c("greenyellow", "sienna3", "pink", "royalblue")

expression_plot_1 <- plot_expression_profile(
    exp = final_exp, 
    net = net, 
    metadata = dfg,
    plot_module = TRUE, 
    modulename = cor_sets[1]
) + scale_fill_brewer(palette = "Set1") +
    ggtitle("greenyellow")

expression_plot_2 <- plot_expression_profile(
    exp = final_exp, 
    net = net, 
    metadata = dfg,
    plot_module = TRUE, 
    modulename = cor_sets[2]
) + scale_fill_brewer(palette = "Set1") +
    ggtitle("sienna3")

expression_plot_3 <- plot_expression_profile(
    exp = final_exp, 
    net = net, 
    metadata = dfg,
    plot_module = TRUE, 
    modulename = cor_sets[3]
) + scale_fill_brewer(palette = "Set1") +
    ggtitle("pink")

expression_plot_4 <- plot_expression_profile(
    exp = final_exp, 
    net = net, 
    metadata = dfg,
    plot_module = TRUE, 
    modulename = cor_sets[4]
) + scale_fill_brewer(palette = "Set1") +
    ggtitle("royalblue")

expression_plot <- ggarrange(expression_plot_1, expression_plot_2, expression_plot_3, expression_plot_4, ncol = 1, common.legend = T)

# Genes and modules
genes_and_modules <- net$genes_and_modules
selected_modules <- genes_and_modules[genes_and_modules$Modules %in% cor_sets,]
selected_modules <- selected_modules[order(selected_modules$Modules, decreasing = T),]
rownames(selected_modules) <- 1:nrow(selected_modules)

selected_modules_add <- merge(convertID, selected_modules, by = 1)
selected_modules_add <- selected_modules_add[order(selected_modules_add$Modules, decreasing = T),]

colnames(selected_modules_add)[4] <- "modules"

# ORA
module_fix <- function(module, module_name){
    if (nrow(module) > 0){
        module$Category <- module_name
        module <- module[c("Category", "ID",
                           "Description", 
                           "RichFactor",
                           "qvalue", 
                           "Count")]
    } else {
        
        module <- NULL    
    }
    
    # module <- module[module$RichFactor > 0.1,]
    return(module)
}

enrich_all_ent <- function(genes, padj){
    KEGG.module <- enrichKEGG(genes, organism = "mmu", qvalueCutoff = padj)
    KEGG.module <- as.data.frame(KEGG.module)
    KEGG.module <- module_fix(KEGG.module, "KEGG")
    
    Reactome.module <- enrichPathway(genes, organism = "mouse", qvalueCutoff = padj)
    Reactome.module <- as.data.frame(Reactome.module)
    Reactome.module <- module_fix(Reactome.module, "Reactome")
    
    WP.module <- enrichWP(genes, organism = "Mus musculus")
    WP.module <- as.data.frame(WP.module)
    WP.module <- WP.module[WP.module$qvalue < padj,]
    WP.module <- module_fix(WP.module, "WikiPathways")
    
    ent_res <- rbind(KEGG.module, WP.module,  Reactome.module)
    
    return(ent_res)
}

enrich_all_ens <- function(genes, padj){
    
    GO.module <- enrichGO(genes, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENSEMBL", 
                          ont = "ALL", 
                          qvalueCutoff = padj, minGSSize = 10)
    GO.module <- as.data.frame(GO.module)
    GO.module <- module_fix(GO.module, "GO")
    
    H.module <- enricher(genes, 
                         TERM2GENE = pathwaysH[c("gs_name", "ensembl_gene")], 
                         qvalueCutoff = padj, minGSSize = 10)
    H.module <- as.data.frame(H.module)
    H.module <- module_fix(H.module, "MSigDB: HALLMARK")
    
    enrich.df <- rbind(GO.module, H.module)
    
    return(enrich.df)
}

ens_greenyellow_1 <- selected_modules_add$ENSEMBL[selected_modules_add$modules == "greenyellow"]
ent_greenyellow_2 <- selected_modules_add$ENTREZID[selected_modules_add$modules == "greenyellow"]

ens_sienna3_1 <- selected_modules_add$ENSEMBL[selected_modules_add$modules == "sienna3"]
ent_sienna3_2 <- selected_modules_add$ENTREZID[selected_modules_add$modules == "sienna3"]

ens_pink_1 <- selected_modules_add$ENSEMBL[selected_modules_add$modules == "pink"]
ent_pink_2 <- selected_modules_add$ENTREZID[selected_modules_add$modules == "pink"]

ens_royalblue_1 <- selected_modules_add$ENSEMBL[selected_modules_add$modules == "royalblue"]
ent_royalblue_2 <- selected_modules_add$ENTREZID[selected_modules_add$modules == "royalblue"]

enrich_1 <- rbind(enrich_all_ens(c(ens_greenyellow_1, ens_sienna3_1), 0.05), 
                              enrich_all_ent(c(ent_greenyellow_2, ent_sienna3_2), 0.05))
rownames(enrich_1) <- 1:nrow(enrich_1)

enrich_2 <- rbind(enrich_all_ens(c(ens_pink_1, ens_royalblue_1), 0.05), 
                  enrich_all_ent(c(ent_pink_2, ent_royalblue_2), 0.05))
rownames(enrich_2) <- 1:nrow(enrich_2)

## Hub genes
hubs <- get_hubs_gcn(final_exp, net)
selected_hubs <- hubs[hubs$Module %in% cor_sets,]
rownames(selected_hubs) <- 1:nrow(selected_hubs)

selected_hubs_add <- merge(convertID[-3], selected_hubs, by = 1)
selected_hubs_add <- selected_hubs_add[order(selected_hubs_add$Module, decreasing = T),]
rownames(selected_hubs_add) <- 1:nrow(selected_hubs_add)

## Net plot with gene symbol names
sbs <- convertID[convertID$ENSEMBL %in% rownames(final_exp),][c(1,2)]
sbs <- sbs[!is.na(sbs$SYMBOL),]

final_exp_2 <- as.data.frame(final_exp)
final_exp_2 <- final_exp_2[rownames(final_exp_2) %in% sbs$ENSEMBL,]

sbs <- sbs[match(rownames(final_exp_2), sbs$ENSEMBL),]

final_exp_2 <- final_exp
rownames(final_exp_2) <- sbs$SYMBOL

sft <- SFT_fit(final_exp_2, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power

net <- exp2gcn(
    final_exp_2, net_type = "signed hybrid", SFTpower = power, 
    cor_method = "pearson"
)

MEtrait <- module_trait_cor(exp = final_exp_2, MEs = net$MEs, metadata = dfg)
module_trait_cor_plot <- plot_module_trait_cor(MEtrait)

hubs <- get_hubs_gcn(final_exp_2, net)

edges_filtered_1 <- get_edge_list(net, 
                                  module = cor_sets[1],
                                  filter = TRUE,
                                  method = "optimalSFT",
                                  r_optimal_test = 0.7)

net_plot_1 <- plot_gcn(
    edgelist_gcn = edges_filtered_1, 
    net = net, 
    color_by = "module", 
    hubs = hubs,
    top_n_hubs = 5 
) + ggtitle("greenyellow")

edges_filtered_2 <- get_edge_list(net, 
                                  module = cor_sets[2],
                                  filter = TRUE,
                                  method = "optimalSFT",
                                  r_optimal_test = 0.7)

net_plot_2 <- plot_gcn(
    edgelist_gcn = edges_filtered_2, 
    net = net, 
    color_by = "module", 
    hubs = hubs,
    top_n_hubs = 5 
) + ggtitle("sienna3")

edges_filtered_3 <- get_edge_list(net, 
                                  module = cor_sets[3],
                                  filter = TRUE,
                                  method = "optimalSFT",
                                  r_optimal_test = 0.7)

net_plot_3 <- plot_gcn(
    edgelist_gcn = edges_filtered_3, 
    net = net, 
    color_by = "module", 
    hubs = hubs,
    top_n_hubs = 5 
) + ggtitle("pink")

edges_filtered_4 <- get_edge_list(net, 
                                  module = cor_sets[4],
                                  filter = TRUE,
                                  method = "optimalSFT",
                                  r_optimal_test = 0.7)

net_plot_4 <- plot_gcn(
    edgelist_gcn = edges_filtered_4, 
    net = net, 
    color_by = data.frame(geneid = unique(c(as.character(edges_filtered_4$Var1), as.character(edges_filtered_4$Var2))), color = "black"), 
    hubs = hubs,
    top_n_hubs = 5 
) + ggtitle("royalblue")

net_plot <- ggarrange(net_plot_1, net_plot_2, net_plot_3, net_plot_4)

## STRING analysis
string_db <- STRINGdb$new(version = "12.0", species = 10090, score_threshold = 100)

mapped_1 <- string_db$map(selected_modules_add[selected_modules_add$modules %in% c("greenyellow", "sienna3"),], "ENSEMBL", removeUnmappedRows = T)
mapped_2 <- string_db$map(selected_modules_add[selected_modules_add$modules %in% c("pink", "royalblue"),], "ENSEMBL", removeUnmappedRows = T)

string_db$plot_network(mapped_1$STRING_id, required_score = 400)
string_db$plot_network(mapped_2$STRING_id, required_score = 400)