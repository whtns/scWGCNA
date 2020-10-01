#!/usr/bin/env Rscript 

library(seuratTools)
library(Seurat)
library(tidyverse)
library(fs)
library(plotly)
library(clusterProfiler)

# https://github.com/milescsmith/scWGCNA
library(scWGCNA)
library(WGCNA)
library(ComplexHeatmap)

WGCNA::enableWGCNAThreads(nThreads = 6)

pseudotimes <- "results/20171031-SHL-FACS-Hs_pseudotimes.csv" %>% 
  read_csv() %>% 
  tibble::column_to_rownames("sample_id") %>%
  identity()

proj_dir <- "~/single_cell_projects/sc_RB_devel/20171031-SHL-FACS-Hs_proj/"

seu <- readRDS("~/single_cell_projects/sc_RB_devel/20171031-SHL-FACS-Hs_proj/output/seurat/allfeatures_seu.rds")$gene

pt1_seu <- seu[,!is.na(seu$pt1_clusters)]

pt2_seu <- seu[,!is.na(seu$pt2_clusters)]

pt_cells <- union(colnames(pt1_seu), colnames(pt2_seu))

all_pt_seu <- seu[,pt_cells]

# all_pt ------------------------------

Idents(all_pt_seu) <- all_pt_seu$custom_cluster

# scWGCNA_out <- scWGCNA(all_pt_seu)
# 
# saveRDS(scWGCNA_out, "output/scwgcna/scWGCNA_all_pt_seu.rds")
scWGCNA_out <- readRDS("output/scwgcna/scWGCNA_all_pt_seu.rds")

dissTom <- calculate_tom(scWGCNA_out$adjacency.matrix, scWGCNA_out$exprMatrix)

dandc <- makeDendroandColors(dissTom, min.module.size = 400)

modules <- extract_modules(dandc$dynamicMods, scWGCNA_out$exprMatrix)

module_df <- 
  modules %>% 
	tibble::enframe("module", "gene") %>% 
	tidyr::unnest()

intramodconnect <- intramodularConnectivity(scWGCNA_out$adjacency.matrix, labels2colors(dandc$dynamicMods)) %>% 
	tibble::rownames_to_column("gene") %>% 
	dplyr::left_join(module_df, by = "gene") %>% 
	dplyr::group_by(module) %>%
  dplyr::mutate(scaled_kWithin = kWithin/dplyr::n()) %>% 
	dplyr::arrange(module, desc(kWithin)) %>% 
	dplyr::mutate(gene = forcats::fct_reorder(gene, module)) %>% 
  dplyr::slice_max(order_by = kWithin, prop = 0.05) %>% 
	dplyr::mutate(top_5 = ifelse(kWithin > quantile(kWithin, 0.95), 1, 0)) %>%
	identity()

mod_cols <- unique(module_df$module)
names(mod_cols) <- mod_cols

intramod_plot <- ggplot(intramodconnect, aes(gene, scaled_kWithin, color = module)) +
	geom_point() +
	# gghighlight(top_5 > 0, keep_scales = TRUE) +
	scale_color_manual(values = mod_cols) +
	labs(title = "intramodular connectivity") +
	NULL


module_entrez <- map(modules, ~bitr(.x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

ggo_all <- map(module_entrez, ~enrichGO(gene = as.character(.x$ENTREZID), 
                                        OrgDb    = org.Hs.eg.db,
                                        ont      = "ALL",
                                        readable = TRUE))

plot_module <- function(module_go, module_name){
  dotplot(module_go) + 
    labs(title = module_name)
  
}

module_plots <- imap(ggo_all, plot_module)

MEList <- moduleEigengenes(
  expr = scWGCNA_out$exprMatrix,
  colors = labels2colors(dandc$dynamicMods),
  softPower = 4
)
MEs <- MEList$eigengenes

pdf("results/scwgcna/scwgcna_all_pt_seu.pdf", width = 20)
scWGCNA_out$plots[[1]]
scWGCNA_out$plots[[2]]
intramod_plot
pdc <- plotDendroAndColors(dandc$geneTree,
													 labels2colors(dandc$dynamicMods),
													 "Dynamic Tree Cut",
													 dendroLabels = FALSE,
													 hang = 0.03,
													 addGuide = TRUE,
													 guideHang = 0.05,
													 main = "Gene dendrogram and module colors"
)
# plot module eigengenes
Heatmap(MEs, show_row_names = FALSE)
module_plots
dev.off()

# pt1 ------------------------------

Idents(pt1_seu) <- pt1_seu$pt1_clusters

# scWGCNA_out <- scWGCNA(pt1_seu)
# 
# saveRDS(scWGCNA_out, "output/scwgcna/scWGCNA_pt1_seu.rds")
scWGCNA_out <- readRDS("output/scwgcna/scWGCNA_pt1_seu.rds")

dissTom <- calculate_tom(scWGCNA_out$adjacency.matrix, scWGCNA_out$exprMatrix)

dandc <- makeDendroandColors(dissTom, min.module.size = 50)

modules <- extract_modules(dandc$dynamicMods, scWGCNA_out$exprMatrix)

module_df <- 
  modules %>% 
  tibble::enframe("module", "gene") %>% 
  tidyr::unnest()

intramodconnect <- intramodularConnectivity(scWGCNA_out$adjacency.matrix, labels2colors(dandc$dynamicMods)) %>% 
  tibble::rownames_to_column("gene") %>% 
  dplyr::left_join(module_df, by = "gene") %>% 
  dplyr::group_by(module) %>%
  dplyr::mutate(scaled_kWithin = kWithin/dplyr::n()) %>% 
  dplyr::arrange(module, desc(kWithin)) %>% 
  dplyr::mutate(gene = forcats::fct_reorder(gene, module)) %>% 
  dplyr::slice_max(order_by = kWithin, prop = 0.05) %>% 
  dplyr::mutate(top_5 = ifelse(kWithin > quantile(kWithin, 0.95), 1, 0)) %>%
  identity()

mod_cols <- unique(module_df$module)
names(mod_cols) <- mod_cols

intramod_plot <- ggplot(intramodconnect, aes(gene, scaled_kWithin, color = module)) +
  geom_point() +
  # gghighlight(top_5 > 0, keep_scales = TRUE) +
  scale_color_manual(values = mod_cols) +
  labs(title = "intramodular connectivity") +
  NULL


module_entrez <- map(modules, ~bitr(.x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

ggo_all <- map(module_entrez, ~enrichGO(gene = as.character(.x$ENTREZID), 
                                        OrgDb    = org.Hs.eg.db,
                                        ont      = "ALL",
                                        readable = TRUE))

plot_module <- function(module_go, module_name){
  dotplot(module_go) + 
    labs(title = module_name)
  
}

module_plots <- imap(ggo_all, plot_module)

MEList <- moduleEigengenes(
  expr = scWGCNA_out$exprMatrix,
  colors = labels2colors(dandc$dynamicMods),
  softPower = 3
)
MEs <- MEList$eigengenes

pdf("results/scwgcna/scwgcna_pt1_seu.pdf", width = 20)
scWGCNA_out$plots[[1]]
scWGCNA_out$plots[[2]]
intramod_plot
pdc <- plotDendroAndColors(dandc$geneTree,
                           labels2colors(dandc$dynamicMods),
                           "Dynamic Tree Cut",
                           dendroLabels = FALSE,
                           hang = 0.03,
                           addGuide = TRUE,
                           guideHang = 0.05,
                           main = "Gene dendrogram and module colors"
)
# plot module eigengenes
Heatmap(MEs, show_row_names = FALSE)
module_plots
dev.off()


# pt2 ------------------------------

Idents(pt2_seu) <- pt2_seu$pt2_clusters

# scWGCNA_out <- scWGCNA(pt2_seu)
# 
# saveRDS(scWGCNA_out, "output/scwgcna/scWGCNA_pt2_seu.rds")
scWGCNA_out <- readRDS("output/scwgcna/scWGCNA_pt2_seu.rds")

dissTom <- calculate_tom(scWGCNA_out$adjacency.matrix, scWGCNA_out$exprMatrix)

dandc <- makeDendroandColors(dissTom, min.module.size = 800)

modules <- extract_modules(dandc$dynamicMods, scWGCNA_out$exprMatrix)

module_df <- 
  modules %>% 
  tibble::enframe("module", "gene") %>% 
  tidyr::unnest()

intramodconnect <- intramodularConnectivity(scWGCNA_out$adjacency.matrix, labels2colors(dandc$dynamicMods)) %>% 
  tibble::rownames_to_column("gene") %>% 
  dplyr::left_join(module_df, by = "gene") %>% 
  dplyr::group_by(module) %>%
  dplyr::mutate(scaled_kWithin = kWithin/dplyr::n()) %>% 
  dplyr::arrange(module, desc(kWithin)) %>% 
  dplyr::mutate(gene = forcats::fct_reorder(gene, module)) %>% 
  dplyr::slice_max(order_by = kWithin, prop = 0.05) %>% 
  dplyr::mutate(top_5 = ifelse(kWithin > quantile(kWithin, 0.95), 1, 0)) %>%
  identity()

mod_cols <- unique(module_df$module)
names(mod_cols) <- mod_cols

intramod_plot <- ggplot(intramodconnect, aes(gene, scaled_kWithin, color = module)) +
  geom_point() +
  # gghighlight(top_5 > 0, keep_scales = TRUE) +
  scale_color_manual(values = mod_cols) +
  labs(title = "intramodular connectivity") +
  NULL


module_entrez <- map(modules, ~bitr(.x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))

ggo_all <- map(module_entrez, ~enrichGO(gene = as.character(.x$ENTREZID), 
                                        OrgDb    = org.Hs.eg.db,
                                        ont      = "ALL",
                                        readable = TRUE))

plot_module <- function(module_go, module_name){
  dotplot(module_go) + 
    labs(title = module_name)
  
}

module_plots <- imap(ggo_all, plot_module)

MEList <- moduleEigengenes(
  expr = scWGCNA_out$exprMatrix,
  colors = labels2colors(dandc$dynamicMods),
  softPower = 2
)
MEs <- MEList$eigengenes

pdf("results/scwgcna/scwgcna_pt2_seu.pdf", width = 20)
scWGCNA_out$plots[[1]]
scWGCNA_out$plots[[2]]
intramod_plot
pdc <- plotDendroAndColors(dandc$geneTree,
                           labels2colors(dandc$dynamicMods),
                           "Dynamic Tree Cut",
                           dendroLabels = FALSE,
                           hang = 0.03,
                           addGuide = TRUE,
                           guideHang = 0.05,
                           main = "Gene dendrogram and module colors"
)
# plot module eigengenes
Heatmap(MEs, show_row_names = FALSE)
module_plots
dev.off()


