###### Load packages ######
# Load necessary packages for single cell RNA-Seq analysis including packages for downstream Gene Ontology Analysis
suppressPackageStartupMessages({
    library(devtools)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(reshape2)
    library(pheatmap)
    library(piano)
    library(colorRamps)
    library(viridis)
    library(stringr)
    library(SFS)
    library(dtw)
    library(monocle)
    library(ggpubr)
    })

sessionInfo()

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

#### Load data ####
pseudospace_cds <- readRDS("HuMEC_pseudospace_cds.rds")

# Create a cds subset for each stimulation condition that contains spatially isolated cells
cds.list <- list()

cds.list[["Mock"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "Mock"]
cds.list[["TGFB"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "TGFB"]

dim(cds.list[["Mock"]])
dim(cds.list[["TGFB"]])

for(sample in names(cds.list)){
    
    print(pData(cds.list[[sample]]) %>%
         group_by(sample) %>% summarize(n = n()))
    
}

# Identify genes that are expressed in at least 50 of cells
expressed_genes.list <- list()

expressed_genes.list[["Mock"]] <- row.names(fData(cds.list[["Mock"]])[Matrix::rowSums(Biobase::exprs(cds.list[["Mock"]]) > 0) > 50 ,])
length(expressed_genes.list[["Mock"]])
expressed_genes.list[["TGFB"]] <- row.names(fData(cds.list[["TGFB"]])[Matrix::rowSums(Biobase::exprs(cds.list[["TGFB"]]) > 0) > 50 ,])
length(expressed_genes.list[["TGFB"]])

for(sample in names(cds.list)){
    
    cds.list[[sample]] <- preprocess_cds(cds.list[[sample]])
    
}

# Identify genes that vary significantly between inner and outer cells fractions
Spatial.DEG.test.list <- list()

for(sample in names(cds.list)){
    
    Spatial.DEG.test.list[[sample]] <- differentialGeneTest(cds.list[[sample]][expressed_genes.list[[sample]]],
                                                      fullModelFormulaStr = "~spatial_id",
                                                      reducedModelFormulaStr = "~1",
                                                      cores = 1)
    print(paste0("done with ",sample))
}

# Calculate fold change in expression levels of significant genes between cell fractions isolated by space
for(sample in names(Spatial.DEG.test.list)){
    
    diff_test_genes <- row.names(Spatial.DEG.test.list[[sample]])
    diff_cds <- cds.list[[sample]][diff_test_genes]
    diff_FC <- diff_foldChange(diff_cds, "spatial_id","inner")
    Spatial.DEG.test.list[[sample]]$log2_foldChange <- diff_FC$log2FC_outer
    
    rm(diff_test_genes,diff_cds,diff_FC)
}

Spatial.sig.genes.list <- list()

for(sample in names(Spatial.DEG.test.list)){
    
    Spatial.sig.genes.list[[sample]] <- row.names(subset(Spatial.DEG.test.list[[sample]], qval <= 1e-6 &
                                                        abs(log2_foldChange) >= 1))
    print(length(Spatial.sig.genes.list[[sample]]))
}

# Create pseudospatial trajectories and examine the distribution of inner and outer cells within them
for(sample in names(cds.list)){
    
    cds.list[[sample]] <- process_cds_Trajectory(cds.list[[sample]], Spatial.sig.genes.list[[sample]], ncenter = 500)
}

cds.list[["Mock"]] <- orderCells(cds.list[["Mock"]], reverse = T)
cds.list[["TGFB"]] <- orderCells(cds.list[["TGFB"]], reverse = T)

plot_cell_trajectory(cds.list[["Mock"]], color_by = "spatial_id",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=24), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#0075F2", "#D62828"), 
                   name = "Spatial Context") +
ggsave(file = "HuMEC_Mock_PseudospatialTrajectory.png", height = 6, width = 6)

ggplot(pData(cds.list[["Mock"]]), aes(x = Pseudotime, fill = spatial_id, color = spatial_id)) + 
  geom_density() + 
  facet_wrap(~spatial_id, ncol = 1) +
  theme_classic() +
  scale_color_manual("Spatial Context", labels = c("inner colony", "outer colony"),
                     values = c("#000000","#000000")) +
  scale_fill_manual("Spatial Context", labels = c("inner colony", "outer colony")
                    , values = c("#0075F2","#D62828")) +
  xlab("Pseudospace") +
  ylab("Cell density") +
  monocle:::monocle_theme_opts() +
  theme(legend.position = "top", legend.direction = "vertical",  text=element_text(size=20)) +
  ggsave("HuMEC_Mock_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

plot_cell_trajectory(cds.list[["TGFB"]], color_by = "spatial_id",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=24), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#70163C","#38726C"), 
                   name = "Spatial Context") +
ggsave(file = "HuMEC_TGFB_PseudospatialTrajectory.png", height = 6, width = 6)

ggplot(pData(cds.list[["TGFB"]]), aes(x = Pseudotime, fill = spatial_id, color = spatial_id)) + 
  geom_density() + 
  facet_wrap(~spatial_id, ncol = 1) +
  theme_classic() +
  scale_color_manual("Spatial Context", labels = c("inner colony", "outer colony"),
                     values = c("#000000","#000000")) +
  scale_fill_manual("Spatial Context", labels = c("inner colony", "outer colony")
                    , values = c("#70163C", "#38726C")) +
  xlab("Pseudospace") +
  ylab("Cell density") +
  monocle:::monocle_theme_opts() +
  theme(legend.position = "top", legend.direction = "vertical",  text=element_text(size=20)) +
  ggsave("HuMEC_TGFB_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

# Plot the expression of known EMT markers across pseudospace
Mock_Figure1_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]

plot_genes_in_pseudotime(Mock_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("HuMEC_Mock_Figure1Markers_byPseudospace.png", width = 6, height =5)

TGFB_Figure1_Mar <- cds.list[["TGFB"]][row.names(subset(fData(cds.list[["TGFB"]]), gene_short_name %in% c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]

plot_genes_in_pseudotime(TGFB_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#70163C","outer"="#38726C")) + 
ggsave("HuMEC_TGFB_Figure1Markers_byPseudospace.png", width = 6, height =5)

# Bin cells across pseudospace and determine the fraction of CDH1 single positive, CDH1/VIM double positive and VIM single positive cells upon confluence dependent EMT
mock_pseudospace_max <- max(pData(cds.list[["Mock"]])$Pseudotime)
mock_pseudospace_quantile <- cut(pData(cds.list[["Mock"]])$Pseudotime, 
                           breaks = seq(0, mock_pseudospace_max, mock_pseudospace_max/5), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5))
pData(cds.list[["Mock"]])$pseudospace_quantile <- mock_pseudospace_quantile

tgfb_pseudospace_max <- max(pData(cds.list[["TGFB"]])$Pseudotime)
tgfb_pseudospace_quantile <- cut(pData(cds.list[["TGFB"]])$Pseudotime, 
                           breaks = seq(0, tgfb_pseudospace_max, tgfb_pseudospace_max/5), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5))
pData(cds.list[["TGFB"]])$pseudospace_quantile <- tgfb_pseudospace_quantile

CDH1_VIM_cds_subset <- cds.list[["Mock"]][fData(cds.list[["Mock"]])$gene_short_name %in% c("CDH1","VIM"),!is.na(pData(cds.list[["Mock"]])$pseudospace_quantile)]
CDH1_VIM_cds_exprs <- exprs(CDH1_VIM_cds_subset)
CDH1_VIM_cds_exprs <- Matrix::t(Matrix::t(CDH1_VIM_cds_exprs)/sizeFactors(CDH1_VIM_cds_subset))

CDH1_expression_cutoff <- mean(exprs(CDH1_VIM_cds_subset[fData(CDH1_VIM_cds_subset)$gene_short_name == "CDH1",]))
VIM_expression_cutoff <- mean(exprs(CDH1_VIM_cds_subset[fData(CDH1_VIM_cds_subset)$gene_short_name == "VIM",]))

CDH1_expression_cutoff
VIM_expression_cutoff

CDH1_VIM_double_positive_cells <- colnames(CDH1_VIM_cds_exprs[,CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff])
CDH1_positive_cells <- colnames(CDH1_VIM_cds_exprs[,CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff & !(colnames(CDH1_VIM_cds_exprs) %in% CDH1_VIM_double_positive_cells)])
VIM_positive_cells <- colnames(CDH1_VIM_cds_exprs[,CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & !(colnames(CDH1_VIM_cds_exprs) %in% CDH1_VIM_double_positive_cells)])

mock_pData <- pData(cds.list[["Mock"]][,!is.na(pData(cds.list[["Mock"]])$pseudospace_quantile)])
mock_pData$positive_marker <- sapply(mock_pData$cell, function(x){
    
    if(x %in% CDH1_VIM_double_positive_cells){
        return("CDH1/VIM double positive")
    }
    
    if(x %in% CDH1_positive_cells){
        return("CDH1 single positive")
    }
    
    if(x %in% VIM_positive_cells){
        return("VIM single positive")
    }
    
    return(NA)
    
})

ggplot(pData(cds.list[["Mock"]])[CDH1_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1 single positive cells") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("HuMEC_Mock_cell_density_accross_pseudospace_CDH1_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.list[["Mock"]])[CDH1_VIM_double_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1/VIM double positive cells") +
  theme(legend.position = "none",  
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("HuMEC_Mock_cell_density_accross_pseudospace_CDH1_VIM_double_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.list[["Mock"]])[VIM_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("VIM single positive cells") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("HuMEC_Mock_cell_density_accross_pseudospace_VIM_positive_cells_geom_density.png", height = 3, width = 6)

### Add CDH2 manually for marker plotting
expressed_genes <- unique(c(union(expressed_genes.list[["Mock"]], expressed_genes.list[["TGFB"]]), "ENSG00000170558"))
length(expressed_genes)

# Use dynamic time warping to align Mock and TGFB pseudospatial trajectories and create a cds object of aligned trajectories
TGFB.to.Mock.aligned.cds <- getDTWcds(cds.list[["TGFB"]],cds.list[["Mock"]], 
                                      ref = "Mock", query = "TGFB", 
                                      expressed_genes = expressed_genes, cores = 10)


cds.aligned.list <- list()
cds.aligned.list[["TGFB to Mock"]] <- TGFB.to.Mock.aligned.cds

for(alignment in names(cds.aligned.list)){
    
    cds.aligned.list[[alignment]] <- preprocess_cds(cds.aligned.list[[alignment]])
    
}

# Examine EMT marker expression across aligned Pseudospace
compare_cell_types_in_pseudospace(cds.aligned.list[["TGFB to Mock"]][row.names(subset(fData(cds.aligned.list[["TGFB to Mock"]]), 
                                                                                    gene_short_name %in% c("CDH1","CDH2","CRB3","FN1","DSP","VIM"))),], 
                                  color_by="Cell.Type", df=3, min_expr=0.1, cell_alpha = 0.01, line_size = 1, ncol = 2,
                                 panel_order = c("CDH1","CDH2","CRB3","FN1","DSP", "VIM")) + 
    scale_color_manual(values = c("Mock" = "#0075F2", "TGFB" = "#70163C"), name = "Treatment") +
    theme(legend.position="none", text=element_text(size=20, family="Arial")) +
    ggsave("HuMEC_EMT_markers_Mock_TGFB_aligned_pseudospace.png", width = 6, height = 5)

saveRDS(cds.list, "HuMEC_pseudospace_processed_trajectories_cds.list.rds")

