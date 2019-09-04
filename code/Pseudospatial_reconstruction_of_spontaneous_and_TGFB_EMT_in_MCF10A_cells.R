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
    library(monocle)
    })

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

#### Load data #####
pseudospace_cds <- readRDS("pseudospace_cds.rds")

# Create a cds subset for each stimulation condition that contains spatially isolated cells
cds.list <- list()

cds.list[["Mock"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "Mock"]
cds.list[["TGFB"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "TGFB"]

ggplot(pData(pseudospace_cds), aes(x = sample, y = total_umis)) + geom_violin() +
geom_jitter(alpha = 0.02, color = "brown4") +
stat_summary(fun.data=mean_sdl, geom="pointrange", color="black") +
ylim(0,40000) +
xlab("Sample") +
ylab("Number of UMIs") +
monocle:::monocle_theme_opts() +
theme(text = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust = 1)) +
ggsave("MCF10A_umis_per_sample.png", width = 5, height = 5)

pseudospace_cds <- detectGenes(pseudospace_cds, min_expr = 0.5)

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

for(sample in c("Mock","TGFB")){
    
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
for(sample in c("Mock","TGFB")){
    
    cds.list[[sample]] <- process_cds_Trajectory(cds.list[[sample]], Spatial.sig.genes.list[[sample]], ncenter = 500)
}

cds.list[["Mock"]] <- orderCells(cds.list[["Mock"]], reverse = T)
cds.list[["TGFB"]] <- orderCells(cds.list[["TGFB"]], reverse = T)

plot_cell_trajectory(cds.list[["Mock"]], color_by = "spatial_id",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=24), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#0075F2", "#D62828"), 
                   name = "Spatial Context") +
ggsave(file = "MCF10A_Mock_PseudospatialTrajectory.png", height = 6, width = 6)

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
  ggsave("Mock_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

ggplot(pData(cds.list[["Mock"]]), aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 24),
       axis.text.y = element_text(size = 14)) +
  ggsave("Mock_cell_density_accross_pseudospace_noFacet_geom_density.png", height = 3, width = 6)

plot_cell_trajectory(cds.list[["TGFB"]], color_by = "spatial_id",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=24), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#70163C", "#38726C"), 
                   name = "Spatial Context") +
ggsave(file = "MCF10A_TGFB_PseudospatialTrajectory.png", height = 6, width = 6)

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
  ggsave("TGFB_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

Mock_Figure1_Epi_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% c("CDH1","CRB3","DSP"))),]
Mock_Figure1_Mes_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% c("CDH2","FN1","VIM"))),]

plot_genes_violin(Mock_Figure1_Epi_Mar, grouping = "spatial_id", fill_by ="spatial_id", 
                  label_by_short_name = TRUE, ncol = 3, log_scale = TRUE, plot_trend = T, 
                  min_expr = 0.1) +                  
theme(axis.title.x=element_blank(), legend.position = "none", text=element_text(size=16)) + 
scale_fill_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave(file = "Mock_Figure1EpiMar_InnervsOuter.png", height = 3, width = 6)

plot_genes_violin(Mock_Figure1_Mes_Mar, grouping = "spatial_id", fill_by ="spatial_id", 
                  label_by_short_name = TRUE, ncol = 3, log_scale = TRUE, plot_trend = T, 
                  min_expr = 0.1) +                  
theme(axis.title.x=element_blank(), legend.position = "none", text=element_text(size=16)) + 
scale_fill_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave(file = "Mock_Figure1MesMar_InnervsOuter.png", height = 3, width = 6)

# Plot the expression of known EMT markers across pseudospace
Mock_Figure1_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]

plot_genes_in_pseudotime(Mock_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("MCF10A_Mock_Figure1Markers_byPseudospace..png", width = 6, height =5)

TGFB_Figure1_Mar <- cds.list[["TGFB"]][row.names(subset(fData(cds.list[["TGFB"]]), gene_short_name %in% c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]

plot_genes_in_pseudotime(TGFB_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#70163C","outer"="#38726C")) + 
ggsave("MCF10A_TGFB_Figure1Markers_byPseudospace..png", width = 6, height =5)

# Bin cells across pseudospace and determine the fraction of CDH1 single positive, CDH1/VIM double positive and VIM single positive cells upon confluence dependent EMT
mock_pseudospace_max <- max(pData(cds.list[["Mock"]])$Pseudotime)
mock_pseudospace_quantile <- cut(pData(cds.list[["Mock"]])$Pseudotime, 
                           breaks = seq(0, mock_pseudospace_max, mock_pseudospace_max/5), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5))
pData(cds.list[["Mock"]])$pseudospace_quantile <- mock_pseudospace_quantile

# Bin cells across pseudospace and determine the fraction of CDH1 single positive, CDH1/VIM double positive and VIM single positive cells upon confluence dependent EMT
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
  ggsave("Mock_cell_density_accross_pseudospace_CDH1_positive_cells_geom_density.png", 
    height = 3, width = 6)

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
  ggsave("Mock_cell_density_accross_pseudospace_CDH1_VIM_double_positive_cells_geom_density.png", 
    height = 3, width = 6)

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
  ggsave("Mock_cell_density_accross_pseudospace_VIM_positive_cells_geom_density.png", 
    height = 3, width = 6)

# Generate subsets from the ends of the mock and tgfb trajectories top determine whether 
# monocle always finds a continuously populated linear trajectory

mock_pseudospace_max <- max(pData(cds.list[["Mock"]])$Pseudotime)
mock_pseudospace_quantile <- cut(pData(cds.list[["Mock"]])$Pseudotime, 
                           breaks = seq(0, mock_pseudospace_max, mock_pseudospace_max/10), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5,6,7,8,9,10))
pData(cds.list[["Mock"]])$pseudospace_quantile <- mock_pseudospace_quantile

tgfb_pseudospace_max <- max(pData(cds.list[["TGFB"]])$Pseudotime)
tgfb_pseudospace_quantile <- cut(pData(cds.list[["TGFB"]])$Pseudotime, 
                           breaks = seq(0, tgfb_pseudospace_max, tgfb_pseudospace_max/10), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5,6,7,8,9,10))
pData(cds.list[["TGFB"]])$pseudospace_quantile <- tgfb_pseudospace_quantile

cds.list[["Mock subset"]] <- cds.list[["Mock"]][,pData(cds.list[["Mock"]])$pseudospace_quantile %in% c(1,10)]

cds.list[["Mock subset"]] <- preprocess_cds(cds.list[["Mock subset"]])


cds.list[["Mock subset"]] <- process_cds_Trajectory(cds.list[["Mock subset"]], Spatial.sig.genes.list[["Mock"]], ncenter = 500)

cds.list[["Mock subset"]] <- orderCells(cds.list[["Mock subset"]], reverse = F)

plot_cell_trajectory(cds.list[["Mock subset"]], color_by = "spatial_id",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=24), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#0075F2", "#D62828"), 
                   name = "Spatial Context") +
ggsave(file = "MCF10A_Mock_CDH1_VIM_subset_PseudospatialTrajectory.png", height = 6, width = 6)

ggplot(pData(cds.list[["Mock subset"]]), 
       aes(x = Pseudotime, fill = spatial_id, color = spatial_id)) + 
  geom_density() + 
  facet_wrap(~spatial_id, ncol = 1) +
  theme_classic() +
  scale_color_manual("Spatial Context", labels = c("inner colony", "outer colony"),
                     values = c("#000000","#000000")) +
  scale_fill_manual("Spatial Context", labels = c("inner colony", "outer colony")
                    , values = c("#0075F2","#D62828")) +
  xlab("Pseudospace") +
  monocle:::monocle_theme_opts() +
  theme(legend.position = "top", legend.direction = "vertical",  text=element_text(size=20)) +
ggsave("Mock_subset_cell_density_accross_pseudospace.png", height = 6, width = 5)

## Identify differentially expressed genes accross pseudospace
Pseudospatial.DEG.test.list <- list()

for(sample in c("Mock","TGFB")){
    
    Pseudospatial.DEG.test.list[[sample]] <- differentialGeneTest(cds.list[[sample]][expressed_genes.list[[sample]]],
                                                      cores = 1)
}

Pseudospatial.sig.genes.list <- list()

for(sample in names(Pseudospatial.DEG.test.list)){
    
    Pseudospatial.sig.genes.list[[sample]] <- row.names(subset(Pseudospatial.DEG.test.list[[sample]], qval < 1e-10))
    print(length(Pseudospatial.sig.genes.list[[sample]]))
    
}

# Calculate the area under the curve for the fitted splines of gene expression over pseudospace to focus on those with the highest fold change
Pseudospatial.auc.list <- list()

for(sample in "Mock"){
    
    Pseudospatial.auc.list[[sample]] <- auc_in_pseudospace(cds.list[[sample]][Pseudospatial.sig.genes.list[[sample]]],
                                                          cores = 20)
}

for(sample in names(Pseudospatial.auc.list)){
    
    Pseudospatial.auc.list[[sample]] <- Pseudospatial.auc.list[[sample]][,-7]
    
}

for(sample in names(Pseudospatial.auc.list)){
    
    Pseudospatial.auc.list[[sample]]$max_auc <- pmax(Pseudospatial.auc.list[[sample]]$'1',Pseudospatial.auc.list[[sample]]$'2',
                                                     Pseudospatial.auc.list[[sample]]$'3',Pseudospatial.auc.list[[sample]]$'4',
                                                     Pseudospatial.auc.list[[sample]]$'5')
    
}

# Identify a reasonable cutoff to filter DEGs by area under the curve
ggplot(Pseudospatial.auc.list[["Mock"]], aes(x = max_auc)) + 
geom_histogram() +
geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
xlim(0,100) +
monocle:::monocle_theme_opts() +
ggsave("X2.png")

Pseudospatial.sig.genes.high.auc.list <- list()

for(sample in c("Mock","TGFB")){
    
    Pseudospatial.sig.genes.high.auc.list[[sample]] <- Pseudospatial.auc.list[[sample]][Pseudospatial.auc.list[[sample]]$max_auc > 10,]$id

    print(length(Pseudospatial.sig.genes.high.auc.list[[sample]]))
    
}

Pseudospatial.sig.genes.high.auc.list[["Mock"]] <- Pseudospatial.auc.list[["Mock"]][Pseudospatial.auc.list[["Mock"]]$max_auc > 10,]$id

print(length(Pseudospatial.sig.genes.high.auc.list[["Mock"]]))

Pseudospace.cds.list <- list()

for(sample in c("Mock","TGFB")){
    
    Pseudospace.cds.list[[sample]] <- cds.list[[sample]][Pseudospatial.sig.genes.high.auc.list[[sample]],]
    
}

ann_colors <- list(Cluster = c("1"="brown3", "2"="darkgreen", "3"="darkblue", "4"="darkorchid4", "5"="firebrick1", "6"="goldenrod1"))

jpeg(file="MCF10A_Mock_PseudospaceHeatmap_q1e-10_aucRank10.jpeg", width=3, height=5, units="in",res=1000, quality=100)
Mock_ph <- plot_pseudotime_heatmap_binSize(Pseudospace.cds.list[["Mock"]], num_clusters = 6, 
  bin_size = 500, hmcols = magma(50),
  cores = 1, main = "Mock",return_heatmap=TRUE)
dev.off()

###### Gene Set Analysis of clustered genes that vary by spatial and pseudospatial context ######

# Load necessary functions
source("loadGSCSafe.R")
source('plot_gsea_go.R')
source("GSA_hyper_helper_functions.R")

## Load Gene Set Collections
GOGSC<-loadGSCSafe(file="Human_GO_bp_no_GO_iea_symbol.gmt")
hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")

Ensembl_GSAlist <- as.matrix(fData(cds.list[["Mock"]][expressed_genes.list[["Mock"]]])$gene_short_name)
rownames(Ensembl_GSAlist)<-row.names(fData(cds.list[["Mock"]][expressed_genes.list[["Mock"]]]))
colnames(Ensembl_GSAlist) <- c("gene_short_name")
Ensembl_GSAlist<-Ensembl_GSAlist[,1]
Ensembl_GSAlist<-toupper(Ensembl_GSAlist)
length(Ensembl_GSAlist)

mock.pseudo.gene.clusters <- cutree(Mock_ph$tree_row, 6)[row.names(fData(Pseudospace.cds.list[["Mock"]]))]

table(mock.pseudo.gene.clusters)

Mock_GOBP_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                    replace_gene_names_vec(mock.pseudo.gene.clusters
                                                                           ,Ensembl_GSAlist),
                                                    GOGSC)

Mock_hallmarks_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                    replace_gene_names_vec(mock.pseudo.gene.clusters
                                                                           ,Ensembl_GSAlist),
                                                    hallmarksGSC)

gsea_bar_plots(Mock_GOBP_GSAhyper, qval_cutoff = 0.0001, pattern = "%GOBP%GO:", width = 8, height = 10, sample = "Mock", gsc = "GOBP")


gsea_bar_plots(Mock_hallmarks_GSAhyper, qval_cutoff = 0.05, pattern = "HALLMARK_", width = 8, height = 10, sample = "Mock", gsc = "Hallmarks")

get_gsea_sig_results <- function(GSAhyper_list, qval_cutoff){
  
  GSA_hyper_results.list <- list()
  
  for(cluster in names(GSAhyper_list)){
    
    GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
    GSAhyper_df$gene_set <- row.names(GSAhyper_df)
    colnames(GSAhyper_df) <- c("qval","gene_set")
    
    GSA_hyper_results.list[[cluster]] <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
      mutate(gene_set = factor(gene_set, levels = gene_set))
    
  }
  
  return(GSA_hyper_results.list)
  
}

GOBP_sig_results <- get_gsea_sig_results(Mock_GOBP_GSAhyper, 0.05)
Hallmarks_sig_results <- get_gsea_sig_results(Mock_hallmarks_GSAhyper, 0.05)

### Rename clusters according to Figure 1 labels
names(GOBP_sig_results) <- c("5","1","6","4","2","3")
names(Hallmarks_sig_results) <- c("5","1","6","4","2","3")

GOBP_sig_results <- do.call("rbind",GOBP_sig_results)
GOBP_sig_results$cluster <- sapply(row.names(GOBP_sig_results), function(x){substring(x, 1, 1)})

Hallmarks_sig_results <- do.call("rbind",Hallmarks_sig_results)
Hallmarks_sig_results$cluster <- sapply(row.names(Hallmarks_sig_results), function(x){substring(x, 1, 1)})

# Plot expressionm of KRAS down-regulated regulated geneset from the cluster enriched for KRAS downregulation
mock_downregulated_cluster_2 <- names(mock.pseudo.gene.clusters[mock.pseudo.gene.clusters == 2])

mock_downregulated_cluster_2_genes <- fData(cds.list[["Mock"]])[mock_downregulated_cluster_2,]$gene_short_name

mock_KRAS_downregulated_cluster_2_genes <- intersect(hallmarksGSC$gsc$HALLMARK_KRAS_SIGNALING_DN,mock_downregulated_cluster_2_genes)

Mock_Figure1_KRAS_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% mock_KRAS_downregulated_cluster_2_genes)),]

plot_genes_in_pseudotime(Mock_Figure1_KRAS_Mar, color_by = "spatial_id", ncol = 3, min_expr = 0.1,) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20, family="Arial")) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("MCF10A_Mock_KRAS_dowregulated_genes_byPseudospace.png", width = 9, height =5)

# Plot expressionm of UPR geneset from the cluster enriched for UPR 
mock_upregulated_cluster_1 <- names(mock.pseudo.gene.clusters[mock.pseudo.gene.clusters == 1])


mock_upregulated_cluster_1_genes <- fData(cds.list[["Mock"]])[mock_upregulated_cluster_1,]$gene_short_name



mock_UPR_upregulated_cluster_1_genes <- intersect(hallmarksGSC$gsc$'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
                                                        mock_upregulated_cluster_1_genes)
mock_UPR_upregulated_cluster_1_genes <- mock_UPR_upregulated_cluster_1_genes[!(mock_UPR_upregulated_cluster_1_genes %in% "CKS1B")]



Mock_Figure1_UPR_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% mock_UPR_upregulated_cluster_1_genes)),]

plot_genes_in_pseudotime(Mock_Figure1_UPR_Mar, color_by = "spatial_id", ncol = 3, min_expr = 0.1,) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20, family="Arial")) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("MCF10A_Mock_UPR_upregulated_genes_byPseudospace.png", width = 9, height = 5)


# Plot expression of UPR geneset from the cluster enriched for EMT 
mock_upregulated_cluster_3 <- names(mock.pseudo.gene.clusters[mock.pseudo.gene.clusters == 3])


mock_upregulated_cluster_3_genes <- fData(cds.list[["Mock"]])[mock_upregulated_cluster_3,]$gene_short_name



mock_UPR_upregulated_cluster_3_genes <- intersect(hallmarksGSC$gsc$'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
                                                        mock_upregulated_cluster_3_genes)



Mock_Figure1_UPR_EMTClusterMar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% mock_UPR_upregulated_cluster_3_genes)),]

plot_genes_in_pseudotime(Mock_Figure1_UPR_EMTClusterMar, color_by = "spatial_id", ncol = 3, min_expr = 0.1,) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20, family="Arial")) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("MCF10A_Mock_UPR_EMTCluster_upregulated_genes_byPseudospace.png", width = 9, height =3.5)

mock_Epidermis_downregulated_cluster_2_genes <- intersect(GOGSC$gsc$`EPIDERMIS DEVELOPMENT%GOBP%GO:0008544`,
                                                        mock_downregulated_cluster_2_genes)

### Create tSNE embeddings of Mock and TGFB exposed cells to examine the dynamics of EMT across clusters

Mock_cds <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "Mock"]

Mock_cds <- preprocess_cds(Mock_cds)

plot_pc_variance_explained(Mock_cds[expressed_genes.list[["Mock"]]], norm_method = "log", pseudo_expr = 1, 
                           return_all = F, max_components = 100) + 
ggsave("Mock_PCA.png")

Mock_cds <- reduceDimension(Mock_cds[expressed_genes.list[["Mock"]]], max_components=2, 
                               norm_method = 'log', num_dim = 25, reduction_method = 'tSNE',
                               verbose = T, cores = 1)

Mock_cds <- clusterCells(Mock_cds, verbose = T)

plot_rho_delta(Mock_cds, rho_threshold = 60, delta_threshold = 5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

Mock_cds <- clusterCells(Mock_cds, verbose = T, 
                                 rho_threshold = 50, 
                                 delta_threshold = 5,
                                 skip_rho_sigma = T)

mock_CDH1_subset <- Mock_cds[fData(Mock_cds)$gene_short_name == "CDH1",]
mock_CDH1_expr <- as.matrix(t(exprs(mock_CDH1_subset))/pData(mock_CDH1_subset)$Size_Factor)
#mock_CDH1_expr <- log10(mock_CDH1_expr + 0.1)

mock_DSP_subset <- Mock_cds[fData(Mock_cds)$gene_short_name == "DSP",]
mock_DSP_expr <- as.matrix(t(exprs(mock_DSP_subset))/pData(mock_DSP_subset)$Size_Factor)
#mock_DSP_expr <- log10(mock_DSP_expr + 0.1)

mock_VIM_subset <- Mock_cds[fData(Mock_cds)$gene_short_name == "VIM",]
mock_VIM_expr <- as.matrix(t(exprs(mock_VIM_subset))/pData(mock_VIM_subset)$Size_Factor)
#mock_VIM_expr <- log10(mock_VIM_expr + 0.1)

pData(Mock_cds)$CDH1_expr <- mock_CDH1_expr
pData(Mock_cds)$DSP_expr <- mock_DSP_expr
pData(Mock_cds)$VIM_expr <- mock_VIM_expr

plot_cell_clusters(Mock_cds, color_by = "log10(CDH1_expr + 0.1)", cell_size = 0.5) +
scale_color_viridis(option = "magma") +
ggtitle("CDH1") +
 theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 20),
              legend.title=element_blank(),
               axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 10),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("spontaneous_EMT_CDH1_tSNE.png", width = 4, height = 4)

plot_cell_clusters(Mock_cds, color_by = "log10(VIM_expr + 0.1)", cell_size = 0.5)+
scale_color_viridis(option = "magma") +
ggtitle("VIM") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("spontaneous_EMT_VIM_tSNE.png", width = 4, height = 4.25)

plot_cell_clusters(Mock_cds, color_by = "log10(DSP_expr + 0.1)", cell_size = 0.5)+
scale_color_viridis(option = "magma") +
ggtitle("DSP") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("spontaneous_EMT_DSP_tSNE.png", width = 4, height = 4.25)

plot_cell_clusters(Mock_cds, color_by = "spatial_id", cell_size = 0.5)+
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828"), labels = c("inner" = "inner\ncolony", "outer" = "outer\ncolony")) + 
ggtitle("Spatial Context") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("spontaneous_EMT_spatial_context_tSNE.png", width = 4, height = 4.25)


TGFB_cds <- pseudospace_cds[,pData(pseudospace_cds)$treatment_id == "TGFB"]

TGFB_cds <- preprocess_cds(TGFB_cds)

plot_pc_variance_explained(TGFB_cds[expressed_genes.list[["TGFB"]]], norm_method = "log", pseudo_expr = 1, 
                           return_all = F, max_components = 100) + 
ggsave("TGFB_PCA.png")

TGFB_cds <- reduceDimension(TGFB_cds[expressed_genes.list[["TGFB"]]], max_components=2, 
                               norm_method = 'log', num_dim = 25, reduction_method = 'tSNE',
                               verbose = T, cores = 1)

TGFB_cds <- clusterCells(TGFB_cds, verbose = T)

plot_rho_delta(TGFB_cds, rho_threshold = 60, delta_threshold = 5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

TGFB_cds <- clusterCells(TGFB_cds, verbose = T, 
                                 rho_threshold = 50, 
                                 delta_threshold = 5,
                                 skip_rho_sigma = T)

tgfb_VIM_subset <- TGFB_cds[fData(TGFB_cds)$gene_short_name == "VIM",]
tgfb_VIM_expr <- as.matrix(t(exprs(tgfb_VIM_subset))/pData(tgfb_VIM_subset)$Size_Factor)

tgfb_FN1_subset <- TGFB_cds[fData(TGFB_cds)$gene_short_name == "FN1",]
tgfb_FN1_expr <- as.matrix(t(exprs(tgfb_FN1_subset))/pData(tgfb_FN1_subset)$Size_Factor)

pData(TGFB_cds)$VIM_expr <- tgfb_VIM_expr
pData(TGFB_cds)$FN1_expr <- tgfb_FN1_expr

plot_cell_clusters(TGFB_cds, color_by = "log10(VIM_expr + 0.1)", cell_size = 0.5)+
scale_color_viridis(option = "magma") +
ggtitle("VIM") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("TGFB_driven_EMT_VIM_tSNE.png", width = 4, height = 4.25)

plot_cell_clusters(TGFB_cds, color_by = "log10(FN1_expr + 0.1)", cell_size = 0.5)+
scale_color_viridis(option = "magma") +
ggtitle("FN1") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("TGFB_driven_EMT_FN1_tSNE.png", width = 4, height = 4.25)

plot_cell_clusters(TGFB_cds, color_by = "spatial_id", cell_size = 0.5)+
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828"), labels = c("inner" = "inner\ncolony", "outer" = "outer\ncolony")) + 
ggtitle("Spatial Context") +
theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5, size = 22),
              legend.title=element_blank(),
              axis.title.x = element_text(size = 20),
               axis.title.y = element_text(size = 20),
              legend.text=element_text(size = 12),
              legend.key.width = unit(0.25, "cm"),
             legend.key.height = unit(0.75, "cm")) +
ggsave("TGFB_EMT_spatial_context_tSNE.png", width = 4, height = 4.25)

saveRDS(cds.list,"pseudospace_processed_trajectories_cds.list.rds")

