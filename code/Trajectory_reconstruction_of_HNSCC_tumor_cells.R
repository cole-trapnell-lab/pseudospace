
###### Load packages ######
# Load necessary packages for single cell RNA-Seq analysis including packages for downstream Gene Ontology Analysis
suppressPackageStartupMessages({
library(devtools)
library(stringr)
library(scales)
library(dtw)
library(monocle)
library(reshape2)
library(GSA)
library(limma)
library(DBI)
library(MASS)
library(plyr)
library(dplyr)
library(matrixStats)
library(piano)
library(cluster)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(viridis)})

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

getDTW <- function(query_cds, ref_cds, ref, query, cores) 
{
    alignment_genes <- row.names(subset(fData(ref_cds),use_for_ordering))
    
    ref_align_cds <- ref_cds[alignment_genes,]
    query_align_cds <- query_cds[alignment_genes,]
    pData(ref_align_cds)$cell_id <- row.names(pData(ref_align_cds))
    pData(ref_align_cds)$Pseudotime <- 100 * pData(ref_align_cds)$Pseudotime/max(pData(ref_align_cds)$Pseudotime)
    ref_align_cds <- ref_align_cds[alignment_genes, as.character(arrange(pData(ref_align_cds), 
        Pseudotime)$cell_id)]
    pData(query_align_cds)$cell_id <- row.names(pData(query_align_cds))
    pData(query_align_cds)$Pseudotime <- 100 * pData(query_align_cds)$Pseudotime/max(pData(query_align_cds)$Pseudotime)
    query_align_cds <- query_align_cds[alignment_genes, as.character(arrange(pData(query_align_cds), 
        Pseudotime)$cell_id)]
    smoothed_ref_exprs <- genSmoothCurves(ref_align_cds[alignment_genes], 
        data.frame(Pseudotime = seq(0, 100, by = 1)), cores = cores)
    smoothed_ref_exprs <- smoothed_ref_exprs[rowSums(is.na(smoothed_ref_exprs)) == 
        0, ]
    vst_smoothed_ref_exprs <- vstExprs(ref_cds, expr_matrix = smoothed_ref_exprs)
    smoothed_query_exprs <- genSmoothCurves(query_align_cds[alignment_genes], 
        data.frame(Pseudotime = seq(0, 100, by = 1)), cores = cores)
    smoothed_query_exprs <- smoothed_query_exprs[rowSums(is.na(smoothed_query_exprs)) == 
        0, ]
    vst_smoothed_query_exprs <- vstExprs(query_cds, expr_matrix = smoothed_query_exprs)
    alignment_genes <- intersect(row.names(vst_smoothed_ref_exprs), 
        row.names(vst_smoothed_query_exprs))
    ref_matrix <- t(scale(t(vst_smoothed_ref_exprs[alignment_genes, 
        ])))
    query_matrix <- t(scale(t(vst_smoothed_query_exprs[alignment_genes, 
        ])))

    ref_query_dtw <- align_cells(query_matrix, ref_matrix, step_pattern = rabinerJuangStepPattern(3, 
        "c"), open.begin = T, open.end = T)
    
    return(ref_query_dtw)
}

#### Load data ####
HNSCC_cds <- readRDS("HSNCC_cds.rds")

expressed_genes <- row.names(fData(HNSCC_cds)[Matrix::rowSums(Biobase::exprs(HNSCC_cds) > 0) > 50 ,])

# Filter HNSCC data for samples with more than 100 cells, non-lymph node, cancer cells processesd with Maxima enzyme
# Exclude HNSCC17 which gives a super funky trajectory that we can't compare to MCF10A. Would be cool to
# to see if there are sub-clones 

metadata_summary <- pData(HNSCC_cds) %>% 
filter(clasified_as_cancer_cell == "1" & Maxima_enzyme == "0" & lymph_node == "0") %>% 
group_by(patient_id) %>% summarize(n = n()) %>% arrange(desc(n))

metadata_summary

HNSCC_cds.list <- list()

HNSCC_cds.list[["HNSCC20"]] <- HNSCC_cds[,pData(HNSCC_cds)$patient_id == "HNSCC20" & 
                                         pData(HNSCC_cds)$clasified_as_cancer_cell == "1" & 
                                         pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                                         pData(HNSCC_cds)$lymph_node == "0"]
HNSCC_cds.list[["HNSCC18"]] <- HNSCC_cds[,pData(HNSCC_cds)$patient_id == "HNSCC18" & 
                                         pData(HNSCC_cds)$clasified_as_cancer_cell == "1" & 
                                         pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                                         pData(HNSCC_cds)$lymph_node == "0"]
HNSCC_cds.list[["HNSCC22"]] <- HNSCC_cds[,pData(HNSCC_cds)$patient_id == "HNSCC22"  & 
                                         pData(HNSCC_cds)$clasified_as_cancer_cell == "1" & 
                                         pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                                         pData(HNSCC_cds)$lymph_node == "0"]

pData(HNSCC_cds.list[["HNSCC20"]])$cell <- row.names(pData(HNSCC_cds.list[["HNSCC20"]]))
pData(HNSCC_cds.list[["HNSCC22"]])$cell <- row.names(pData(HNSCC_cds.list[["HNSCC22"]]))
pData(HNSCC_cds.list[["HNSCC18"]])$cell <- row.names(pData(HNSCC_cds.list[["HNSCC18"]]))

for(patient in names(HNSCC_cds.list)){
    
names(pData(HNSCC_cds.list[[patient]])$Maxima_enzyme) <- NULL
names(pData(HNSCC_cds.list[[patient]])$lymph_node) <- NULL
names(pData(HNSCC_cds.list[[patient]])$clasified_as_cancer_cell) <- NULL
names(pData(HNSCC_cds.list[[patient]])$clasified_as_non_cancer_cell) <- NULL
names(pData(HNSCC_cds.list[[patient]])$non_cancer_cell_type) <- NULL
    
}

for(patient in names(HNSCC_cds.list)){
  
  print(dim(HNSCC_cds.list[[patient]]))
  
}

for(patient in names(HNSCC_cds.list)){

  HNSCC_cds.list[[patient]] <- preprocess_cds(HNSCC_cds.list[[patient]])

}

for(patient in names(HNSCC_cds.list)){
  
  HNSCC_cds.list[[patient]] <- reduceDimension(HNSCC_cds.list[[patient]][expressed_genes,],
                                               max_components=2, 
                                               norm_method = 'log', num_dim = 25, reduction_method = 'tSNE',
                                               verbose = T, cores = detectCores()-1)
  
}

for(patient in names(HNSCC_cds.list)){
  
  HNSCC_cds.list[[patient]] <- setOrderingFilter(HNSCC_cds.list[[patient]], expressed_genes)
  HNSCC_cds.list[[patient]] <- reduceDimension(HNSCC_cds.list[[patient]], norm_method = "log")
  
}

HNSCC_cds.list[["HNSCC20"]] <- orderCells(HNSCC_cds.list[["HNSCC20"]])
HNSCC_cds.list[["HNSCC22"]] <- orderCells(HNSCC_cds.list[["HNSCC22"]])
HNSCC_cds.list[["HNSCC18"]] <- orderCells(HNSCC_cds.list[["HNSCC18"]], reverse = FALSE)

plot_cell_trajectory(HNSCC_cds.list[["HNSCC20"]], color_by = "Pseudotime", show_branch_points = FALSE) +
theme(legend.position="right", 
      text=element_text(size=24, family="Arial"), 
      legend.direction = "vertical",
      legend.title = element_text(size = 14),
     legend.text = element_text(size = 11),
     legend.key.width = unit(0.2, "in"),
     legend.key.height = unit(0.4, "in")) +
scale_color_viridis(option = "magma") +
ggsave(file = "HNSCC20_Trajectory.png", height = 4, width = 5)

plot_genes_in_pseudotime(HNSCC_cds.list[["HNSCC20"]][fData(HNSCC_cds.list[["HNSCC22"]])$gene_short_name == "SNAI2",],
                        min_expr = 1, color_by = "Pseudotime") +
scale_color_viridis(option = "magma") + 
theme(text = element_text(size = 24, family = "Arial")) +
ggsave("SNAI2_expression_HNSCC20_pseudotime.png", height = 4, width = 6)

plot_cell_trajectory(HNSCC_cds.list[["HNSCC22"]], color_by = "Pseudotime", show_branch_points = FALSE) +
theme(legend.position="right", 
      text=element_text(size=24, family="Arial"), 
      legend.direction = "vertical",
      legend.title = element_text(size = 14),
     legend.text = element_text(size = 11),
     legend.key.width = unit(0.2, "in"),
     legend.key.height = unit(0.4, "in")) +
scale_color_viridis(option = "magma") +
ggsave(file = "HNSCC22_Trajectory.png", height = 4, width = 5)

plot_genes_in_pseudotime(HNSCC_cds.list[["HNSCC22"]][fData(HNSCC_cds.list[["HNSCC22"]])$gene_short_name == "SNAI2",],
                        min_expr = 1, color_by = "Pseudotime") +
scale_color_viridis(option = "magma") +
theme(text = element_text(size = 24, family = "Arial")) +
ggsave("SNAI2_expression_HNSCC22_pseudotime.png", height = 4, width = 6)

plot_cell_trajectory(HNSCC_cds.list[["HNSCC18"]], color_by = "Pseudotime") +
theme(legend.position="right", 
      text=element_text(size=24, family="Arial"), 
      legend.direction = "vertical",
      legend.title = element_text(size = 14),
     legend.text = element_text(size = 11),
     legend.key.width = unit(0.2, "in"),
     legend.key.height = unit(0.4, "in")) +
scale_color_viridis(option = "magma") +
ggsave(file = "HNSCC18_Trajectory.png", height = 4, width = 5)

plot_genes_in_pseudotime(HNSCC_cds.list[["HNSCC18"]][fData(HNSCC_cds.list[["HNSCC18"]])$gene_short_name == "SNAI2",],
                        min_expr = 1, color_by = "Pseudotime") +
scale_color_viridis(option = "magma") +
theme(text = element_text(size = 24, family = "Arial")) +
ggsave("SNAI2_expression_HNSCC18_pseudotime.png", height = 4, width = 6)

# Load Mock and TGFB cds objects created in Figure1 code to align MF10A and HNSCC trajectories
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

HNSCC_genes <- unique(intersect(expressed_genes,row.names(fData(HNSCC_cds))))

cds.list[["Mock"]] <- cds.list[["Mock"]][HNSCC_genes,]
cds.list[["TGFB"]] <- cds.list[["TGFB"]][HNSCC_genes,]


HNSCC20.to.Mock.dtw <- getDTW(HNSCC_cds.list[["HNSCC20"]], cds.list[["Mock"]],
                                ref = "Mock", query = "HNSCC20",
                                cores = 1)

HNSCC20.to.TGFB.dtw <- getDTW(HNSCC_cds.list[["HNSCC20"]], cds.list[["TGFB"]],
                                ref = "TGFB", query = "HNSCC20",
                                cores = 1)

HNSCC22.to.Mock.dtw <- getDTW(HNSCC_cds.list[["HNSCC22"]], cds.list[["Mock"]],
                                ref = "Mock", query = "HNSCC22",
                                cores = 1)

HNSCC22.to.TGFB.dtw <- getDTW(HNSCC_cds.list[["HNSCC22"]], cds.list[["TGFB"]],
                                ref = "TGFB", query = "HNSCC22",
                                cores = 1)

HNSCC18.to.Mock.dtw <- getDTW(HNSCC_cds.list[["HNSCC18"]], cds.list[["Mock"]],
                                ref = "Mock", query = "HNSCC18",
                                cores = 1)

HNSCC18.to.TGFB.dtw <- getDTW(HNSCC_cds.list[["HNSCC18"]], cds.list[["TGFB"]],
                                ref = "TGFB", query = "HNSCC18",
                                cores = 1)

pdf("HNSCC20.to.Mock.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC20.to.Mock.dtw, normalize=T, xlab="Spontaneous EMT", ylab="HNSCC20 Pseudotime")
dev.off()

pdf("HNSCC22.to.Mock.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC22.to.Mock.dtw, normalize=T, xlab="Spontaneous EMT", ylab="HNSCC22 Pseudotime")
dev.off()

pdf("HNSCC18.to.Mock.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC18.to.Mock.dtw, normalize=T, xlab="Spontaneous EMT", ylab="HNSCC18 Pseudotime")
dev.off()

pdf("HNSCC20.to.TGFB.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC20.to.TGFB.dtw, normalize=T, xlab="TGFB EMT", ylab="HNSCC20 Pseudotime")
dev.off()

pdf("HNSCC22.to.TGFB.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC22.to.TGFB.dtw, normalize=T, xlab="TGFB EMT", ylab="HNSCC22 Pseudotime")
dev.off()


pdf("HNSCC18.to.TGFB.dtw.openEnd.openStart.pdf")
dtwPlotDensity(HNSCC18.to.TGFB.dtw, normalize=T, xlab="TGFB EMT", ylab="HNSCC18 Pseudotime")
dev.off()


