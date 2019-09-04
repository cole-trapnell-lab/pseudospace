
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

#### Load clean and processed data ####
MCF10A_pseudospace_cds <- readRDS("pseudospace_cds.rds")
HuMEC_pseudospace_cds <- readRDS("HuMEC_pseudospace_cds.rds")

# Load Mock and TGFB cds objects created in Figure1 code
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

Mock_inner_subset_cells <- row.names(pData(cds.list[["Mock"]])[pData(cds.list[["Mock"]])$pseudospace_quantile == 1,])
length(Mock_inner_subset_cells)
Mock_outer_subset_cells <- row.names(pData(cds.list[["Mock"]])[pData(cds.list[["Mock"]])$pseudospace_quantile == 10,])
length(Mock_outer_subset_cells)
TGFB_inner_subset_cells <- row.names(pData(cds.list[["TGFB"]])[pData(cds.list[["TGFB"]])$pseudospace_quantile == 1,])
length(TGFB_inner_subset_cells)
TGFB_outer_subset_cells <- row.names(pData(cds.list[["TGFB"]])[pData(cds.list[["TGFB"]])$pseudospace_quantile == 10,])
length(TGFB_outer_subset_cells)


MCF10A_pseudospace_cds <- MCF10A_pseudospace_cds[,c(Mock_inner_subset_cells,Mock_outer_subset_cells,
                                      TGFB_inner_subset_cells,TGFB_outer_subset_cells)]


HuMEC_cds.list <- readRDS("HuMEC_pseudospace_processed_trajectories_cds.list.rds")

# Identify MCF10A and HuMEC cells expressed at the beginning and end of spontaneous (mock) and TGFB-driven (tgfb) pseudospatial trajectories
HuMEC_mock_pseudospace_max <- max(pData(HuMEC_cds.list[["Mock"]])$Pseudotime)
HuMEC_mock_pseudospace_quantile <- cut(pData(HuMEC_cds.list[["Mock"]])$Pseudotime, 
                           breaks = seq(0, HuMEC_mock_pseudospace_max, HuMEC_mock_pseudospace_max/10), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5,6,7,8,9,10))
pData(HuMEC_cds.list[["Mock"]])$pseudospace_quantile <- HuMEC_mock_pseudospace_quantile

HuMEC_tgfb_pseudospace_max <- max(pData(HuMEC_cds.list[["TGFB"]])$Pseudotime)
HuMEC_tgfb_pseudospace_quantile <- cut(pData(HuMEC_cds.list[["TGFB"]])$Pseudotime, 
                           breaks = seq(0, HuMEC_tgfb_pseudospace_max, HuMEC_tgfb_pseudospace_max/10), 
                           include.lowest = TRUE, labels = c(1,2,3,4,5,6,7,8,9,10))
pData(HuMEC_cds.list[["TGFB"]])$pseudospace_quantile <- HuMEC_tgfb_pseudospace_quantile

HuMEC_Mock_inner_subset_cells <- row.names(pData(HuMEC_cds.list[["Mock"]])[pData(HuMEC_cds.list[["Mock"]])$pseudospace_quantile == 1,])
length(HuMEC_Mock_inner_subset_cells)
HuMEC_Mock_outer_subset_cells <- row.names(pData(HuMEC_cds.list[["Mock"]])[pData(HuMEC_cds.list[["Mock"]])$pseudospace_quantile == 10,])
length(HuMEC_Mock_outer_subset_cells)
HuMEC_TGFB_inner_subset_cells <- row.names(pData(HuMEC_cds.list[["TGFB"]])[pData(HuMEC_cds.list[["TGFB"]])$pseudospace_quantile == 1,])
length(HuMEC_TGFB_inner_subset_cells)
HuMEC_TGFB_outer_subset_cells <- row.names(pData(HuMEC_cds.list[["TGFB"]])[pData(HuMEC_cds.list[["TGFB"]])$pseudospace_quantile == 10,])
length(HuMEC_TGFB_outer_subset_cells)


HuMEC_pseudospace_cds <- HuMEC_pseudospace_cds[,c(HuMEC_Mock_inner_subset_cells,HuMEC_Mock_outer_subset_cells,
                                      HuMEC_TGFB_inner_subset_cells,HuMEC_TGFB_outer_subset_cells)]

MCF10A_pData <- pData(MCF10A_pseudospace_cds)
HuMEC_pData <- pData(HuMEC_pseudospace_cds)

new_fData <- fData(MCF10A_pseudospace_cds)

new_pData <- rbind(MCF10A_pData,HuMEC_pData)

MCF10A_exprs <- Biobase::exprs(MCF10A_pseudospace_cds)
HuMEC_exprs <- Biobase::exprs(HuMEC_pseudospace_cds)

new_exprs <- cbind(MCF10A_exprs,HuMEC_exprs)

new_fData <- new("AnnotatedDataFrame", data = new_fData)
new_pData <- new("AnnotatedDataFrame", data = new_pData)

HuMEC_MCF10A_cds <- newCellDataSet(as(as.matrix(new_exprs), "sparseMatrix"),
                       phenoData = new_pData,
                       featureData = new_fData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HuMEC_MCF10A_cds <- estimateSizeFactors(HuMEC_MCF10A_cds)

sample_id <- sapply(pData(HuMEC_MCF10A_cds)$cell, function(x){
    
    if(x %in% Mock_inner_subset_cells)return("MCF10A_Mock_early")
    if(x %in% Mock_outer_subset_cells)return("MCF10A_Mock_late")
    if(x %in% TGFB_inner_subset_cells)return("MCF10A_TGFB_early")
    if(x %in% TGFB_outer_subset_cells)return("MCF10A_TGFB_late")
    if(x %in% HuMEC_Mock_inner_subset_cells)return("HuMEC_Mock_early")
    if(x %in% HuMEC_Mock_outer_subset_cells)return("HuMEC_Mock_late")
    if(x %in% HuMEC_TGFB_inner_subset_cells)return("HuMEC_TGFB_early")
    if(x %in% HuMEC_TGFB_outer_subset_cells)return("HuMEC_TGFB_late")
    
})


pData(HuMEC_MCF10A_cds)$Cell.Type <- sample_id

# Score cells by expression of genes that define EMT and epithelial cells according to established geneset collections
source("loadGSCSafe.R")

hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")
GOGSC<-loadGSCSafe(file="Human_GO_bp_no_GO_iea_symbol.gmt")

EMT_marker_gene_names <- hallmarksGSC$gsc$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
EMT_marker_id <- unique(fData(HuMEC_MCF10A_cds)[fData(HuMEC_MCF10A_cds)$gene_short_name %in% EMT_marker_gene_names,]$id)

epithelial_marker_gene_names <- GOGSC$gsc$'EPIDERMIS DEVELOPMENT%GOBP%GO:0008544'
epithelial_marker_id <- unique(fData(HuMEC_MCF10A_cds)[fData(HuMEC_MCF10A_cds)$gene_short_name %in% epithelial_marker_gene_names,]$id)

Epithelial_marker_genes <- epithelial_marker_id
Canonical_EMT_marker_genes <- EMT_marker_id


EMT_signature_df <- as.data.frame(matrix(c(Epithelial_marker_genes,
                                                 Canonical_EMT_marker_genes,
                                                 rep("Epithelial_score",
                                                     length(Epithelial_marker_genes)),
                                                 rep("Canonical_EMT_score", 
                                                     length(Canonical_EMT_marker_genes))), 
                                               ncol = 2))
colnames(EMT_signature_df) <- c("id","signature")
EMT_signature_df$id <- as.character(EMT_signature_df$id)
EMT_signature_df$signature <- as.character(EMT_signature_df$signature)

HuMEC_MCF10A_EMT_scores <- calculate_signature_scores(HuMEC_MCF10A_cds,EMT_signature_df)

HuMEC_MCF10A_EMT_scores$Sample <- pData(HuMEC_MCF10A_cds)$Cell.Type

HuMEC_MCF10A_EMT_scores_scaled <- HuMEC_MCF10A_EMT_scores %>% mutate_at(c(1,2), funs(c(scale(.))))

HuMEC_MCF10A_EMT_scores_scaled$Sample <- factor(HuMEC_MCF10A_EMT_scores_scaled$Sample, 
                                   levels = c("MCF10A_Mock_early","MCF10A_Mock_late",
                                              "MCF10A_TGFB_early","MCF10A_TGFB_late",
                                              "HuMEC_Mock_early","HuMEC_Mock_late",
                                              "HuMEC_TGFB_early","HuMEC_TGFB_late"))

HuMEC_MCF10A_EMT_scores_scaled %>% 
filter(Sample %in% c("MCF10A_Mock_early","MCF10A_TGFB_late")) %>% 
group_by(Sample) %>% 
summarize(median_EMT_score = median(Canonical_EMT_score), median_Epithelial_score = median(Epithelial_score))

ggplot(HuMEC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Canonical_EMT_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = -1.479929, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = 1.163579, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24, family="Arial"),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
ylab("Hallmarks EMT\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HuMEC_MCF10A_Canonical_EMT_score.png", height = 6, width = 6)

ggplot(HuMEC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Epithelial_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = 0.7163635, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = -1.0615395, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24, family="Arial"),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
ylim(-3,5) +
ylab("GO-BP epithelial\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HuMEC_MCF10A_Epithelial_score.png", height = 6, width = 6)

compare_means(Canonical_EMT_score~Sample,
              data = HuMEC_MCF10A_EMT_scores_scaled[HuMEC_MCF10A_EMT_scores_scaled$Sample %in%
                                       c("MCF10A_Mock_early","MCF10A_Mock_late",
                                              "MCF10A_TGFB_early","MCF10A_TGFB_late"),], 
              method= "wilcox.test")

compare_means(Canonical_EMT_score~Sample,
              data = HuMEC_MCF10A_EMT_scores_scaled[HuMEC_MCF10A_EMT_scores_scaled$Sample %in%
                                       c("HuMEC_Mock_early","HuMEC_Mock_late",
                                        "HuMEC_TGFB_early","HuMEC_TGFB_late"),], 
              method= "wilcox.test")

compare_means(Epithelial_score~Sample,
              data = HuMEC_MCF10A_EMT_scores_scaled[HuMEC_MCF10A_EMT_scores_scaled$Sample %in%
                                       c("MCF10A_Mock_early","MCF10A_Mock_late",
                                              "MCF10A_TGFB_early","MCF10A_TGFB_late"),], 
              method= "wilcox.test")

compare_means(Epithelial_score~Sample,
              data = HuMEC_MCF10A_EMT_scores_scaled[HuMEC_MCF10A_EMT_scores_scaled$Sample %in%
                                       c("HuMEC_Mock_early","HuMEC_Mock_late",
                                        "HuMEC_TGFB_early","HuMEC_TGFB_late"),], 
              method= "wilcox.test")



