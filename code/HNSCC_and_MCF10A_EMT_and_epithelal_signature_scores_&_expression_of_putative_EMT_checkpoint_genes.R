###### Load packages ######
# Load necessary packages for single cell RNA-Seq analysis including packages for downstream Gene Ontology Analysis
suppressPackageStartupMessages({
library(devtools)
library(stringr)
library(scales)
library(dtw)
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
library(viridis)
library(dendextend)
library(monocle)})

sessionInfo()

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

calculate_signature_scores <- function (cds, gene_list) 
{
    vst = as.data.frame(as.matrix(log2((Matrix::t(Matrix::t(Biobase:::exprs(cds))/pData(cds)$Size_Factor)) + 1)))

    vst = vst[gene_list$id, ]
    
    vst = merge(vst, gene_list[, names(gene_list) %in% c("id", 
        "signature")], by.x = "row.names", by.y = "id")
    
    vst = melt(vst, id.vars = c("Row.names", "signature"), variable.name = "cell", 
        value.name = "vst")
    scores = vst %>% group_by(cell, signature) %>% summarize(score = mean(vst))
    
    final_scores = recast(scores, signature ~ cell)
    final_scores = t(final_scores)
    
    colnames(final_scores) = final_scores[1, ]
    final_scores = as.data.frame(final_scores[-1, ])
    
    for (column in colnames(final_scores)) {
        final_scores[, column] <- as.numeric(as.character(final_scores[, 
            column]))
    }
    final_scores = final_scores
    return(final_scores)
}

get_expression_matrix <- function(cds, gene_list) 
{
    vst = as.data.frame(as.matrix(log2((Matrix::t(Matrix::t(Biobase:::exprs(cds))/pData(cds)$Size_Factor)) + 1)))

    vst = vst[gene_list$id, ]
    
    vst = merge(vst, gene_list[, names(gene_list) %in% c("id", 
        "gene_short_name")], by.x = "row.names", by.y = "id")
    
    vst = melt(vst, id.vars = c("Row.names", "gene_short_name"), variable.name = "cell", 
        value.name = "vst")
#     scores = vst %>% group_by(cell, signature) %>% summarize(score = mean(vst))
    
    return(vst)

}


# Load Mock and TGFB cds objects created in Figure1 code
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

# Identify genes expressed in at least 50 cells
expressed_genes.list <- list()

expressed_genes.list[["Mock"]] <- row.names(fData(cds.list[["Mock"]])[Matrix::rowSums(Biobase::exprs(cds.list[["Mock"]]) > 0) > 50 ,])
length(expressed_genes.list[["Mock"]])
expressed_genes.list[["TGFB"]] <- row.names(fData(cds.list[["TGFB"]])[Matrix::rowSums(Biobase::exprs(cds.list[["TGFB"]]) > 0) > 50 ,])
length(expressed_genes.list[["TGFB"]])

# Identify MCF10A cells expressed at the beginning and end of spontaneous (mock) and TGFB-driven (tgfb) pseudospatial trajectories
Mock_inner_subset_cells <- row.names(pData(cds.list[["Mock"]])[pData(cds.list[["Mock"]])$pseudospace_quantile == 1,])
length(Mock_inner_subset_cells)
Mock_outer_subset_cells <- row.names(pData(cds.list[["Mock"]])[pData(cds.list[["Mock"]])$pseudospace_quantile == 10,])
length(Mock_outer_subset_cells)
TGFB_inner_subset_cells <- row.names(pData(cds.list[["TGFB"]])[pData(cds.list[["TGFB"]])$pseudospace_quantile == 1,])
length(TGFB_inner_subset_cells)
TGFB_outer_subset_cells <- row.names(pData(cds.list[["TGFB"]])[pData(cds.list[["TGFB"]])$pseudospace_quantile == 10,])
length(TGFB_outer_subset_cells)

#### Load clean datasets ####
# Unecessary for MCF10A, change this later
pseudospace_cds <- readRDS("pseudospace_cds.rds")

pseudospace_cds <- pseudospace_cds[,c(Mock_inner_subset_cells,Mock_outer_subset_cells,
                                      TGFB_inner_subset_cells,TGFB_outer_subset_cells)]

# Load HNSCC cds object
HNSCC_cds <- readRDS("HSNCC_cds.rds")

# Filter HNSCC data for samples with more than 40, non-lymph node, cancer cells processesd with Maxima enzyme
HNSCC_patient_list <- as.data.frame(pData(HNSCC_cds) %>% 
filter(Maxima_enzyme == "0" & clasified_as_cancer_cell == "1" & lymph_node == "0") %>%
group_by(patient_id) %>% summarize(n = n()) %>% arrange(desc(n)))

HNSCC_patients_forAnalysis <- HNSCC_patient_list[HNSCC_patient_list$n > 40,]$patient_id

HNSCC_cds <- HNSCC_cds[,pData(HNSCC_cds)$patient_id %in% HNSCC_patients_forAnalysis & 
                       pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                       pData(HNSCC_cds)$clasified_as_cancer_cell == "1" &
                      pData(HNSCC_cds)$lymph_node == "0"]

# Create a joint MCF10A-HNSCC cds object

HNSCC_fData <- fData(HNSCC_cds)
MCF10A_fData <- fData(pseudospace_cds)[,c("id","gene_short_name")]

pData(HNSCC_cds)$cell <- as.character(row.names(pData(HNSCC_cds)))
HNSCC_pData <- pData(HNSCC_cds)[,c("cell","patient_id")]
MCF10A_pData <- pData(pseudospace_cds)[,c("cell","sample")]
colnames(HNSCC_pData) <- c("cell","sample")

HNSCC_exprs <- as(round(Biobase::exprs(HNSCC_cds)),"dgTMatrix")
MCF10A_exprs <- Biobase::exprs(pseudospace_cds)
MCF10A_exprs <- MCF10A_exprs[row.names(HNSCC_exprs),]

MCF10A_fData <- MCF10A_fData[row.names(HNSCC_fData),]

# Create joint exprs matrices as well as metadata and feature data data.frames

new_pData <- rbind(MCF10A_pData, HNSCC_pData)
new_pData <- new("AnnotatedDataFrame", data = new_pData)
new_fData <- new("AnnotatedDataFrame", data = MCF10A_fData)
new_exprs <- cbind(MCF10A_exprs,HNSCC_exprs)

HNSCC_MCF10A_cds <- newCellDataSet(as(as.matrix(new_exprs), "sparseMatrix"),
                       phenoData = new_pData,
                       featureData = new_fData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

pData(HNSCC_MCF10A_cds)$cell <- as.character(pData(HNSCC_MCF10A_cds)$cell)

# Add a spatial sample id for MCF10A cells
sample_id <- sapply(pData(HNSCC_MCF10A_cds)$cell, function(x){
    
    if(x %in% Mock_inner_subset_cells)return("Mock_early")
    if(x %in% Mock_outer_subset_cells)return("Mock_late")
    if(x %in% TGFB_inner_subset_cells)return("TGFB_early")
    if(x %in% TGFB_outer_subset_cells)return("TGFB_late")
    return(pData(HNSCC_MCF10A_cds)[x,]$sample)
    
})


pData(HNSCC_MCF10A_cds)$Cell.Type <- sample_id

# Pre-process the new cds object
HNSCC_MCF10A_cds <- estimateSizeFactors(HNSCC_MCF10A_cds)
HNSCC_MCF10A_cds <- estimateDispersions(HNSCC_MCF10A_cds)

# Load results of differential gene expression test and clustering pertaining to Figure 2c. 
pseudo.DTW.gene.clusters.aucRank <- readRDS("pseudo.DTW.gene.clusters.aucRank.rds")

# Isolate genes that define an early EMT induction and 
# KRAS ASSOCIATED induction of a mature mesenchymal state 
Early_EMT_genes <- names(pseudo.DTW.gene.clusters.aucRank[pseudo.DTW.gene.clusters.aucRank == 1])
Late_EMT_genes <- names(pseudo.DTW.gene.clusters.aucRank[pseudo.DTW.gene.clusters.aucRank == 5])


# Load MSigDB hallmarks and GO biological process geneset collections
source("loadGSCSafe.R")
hallmarksGSC <- loadGSCSafe(file="h.all.v6.0.symbols.gmt")
GOGSC<-loadGSCSafe(file="Human_GO_bp_no_GO_iea_symbol.gmt")

# Isolate genes from genesets associated with induction of an EMT and epidermis development
EMT_marker_gene_names <- hallmarksGSC$gsc$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
EMT_marker_id <- unique(fData(cds.list[["Mock"]])[fData(cds.list[["Mock"]])$gene_short_name %in% EMT_marker_gene_names,]$id)

epithelial_marker_gene_names <- GOGSC$gsc$'EPIDERMIS DEVELOPMENT%GOBP%GO:0008544'
epithelial_marker_id <- unique(fData(cds.list[["Mock"]])[fData(cds.list[["Mock"]])$gene_short_name %in% epithelial_marker_gene_names,]$id)

Epithelial_marker_genes <- epithelial_marker_id
Canonical_EMT_marker_genes <- EMT_marker_id

# Create a data.frame consisting of all signatures
EMT_signature_df <- as.data.frame(matrix(c(Epithelial_marker_genes,
                                                 Early_EMT_genes,
                                                 Late_EMT_genes,
                                                 Canonical_EMT_marker_genes,
                                                 rep("Epithelial_score",
                                                     length(Epithelial_marker_genes)),
                                                 rep("Early_EMT_score", 
                                                     length(Early_EMT_genes)),
                                                 rep("Late_EMT_score", 
                                                     length(Late_EMT_genes)),
                                                 rep("Canonical_EMT_score", 
                                                     length(Canonical_EMT_marker_genes))), 
                                               ncol = 2))
colnames(EMT_signature_df) <- c("id","signature")
EMT_signature_df$id <- as.character(EMT_signature_df$id)
EMT_signature_df$signature <- as.character(EMT_signature_df$signature)

# Calculate and scale scores for every cell across every signature
HNSCC_MCF10A_EMT_scores <- calculate_signature_scores(HNSCC_MCF10A_cds,EMT_signature_df)
HNSCC_MCF10A_EMT_scores$Sample <- pData(HNSCC_MCF10A_cds)$Cell.Type
HNSCC_MCF10A_EMT_scores_scaled <- HNSCC_MCF10A_EMT_scores %>% mutate_at(c(1,2,3,4), funs(c(scale(.))))

HNSCC_MCF10A_EMT_scores_scaled %>% 
group_by(Sample) %>% 
summarize(median_EMT_score = median(Canonical_EMT_score), median_Epithelial_score = median(Epithelial_score),
         median_Early_EMT_score = median(Early_EMT_score), median_Late_EMT_score = median(Late_EMT_score)) %>%
arrange(median_EMT_score)

HNSCC_MCF10A_EMT_scores_scaled$Sample <- factor(HNSCC_MCF10A_EMT_scores_scaled$Sample, 
                                   levels = c("Mock_early","Mock_late",
                                              "TGFB_early","TGFB_late",
                                              "HNSCC6","HNSCC20","HNSCC5",
                                                    "HNSCC18","HNSCC22","HNSCC25","HNSCC17",
                                                    "HNSCC16"))

HNSCC_MCF10A_EMT_scores_scaled %>% 
filter(Sample %in% c("Mock_early","TGFB_late")) %>% 
group_by(Sample) %>% 
summarize(median_EMT_score = median(Canonical_EMT_score), median_Epithelial_score = median(Epithelial_score),
         median_Early_EMT_score = median(Early_EMT_score), median_Late_EMT_score = median(Late_EMT_score))

ggplot(HNSCC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Canonical_EMT_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = -1.482951, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = 1.143410, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
ylab("Hallmarks EMT\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HNSCC_MCF10A_Canonical_EMT_score.png", height = 4.2, width = 6)

ggplot(HNSCC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Epithelial_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = 1.0393513, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = -0.9142161, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
#ylim(-3,5) +
ylab("GO-BP epithelial\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HNSCC_MCF10A_Epithelial_score.png", height = 4.2, width = 6)

ggplot(HNSCC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Early_EMT_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = -1.313916, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = 0.738604, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
ylab("MCF10A early EMT\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HNSCC_MCF10A_Early_EMT_score.png", height = 4.2, width = 6)

ggplot(HNSCC_MCF10A_EMT_scores_scaled, aes(x = Sample, y = Late_EMT_score)) + 
geom_boxplot(fill = "gray70") +
geom_jitter(size = 0.01, color = "#702632", alpha = 0.1) +
geom_hline(yintercept = -1.589488, "#0075F2", linetype = "dashed") +
geom_hline(yintercept = 1.523222, color = "#70163C", linetype = "dashed") +
theme(text=element_text(size=24),
      axis.text.x = element_text(angle = 45, hjust = 1),
     axis.title.x = element_blank()) +
ylab("MCF10A late EMT\nscore") +
monocle:::monocle_theme_opts() +
ggsave(file = "HNSCC_MCF10A_Late_EMT_score.png", height = 4.2, width = 6)

# Generate heatmaps of the expression of genes whose loss leads to accumulation of MCF10A cells 
# across pseudospace in HNSCC cells and investigate their correlation with partial-EMT rank from Puram et al.
all_enriched_target.list <- list()

all_enriched_target.list[["Mock"]] <- c('MYC','TWIST1','ITGB1','FOXQ1',
  'FOXC2','FOXO3','TWIST2','FOXA1',
  'KLF8','FGFR1','FZD7','IL6R',
  'POU5F1','MET','FZD4','FOXD3',
  'PRRX1','TRIM28','EGFR','GATA6',
  'ITGAV','ZEB2','NOTCH1','PTCH1','SNAI2')

all_enriched_target.list[["TGFB"]] <- c('FGFR1','FOXC2','GATA6','POU5F1',
  'ZEB2','NOTCH1','FZD2','MET',
  'PTCH1','TGFBR1','TGFBR2','ITGAV',
  'ZNF703','ZEB1','PRRX1','TWIST2')

all_enriched_targets <- unique(union(all_enriched_target.list[["Mock"]],all_enriched_target.list[["TGFB"]]))

enriched_targets_receptors <- c("EGFR","FGFR1","FZD2","FZD4",
                                "FZD7","IL6R","ITGAV","ITGB1","MET",
                                "NOTCH1","PTCH1","TGFBR1","TGFBR2")

enriched_targets_transcription_factors <- c("FOXA1","FOXC2","FOXD3","FOXO3",
                                            "FOXQ1","GATA6","KLF8","MYC",
                                            "POU5F1","PRRX1","SNAI2","TRIM28",
                                            "TWIST1","TWIST2","ZEB1","ZEB2",
                                            "ZNF703")

HNSCC_cds <- estimateSizeFactors(HNSCC_cds)

# Filter for genes expressed in HNSCC cells
HNSCC_expressed_genes <- row.names(fData(HNSCC_cds)[Matrix::rowSums(Biobase::exprs(HNSCC_cds[,pData(HNSCC_cds)$patient_id %in% HNSCC_patients_forAnalysis, 
                                                                                             pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                                                                                             pData(HNSCC_cds)$clasified_as_cancer_cell == "1" &
                                                                                             pData(HNSCC_cds)$lymph_node == "0"]) > 0) > 50 ,])
length(HNSCC_expressed_genes)

enriched_targets_receptors <- unique(intersect(as.character(fData(HNSCC_cds)[fData(HNSCC_cds)$gene_short_name %in% enriched_targets_receptors,]$id), HNSCC_expressed_genes))
length(enriched_targets_receptors)

enriched_targets_transcription_factors <- unique(intersect(as.character(fData(HNSCC_cds)[fData(HNSCC_cds)$gene_short_name %in% enriched_targets_transcription_factors,]$id), HNSCC_expressed_genes))
length(enriched_targets_transcription_factors)

enriched_targets_receptors <- fData(HNSCC_cds)[enriched_targets_receptors,]$gene_short_name
enriched_targets_transcription_factors <- fData(HNSCC_cds)[enriched_targets_transcription_factors,]$gene_short_name

receptors_enriched_fData <- fData(HNSCC_cds)[fData(HNSCC_cds)$gene_short_name %in% enriched_targets_receptors,]
transcriptiona_factors_enriched_fData <- fData(HNSCC_cds)[fData(HNSCC_cds)$gene_short_name %in% enriched_targets_transcription_factors,]


receptors_enriched_fData$id <- as.character(receptors_enriched_fData$id)
receptors_enriched_fData$gene_short_name <- as.character(receptors_enriched_fData$gene_short_name)

transcriptiona_factors_enriched_fData$id <- as.character(transcriptiona_factors_enriched_fData$id)
transcriptiona_factors_enriched_fData$gene_short_name <- as.character(transcriptiona_factors_enriched_fData$gene_short_name)

receptors_enriched_targets_HNSCC_counts <- get_expression_matrix(HNSCC_cds, receptors_enriched_fData)
transcription_factors_enriched_targets_HNSCC_counts <- get_expression_matrix(HNSCC_cds, transcriptiona_factors_enriched_fData)

receptors_enriched_targets_HNSCC_counts <- dcast(receptors_enriched_targets_HNSCC_counts, 
                                            cell~gene_short_name, value.var = "vst")
transcription_factors_enriched_targets_HNSCC_counts <- dcast(transcription_factors_enriched_targets_HNSCC_counts, 
                                            cell~gene_short_name, value.var = "vst")

HNSCC_pData_receptors_enriched_counts <- merge(pData(HNSCC_cds), receptors_enriched_targets_HNSCC_counts, by = "cell")
HNSCC_pData_transcription_factors_enriched_counts <- merge(pData(HNSCC_cds), transcription_factors_enriched_targets_HNSCC_counts, by = "cell")


row.names(HNSCC_pData_receptors_enriched_counts) <- HNSCC_pData_receptors_enriched_counts$cell
row.names(HNSCC_pData_transcription_factors_enriched_counts) <- HNSCC_pData_transcription_factors_enriched_counts$cell

head(HNSCC_pData_receptors_enriched_counts[,c("patient_id",as.character(enriched_targets_receptors))])

HNSCC_pData_receptors_enriched_counts <- HNSCC_pData_receptors_enriched_counts[,c("patient_id",as.character(enriched_targets_receptors))]
HNSCC_pData_transcription_factors_enriched_counts <- HNSCC_pData_transcription_factors_enriched_counts[,c("patient_id", as.character(enriched_targets_transcription_factors))]


HNSCC_average_receptor_enriched_counts <- as.data.frame(HNSCC_pData_receptors_enriched_counts %>% 
group_by(patient_id) %>% 
summarize_all(funs(mean)))

HNSCC_average_transcription_factors_enriched_counts <- as.data.frame(HNSCC_pData_transcription_factors_enriched_counts %>% 
group_by(patient_id) %>% 
summarize_all(funs(mean)))

row.names(HNSCC_average_receptor_enriched_counts) <- as.character(HNSCC_average_receptor_enriched_counts$patient_id)
HNSCC_average_receptor_enriched_counts <- HNSCC_average_receptor_enriched_counts[,-1]

row.names(HNSCC_average_transcription_factors_enriched_counts) <- as.character(HNSCC_average_transcription_factors_enriched_counts$patient_id)
HNSCC_average_transcription_factors_enriched_counts <- HNSCC_average_transcription_factors_enriched_counts[,-1]


head(HNSCC_average_receptor_enriched_counts)
head(HNSCC_average_transcription_factors_enriched_counts)

HNSCC_average_receptor_enriched_counts_scaled <- t(scale(t(scale(HNSCC_average_receptor_enriched_counts))))
HNSCC_average_transcription_factors_enriched_counts_scaled <- t(scale(t(scale(HNSCC_average_transcription_factors_enriched_counts))))

HNSCC_EMT_annotation <- data.frame(row.names = c("HNSCC6","HNSCC20","HNSCC5","HNSCC18","HNSCC22","HNSCC25","HNSCC17","HNSCC16"),
                                  patient_id = c("HNSCC6","HNSCC20","HNSCC5","HNSCC18","HNSCC22","HNSCC25","HNSCC17","HNSCC16"),
                                  EMT_rank = c("1","2","3","4","5","6","7","8"))

pheatmap(t(HNSCC_average_receptor_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = FALSE,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_receptors.png",
        width = 3.5, height = 3)

receptors_ph <- pheatmap(t(HNSCC_average_receptor_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = TRUE,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_receptors_clustered.png",
        width = 3.5, height = 4)

# Flip clades around in a manner that is proper and highlights the relative EMT rank
receptors_ph$tree_col$order

new_receptor_tree_col_order <- receptors_ph$tree_col
new_receptor_tree_col_order$order <- c(8,4,3,5,6,7,2,1)
plot(receptors_ph$tree_col)
plot(new_receptor_tree_col_order)

pheatmap(t(HNSCC_average_receptor_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = new_receptor_tree_col_order,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_receptors_clustered.png",
        width = 4.2, height = 4)

pheatmap(t(HNSCC_average_transcription_factors_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = FALSE,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_transcription_factors.png",
        width = 3.5, height = 3)

TFs_ph <- pheatmap(t(HNSCC_average_transcription_factors_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = TRUE,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_transcription_factors_clustered.png",
        width = 3.5, height = 4)

TFs_ph$tree_col$order

new_TF_tree_col_order <- TFs_ph$tree_col
new_TF_tree_col_order$order <- c(8,4,5,6,2,3,7,1)
plot(TFs_ph$tree_col)
plot(new_TF_tree_col_order)

pheatmap(t(HNSCC_average_transcription_factors_enriched_counts_scaled),
        show_colnames = TRUE,
        show_rownames = TRUE,
        cluster_cols = new_TF_tree_col_order,
        annotation_col = data.frame(row.names = row.names(HNSCC_EMT_annotation),
                                    EMT_rank = as.numeric(HNSCC_EMT_annotation$EMT_rank)),
        color = inferno(25),
        clustering_method = "ward.D2",
        file = "Mean_expression_levels_of_enriched_transcription_factors_clustered.png",
        width = 4.2, height = 4)




