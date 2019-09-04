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

getDTWcds <- function (query_cds, ref_cds, ref, query, expressed_genes, cores) 
{
    alignment_genes <- intersect(row.names(subset(fData(ref_cds), 
        use_for_ordering)), row.names(subset(fData(query_cds), 
        use_for_ordering)))
    ref_align_cds <- ref_cds[alignment_genes]
    query_align_cds <- query_cds[alignment_genes]
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

    ref_query_dtw <- align_cells(ref_matrix, query_matrix, step_pattern = rabinerJuangStepPattern(3, 
        "c"), open.begin = F, open.end = F)
    align_res <- warp_pseudotime(ref_align_cds, query_align_cds, 
        ref_query_dtw)
    query_ref_aligned <- align_res$query_cds
    pData(query_ref_aligned)$Pseudotime <- pData(query_ref_aligned)$Alignment_Pseudotime
    ref_aligned_cell_ids <- setdiff(row.names(pData(ref_align_cds)), 
        "duplicate_root")
    query_aligned_cell_ids <- setdiff(row.names(pData(query_align_cds)), 
        "duplicate_root")
    combined_exprs <- cBind(exprs(query_cds[expressed_genes, 
        query_aligned_cell_ids]), exprs(ref_cds[expressed_genes, 
        ref_aligned_cell_ids]))
    pData_ref <- pData(ref_align_cds)[, c("treatment_id", "spatial_id", 
        "Pseudotime")]
    pData_ref$Cell.Type <- ref
    pData_query_aligned <- pData(query_ref_aligned)[, c("treatment_id", 
        "spatial_id", "Pseudotime")]
    pData_query_aligned$Cell.Type <- query
    combined_pData <- rbind(pData_query_aligned, pData_ref)
    combined_pData <- combined_pData[colnames(combined_exprs), 
        ]
    combined_pd <- new("AnnotatedDataFrame", data = combined_pData)
    
    fd <- new("AnnotatedDataFrame", data = fData(ref_cds)[row.names(combined_exprs), 
        1:2])

    ref_queryToRef_combined_cds <- newCellDataSet(combined_exprs, 
                                                  phenoData = combined_pd, 
                                                  featureData = fd, 
                                                  expressionFamily = negbinomial.size(), 
        lowerDetectionLimit = 1)
    pData(ref_queryToRef_combined_cds)$cell_id <- row.names(pData(ref_queryToRef_combined_cds))
    return(ref_queryToRef_combined_cds)
}

compare_cell_types_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
    min_expr = NULL, cell_size = 0.75, cell_alpha = 1, line_size = 1, 
    nrow = NULL, ncol = 1, panel_order = NULL, color_by = "Cell.Type", 
    shade_by = NULL, df = 3, maxit = 300, relative_expr = TRUE, 
    pseudocount = 0) 
{
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- cds_exprs + pseudocount
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- t(t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- melt(ceiling(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- cds_exprs + pseudocount
        cds_exprs <- melt(as.matrix(cds_exprs))
    }
    colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
        by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    if (integer_expression) {
        cds_exprs$expression <- cds_exprs$expression
    }
    else {
        cds_exprs$expression <- log10(cds_exprs$expression)
    }
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
        cds_exprs$gene_label <- as.character(cds_exprs$gene_short_name)
        cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
    }
    else {
        cds_exprs$gene_label <- cds_exprs$gene_id
    }
    cds_exprs$gene_label <- factor(cds_exprs$gene_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, 
        Cell.Type = pData(cds_subset)$Cell.Type)
    model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
        trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
        x[1]])
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    if (is.null(panel_order) == FALSE) {
        cds_exprs$gene_label <- factor(cds_exprs$gene_label, 
            levels = panel_order)
    }
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
    if (is.null(color_by) == FALSE) {
        q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), 
            alpha = I(cell_alpha))
    }
    else {
        q <- q + geom_point(size = I(cell_size))
    }
    q <- q + geom_line(aes(Pseudotime, expectation, color = Cell.Type), 
        size = line_size)
    q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow = nrow, 
        ncol = ncol, scales = "free_y")
    q <- q + ylab("Expression") + xlab("Pseudospace")
    q <- q + monocle:::monocle_theme_opts()
    q
}

delta_auc_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
    min_expr = NULL, df = 3, maxit = 300, relative_expr = TRUE, 
    pseudocount = 0, cores = 1) 
{
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
        "negbinomial.size")) {
        integer_expression <- TRUE
    }
    else {
        integer_expression <- FALSE
        relative_expr <- TRUE
    }
    if (integer_expression) {
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- cds_exprs + pseudocount
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- t(t(cds_exprs)/sizeFactors(cds_subset))
        }
        cds_exprs <- melt(ceiling(as.matrix(cds_exprs)))
    }
    else {
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- cds_exprs + pseudocount
        cds_exprs <- melt(as.matrix(cds_exprs))
    }
    colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
        by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    if (integer_expression) {
        cds_exprs$expression <- cds_exprs$expression
    }
    else {
        cds_exprs$expression <- log10(cds_exprs$expression)
    }
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
        cds_exprs$gene_label <- as.character(cds_exprs$gene_short_name)
        cds_exprs$gene_label[is.na(cds_exprs$gene_label)] <- cds_exprs$gene_id
    }
    else {
        cds_exprs$gene_label <- cds_exprs$gene_id
    }
    cds_exprs$gene_label <- factor(cds_exprs$gene_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime, 
        Cell.Type = pData(cds_subset)$Cell.Type)
    model_expectation <- genSmoothCurves(cds_subset, cores = cores, 
        trend_formula = trend_formula, relative_expr = T, new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
        x[1]])
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    pseudotime_max = max(cds_exprs$Pseudotime)
    cds_exprs.auc = cds_exprs %>% dplyr::select(id, Cell.Type, 
        Pseudotime, expectation) %>% dplyr::mutate(pseudotime_quantile = cut(Pseudotime, 
        breaks = seq(0, pseudotime_max, pseudotime_max/5))) %>% 
        unique() %>% dplyr::group_by(id, Cell.Type, pseudotime_quantile) %>% 
        dplyr::summarize(AUC = MESS::auc(Pseudotime, expectation, 
            type = "spline")) %>% dplyr::group_by(id, pseudotime_quantile) %>% 
        dplyr::summarize(auc_difference = (AUC[1] - AUC[2])/(AUC[1] + 
            AUC[2])) %>% dplyr::arrange(desc(auc_difference))
    cds_exprs.auc$pseudotime_quantile = as.numeric(factor(cds_exprs.auc$pseudotime_quantile, 
        levels = sort(unique(cds_exprs.auc$pseudotime_quantile))))
    cds_exprs.auc_rank = tidyr::spread(cds_exprs.auc, key = pseudotime_quantile, 
        value = auc_difference)
    return(as.data.frame(cds_exprs.auc_rank))
}

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


# Load Mock and TGFB cds objects created in Figure1 code
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

# Identify genes that are expressed in at least 50 of cells 
expressed_genes.list <- list()

expressed_genes.list[["Mock"]] <- row.names(fData(cds.list[["Mock"]])[Matrix::rowSums(Biobase::exprs(cds.list[["Mock"]]) > 0) > 50 ,])
length(expressed_genes.list[["Mock"]])
expressed_genes.list[["TGFB"]] <- row.names(fData(cds.list[["TGFB"]])[Matrix::rowSums(Biobase::exprs(cds.list[["TGFB"]]) > 0) > 50 ,])
length(expressed_genes.list[["TGFB"]])

expressed_genes <- unique(union(expressed_genes.list[["Mock"]], expressed_genes.list[["TGFB"]]))
length(expressed_genes)

# Use dynamic time warping to align Mock and TGFB pseudospatial trajectories and create a cds object of aligned trajectories
TGFB.to.Mock.aligned.cds <- getDTWcds(cds.list[["TGFB"]],cds.list[["Mock"]], 
                                      ref = "Mock", query = "TGFB", 
                                      expressed_genes = expressed_genes, cores = 1)
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
    theme(legend.position="none", text=element_text(size=20)) +
    ggsave("EMT_markers_Mock_TGFB_aligned_pseudospace.png", width = 6, height = 5)

# Identify genes that are differentially expressed across pseudospace as a function of treatment
aligned.pseudospace.DEG.test.list <- list()

for(alignment in names(cds.aligned.list)) {

aligned.pseudospace.DEG.test.list[[alignment]] <- differentialGeneTest(cds.aligned.list[[alignment]][expressed_genes], 
                 fullModelFormulaStr="~sm.ns(Pseudotime, df=3)*Cell.Type", 
                 reducedModelFormulaStr="~sm.ns(Pseudotime, df=3) + Cell.Type", cores = 1)
    
    }


Pseudospatial.aligned.sig.genes.list <- list()

for(sample in names(aligned.pseudospace.DEG.test.list)){
    
    Pseudospatial.aligned.sig.genes.list[[sample]] <- row.names(subset(aligned.pseudospace.DEG.test.list[[sample]], 
                                                                       qval <= 1e-10))
    print(sample)
    print(length(Pseudospatial.aligned.sig.genes.list[[sample]]))
}

#### Obtain all auc differences ####

full_tgfb_to_mock_aucRank <- delta_auc_in_pseudospace(cds.aligned.list[["TGFB to Mock"]][row.names(aligned.pseudospace.DEG.test.list[["TGFB to Mock"]])],
                                                                 pseudocount = 1,
                                                                 cores = 20)

full_tgfb_to_mock_aucRank$'<NA>' <- NULL

row.names(full_tgfb_to_mock_aucRank) <- full_tgfb_to_mock_aucRank$id

full_tgfb_to_mock_aucRank$id <- NULL

max_auc_idx <- as.character(apply(abs(full_tgfb_to_mock_aucRank), 1, which.max))
names(max_auc_idx) <- row.names(full_tgfb_to_mock_aucRank)

max_auc_values <- as.data.frame(sapply(row.names(full_tgfb_to_mock_aucRank), function(x){
    
    return(full_tgfb_to_mock_aucRank[x,max_auc_idx[x]])
    
}))

max_auc_values$id <- row.names(max_auc_values)
colnames(max_auc_values) <- c("max_auc_difference","id")

# Merge max bin auc with deg test results
merged_max_delta_auc_DEG_test <- merge(aligned.pseudospace.DEG.test.list[["TGFB to Mock"]],max_auc_values,
                                                 by.x = "id", by.y = "id")

ggplot(merged_max_delta_auc_DEG_test, aes(x = max_auc_difference, y = -log10(qval), color = (qval < 1e-10 & abs(max_auc_difference) > 0.02))) + 
geom_point(alpha = 0.1) +
scale_color_manual("DEG Filter", 
                   values = c("TRUE" = "brown4","FALSE" = "dimgrey"),
                  labels = c("TRUE" = "Pass","FALSE" = "FAIL"))+
xlab("Delta AUC") +
monocle:::monocle_theme_opts() +
theme(text=element_text(size=20)) +
xlim(-1,1) +
ggsave("Relative_AUC_difference_qval_volcano_plot.png", 
    width = 5, height = 4)

# Check the relationship between genes that fail the au cutoff but pass significance and their expression level
sig_fail_genes <- merged_max_delta_auc_DEG_test[merged_max_delta_auc_DEG_test$qval < 1e-10 & 
                              abs(merged_max_delta_auc_DEG_test$max_auc_difference) <= 0.02,]$id

mock_cells <- row.names(pData(cds.aligned.list[["TGFB to Mock"]])[pData(cds.aligned.list[["TGFB to Mock"]])$Cell.Type == "Mock",])
tgfb_cells <- row.names(pData(cds.aligned.list[["TGFB to Mock"]])[pData(cds.aligned.list[["TGFB to Mock"]])$Cell.Type == "TGFB",])


mean(Matrix::rowMeans(Biobase::exprs(cds.aligned.list[["TGFB to Mock"]])[row.names(Biobase::exprs(cds.aligned.list[["TGFB to Mock"]])) %in% 
                                                                         sig_fail_genes,mock_cells]))
mean(Matrix::rowMeans(Biobase::exprs(cds.aligned.list[["TGFB to Mock"]])[row.names(Biobase::exprs(cds.aligned.list[["TGFB to Mock"]])) %in% 
                                                                         sig_fail_genes,tgfb_cells]))

# Identify a reasonable cutoff to filter DEGs by area under the curve
ggplot(max_auc_values, aes(x = abs(max_auc_difference))) + 
geom_histogram() + 
geom_vline(xintercept = 0.02, color = "red") +
xlab("Maximum difference in AUC across quantiles") +
monocle:::monocle_theme_opts() + 
ggsave("AUC_rank_TGFBtoMock_dynamic_time_warping.png")

pseudospatial_aligned_sig_maxRank_genes <- merged_max_delta_auc_DEG_test[merged_max_delta_auc_DEG_test$qval < 1e-10 &
                                                                         abs(merged_max_delta_auc_DEG_test$max_auc_difference) > 0.02,]$id


length(pseudospatial_aligned_sig_maxRank_genes)

aligned.pseudospace.smoothed.combined.exprs.list <- list()

for(alignment in names(cds.aligned.list)){
    
aligned.pseudospace.smoothed.combined.exprs.list[[alignment]] <- genSmoothCurves(cds.aligned.list[[alignment]][pseudospatial_aligned_sig_maxRank_genes], 
                                           new_data=data.frame(Pseudotime=rep(seq(0,99, by=1), 2), 
                                                      Cell.Type=c(rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[1], 100), 
                                                                  rep(unique(pData(cds.aligned.list[[alignment]])$Cell.Type)[2], 100))), 
                                           trend_formula="~sm.ns(Pseudotime, df=3)*Cell.Type",
                                           cores=1)
    
    }

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.list)){
    
aligned.pseudospace.smoothed.combined.exprs.list[[alignment]] <- aligned.pseudospace.smoothed.combined.exprs.list[[alignment]][rowSums(is.na(aligned.pseudospace.smoothed.combined.exprs.list[[alignment]])) == 0,]
    
}

aligned.pseudospace.smoothed.combined.exprs.scaled.list <- list()

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.list)){
    
    aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] <- scale(t(scale(t(aligned.pseudospace.smoothed.combined.exprs.list[[alignment]]))))
    aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] > 3] <- 3
    aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]] < -3] <- -3
    aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]][is.na(aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]])] <- 0
    
    }

row_dist.list <- list()

for(alignment in names(aligned.pseudospace.smoothed.combined.exprs.scaled.list)){

    row_dist.list[[alignment]] <- as.dist((1 - cor(Matrix::t(aligned.pseudospace.smoothed.combined.exprs.scaled.list[[alignment]])))/2)
    row_dist.list[[alignment]][is.na(row_dist.list[[alignment]])] <- 1
    
    }

dim(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['TGFB to Mock']])

num_clusters <- 10
bk <- seq(-3.1,3.1, by=0.1)
hmcols <-  blue2green2red(length(bk) - 1)
ph_aucRank <- pheatmap(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['TGFB to Mock']], 
         useRaster = T,
         cluster_cols=FALSE, 
         cluster_rows=TRUE, 
         show_rownames=F, 
         show_colnames=F, 
         clustering_distance_rows=row_dist.list[["TGFB to Mock"]],
         clustering_method = "ward.D2",
         cutree_rows=num_clusters,
         color=hmcols,
         filename="DTW_DEGs_Mock_TGFB_heatmap_aucRank.png",
         )


aucRank_annotation_row <- data.frame(Cluster=factor(cutree(ph_aucRank$tree_row, 10)))


ann_colors <- list(Cluster = c("1"="brown3", "2"="darkgreen", "3"="darkblue", "4"="darkorchid4", "5"="firebrick1", "6"="goldenrod1", "7"="darkcyan", "8"="chocolate","9" = "deeppink1", "10" = "darkolivegreen1", "11" = "plum4", "12" = "slategray3"))

num_clusters <- 10
bk <- seq(-3.1,3.1, by=0.1)
hmcols <-  inferno(length(bk) - 1)
aucRank_TGFB_ph_res <- pheatmap(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['TGFB to Mock']][,1:100], 
         useRaster = TRUE,
         cluster_cols=FALSE, 
         cluster_rows=ph_aucRank$tree_row, 
         show_rownames=F, 
         show_colnames=F, 
         cutree_rows=num_clusters,
         annotation_row=aucRank_annotation_row,
         annotation_colors=ann_colors,
         treeheight_row = 1.5, 
         breaks=bk,
         fontsize = 6,
         color=hmcols, 
         silent=TRUE,
         filename=NA
         )

png("DTW_DEGs_MocktoTGFB_TGFB_heatmap_aucRank.png", width=2.5, height=3, units = "in",res = 2000)
aucRank_TGFB_ph_res  
grid::grid.rect(gp=grid::gpar("fill", col=NA))
grid::grid.draw(aucRank_TGFB_ph_res$gtable)
dev.off()

num_clusters <- 10
bk <- seq(-3.1,3.1, by=0.1)
hmcols <-  inferno(length(bk) - 1)

aucRank_mock_ph_res <- pheatmap(aligned.pseudospace.smoothed.combined.exprs.scaled.list[['TGFB to Mock']][,101:200],
         useRaster = TRUE,
         cluster_cols=FALSE, 
         cluster_rows=ph_aucRank$tree_row, 
         show_rownames=F, 
         show_colnames=F, 
         cutree_rows=num_clusters,
         annotation_row=aucRank_annotation_row,
         annotation_colors=ann_colors,
         treeheight_row = 1.5, 
         breaks=bk,
         fontsize = 6,
         color=hmcols, 
         silent=TRUE,
         filename=NA
         )



png("DTW_DEGs_MocktoTGFB_Mock_heatmap_aucRank.png", width=2.5, height=3, units = "in",res = 2000)
aucRank_mock_ph_res
grid::grid.rect(gp=grid::gpar("fill", col=NA))
grid::grid.draw(aucRank_mock_ph_res$gtable)
dev.off()

# Load necessary functions
source("loadGSCSafe.R")
source('plot_gsea_go.R')
source("GSA_hyper_helper_functions.R")

## Load Gene Set Collections
OncogeneSignatureGSC <- loadGSCSafe(file="c6.all.v6.0.OncogenicSignatures.symbols.gmt")
hallmarksGSC <- loadGSCSafe(file="Pathways/h.all.v6.0.symbols.gmt")
GOGSC<-loadGSCSafe(file="Human_GO_bp_no_GO_iea_symbol.gmt")

Ensembl_GSAlist <- as.matrix(fData(cds.aligned.list[["TGFB to Mock"]][expressed_genes])$gene_short_name)
rownames(Ensembl_GSAlist)<-row.names(fData(cds.aligned.list[["TGFB to Mock"]][expressed_genes]))
colnames(Ensembl_GSAlist) <- c("gene_short_name")
Ensembl_GSAlist<-Ensembl_GSAlist[,1]
Ensembl_GSAlist<-toupper(Ensembl_GSAlist)
length(Ensembl_GSAlist)

pseudo.DTW.gene.clusters.aucRank <- cutree(ph_aucRank$tree_row, 10)
saveRDS(pseudo.DTW.gene.clusters.aucRank,"pseudo.DTW.gene.clusters.aucRank.rds")

TGFBtoMock_aucRank_oncogeneSignature_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                    replace_gene_names_vec(pseudo.DTW.gene.clusters.aucRank
                                                                           ,Ensembl_GSAlist),
                                                    OncogeneSignatureGSC)

TGFBtoMock_aucRank_hallmarks_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                    replace_gene_names_vec(pseudo.DTW.gene.clusters.aucRank
                                                                           ,Ensembl_GSAlist),
                                                    hallmarksGSC)

TGFBtoMock_aucRank_GOGSC_GSAhyper <- collect_gsa_hyper_results_clusters(Ensembl_GSAlist,
                                                    replace_gene_names_vec(pseudo.DTW.gene.clusters.aucRank
                                                                           ,Ensembl_GSAlist),
                                                    GOGSC)

gsea_bar_plots(TGFBtoMock_aucRank_oncogeneSignature_GSAhyper, qval_cutoff = 0.05, width = 6, height = 10, pattern = NULL,sample = "TGFB to Mock", gsc = "OncogeneSignature")

gsea_bar_plots(TGFBtoMock_aucRank_hallmarks_GSAhyper, qval_cutoff = 0.05, width = 6, height = 10, pattern = "HALLMARK_",sample = "TGFB to Mock", gsc = "Hallmarks")


gsea_bar_plots(TGFBtoMock_aucRank_GOGSC_GSAhyper, qval_cutoff = 0.0001, pattern = "%GOBP%GO:", width = 6, height = 20, sample = "TGFB to Mock", gsc = "GOBP")


OncogeneSignatures_sig_results <- get_gsea_sig_results(TGFBtoMock_aucRank_oncogeneSignature_GSAhyper, 0.05)
Hallmarks_sig_results <- get_gsea_sig_results(TGFBtoMock_aucRank_hallmarks_GSAhyper, 0.05)
GOBP_sig_results <- get_gsea_sig_results(TGFBtoMock_aucRank_GOGSC_GSAhyper, 0.05)


### Rename clusters according to Figure 2c labels
# Order in geneset enrichment results
names(OncogeneSignatures_sig_results)
# order to map back to heatmap in Figure 2c
c("6","10","3","9","8","7","1","5","4","2")

names(OncogeneSignatures_sig_results) <- c("6","10","3","9","8","7","1","5","4","2")
names(Hallmarks_sig_results) <- c("6","10","3","9","8","7","1","5","4","2")
names(GOBP_sig_results) <- c("6","10","3","9","8","7","1","5","4","2")

OncogeneSignatures_sig_results <- do.call("rbind",OncogeneSignatures_sig_results)
OncogeneSignatures_sig_results$cluster <- sapply(row.names(OncogeneSignatures_sig_results), function(x){substring(x, 1, 1)})

Hallmarks_sig_results <- do.call("rbind",Hallmarks_sig_results)
Hallmarks_sig_results$cluster <- sapply(row.names(Hallmarks_sig_results), function(x){substring(x, 1, 1)})

GOBP_sig_results <- do.call("rbind",GOBP_sig_results)
GOBP_sig_results$cluster <- sapply(row.names(GOBP_sig_results), function(x){substring(x, 1, 1)})


# Plot interesting genes from gene sets that are significantly enriched across clusters
genes_of_interest <- fData(cds.aligned.list[["TGFB to Mock"]])[fData(cds.aligned.list[["TGFB to Mock"]])$gene_short_name %in% 
                                                               c("NNMT","SCG5","TAGLN","DST",
                                                                 "CALD1","ITGA2"),]$id
fData(cds.aligned.list[["TGFB to Mock"]])[names(pseudo.DTW.gene.clusters.aucRank[names(pseudo.DTW.gene.clusters.aucRank) %in% 
                                                                           genes_of_interest]),]

compare_cell_types_in_pseudospace(cds.aligned.list[["TGFB to Mock"]][row.names(subset(fData(cds.aligned.list[["TGFB to Mock"]]), 
                                                                                    gene_short_name %in% 
                                                                                      c("NNMT","DST",
                                                                                        "TAGLN","SCG5",
                                                                                        "CALD1","ITGA2"))),],
                                  panel_order = c("NNMT","DST",
                                                  "CALD1","ITGA2",
                                                  "TAGLN","SCG5"),
                                  color_by="Cell.Type", df=3, min_expr=0.1, cell_alpha = 0.01, line_size = 1, ncol = 2) +
scale_color_manual(values = c("Mock" = "#0075F2", "TGFB" = "#70163C"), name = "Treatment") +
    theme(legend.position="none", text=element_text(size=18, family="Arial")) +
    ggsave("Final_MesenchymalStateHigh_Mock_TGFB_aligned_pseudospace.png", width = 6, height = 5)
    
cluster_5_gene_id <- row.names(subset(aucRank_annotation_row, Cluster == 5))

cluster_5_genes <- fData(cds.aligned.list[["TGFB to Mock"]])[cluster_5_gene_id,]$gene_short_name

KRAS_upregulated_cluster_5_genes <- intersect(OncogeneSignatureGSC$'gsc'$KRAS.600.LUNG.BREAST_UP.V1_UP, cluster_5_genes)

compare_cell_types_in_pseudospace(cds.aligned.list[["TGFB to Mock"]][row.names(subset(fData(cds.aligned.list[["TGFB to Mock"]]), 
                                                                                    gene_short_name %in% KRAS_upregulated_cluster_5_genes)),],
                                  color_by="Cell.Type", df=3, min_expr=0.1, cell_alpha = 0.01, line_size = 1, ncol = 4) +
scale_color_manual(values = c("Mock" = "#0075F2", "TGFB" = "#70163C"), name = "Treatment") +
    theme(legend.position="none", text=element_text(size=18, family="Arial")) +
    ggsave("KRAS_upregulated_Mock_TGFB_aligned_pseudospace.png", width = 13, height = 4.5)
    

cluster_2_gene_id <- row.names(subset(aucRank_annotation_row, Cluster == 2))

cluster_2_genes <- fData(cds.aligned.list[["TGFB to Mock"]])[cluster_2_gene_id,]$gene_short_name

KRAS_upregulated_cluster_2_genes <- intersect(OncogeneSignatureGSC$'gsc'$KRAS.LUNG_UP.V1_DN, cluster_2_genes)

compare_cell_types_in_pseudospace(cds.aligned.list[["TGFB to Mock"]][row.names(subset(fData(cds.aligned.list[["TGFB to Mock"]]), 
                                                                                    gene_short_name %in% KRAS_upregulated_cluster_2_genes)),],
                                  color_by="Cell.Type", df=3, min_expr=0.1, cell_alpha = 0.01, line_size = 1, ncol = 4) +
scale_color_manual(values = c("Mock" = "#0075F2", "TGFB" = "#70163C"), name = "Treatment") +
    theme(legend.position="none", text=element_text(size=18, family="Arial")) +
    ggsave("KRAS_downregulated_Mock_TGFB_aligned_pseudospace.png", width = 13, height = 4.5)
    
