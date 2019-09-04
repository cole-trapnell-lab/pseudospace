preprocess_cds <- function(cds){
    
    cds <- detectGenes(cds, min_expr = 1)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    return(cds)
}

process_cds_Trajectory <- function(cds, sig_genes, ncenter){
    
        cds <- setOrderingFilter(cds, sig_genes)
    
       if(missing(ncenter)){
            
            cds <- reduceDimension(cds, norm_method = "log")
        }
    
            cds <- reduceDimension(cds, norm_method = "log", ncenter = ncenter)
    
    
    return(cds)
}

diff_foldChange <- function(X, id_to_average_over, normalize_to){
  Group_subset <- list()
  Grouping_list <- unique(pData(X)[,id_to_average_over])
  Grouping_list <- setdiff(Grouping_list, normalize_to)
  
    print(Grouping_list) 
    
    Norm_cds <- X[,pData(X)[,id_to_average_over] == normalize_to]
    Norm_Size_Factor <- Matrix::t(Matrix::t(Biobase::exprs(Norm_cds))/pData(Norm_cds)$Size_Factor)

  for (Group in Grouping_list) {
        
        temp_cds <- X[,pData(X)[,id_to_average_over] == Group]
        temp_Size_Factor <- Matrix::t(Matrix::t(Biobase::exprs(temp_cds))/pData(temp_cds)$Size_Factor)
        temp_foldChange <- as.matrix(log2(Matrix:::rowMeans(temp_Size_Factor)/Matrix:::rowMeans(Norm_Size_Factor)))
    
        Group_subset[[which(Grouping_list == Group)]] <- temp_foldChange
        names(Group_subset) <- paste0("log2FC_",Grouping_list[seq_along(Group_subset)],sep ="")
      
        print(Group)
  }
 return(Group_subset)
}

plot_genes_violin <- function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
    nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, fill_by = NULL,
    plot_trend = FALSE, label_by_short_name = TRUE, relative_expr = TRUE, log_scale = FALSE, 
                               horizontal_orientation = FALSE, trend_size = 0.35) 
{
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
        if (relative_expr) {
            if (is.null(sizeFactors(cds_subset))) {
                stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
            }
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/pData(cds_subset)$Size_Factor)
        }
        #cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
        cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    }
    else {
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
    }
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    if (label_by_short_name == TRUE) {
        if (is.null(cds_exprs$gene_short_name) == FALSE) {
            cds_exprs$feature_label <- cds_exprs$gene_short_name
            cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
        }
        else {
            cds_exprs$feature_label <- cds_exprs$f_id
        }
    }
    else {
        cds_exprs$feature_label <- cds_exprs$f_id
    }
    if (is.null(panel_order) == FALSE) {
        cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
            levels = panel_order)
    }
    q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
    if (is.null(fill_by) == FALSE) {
        q <- q + geom_violin(aes_string(fill = fill_by)) 
    }
    else {
        q <- q + geom_violin() 
    }
    if (plot_trend == TRUE) {
        q <- q + stat_summary(aes_string(color = color_by), fun.data = "mean_cl_boot", 
            size = trend_size)
        q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
            color = color_by, group = color_by), fun.data = "mean_cl_boot", 
            size = trend_size, geom = "line")
    }
    
    q <- q + facet_wrap(~feature_label, nrow = nrow, 
        ncol = ncol, scales = "free_y") 
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    
    if (log_scale == TRUE){

    q <- q + scale_y_log10()
   
    }
    
    q <- q + ylab("Expression") + xlab(grouping)
    q <- q +  monocle:::monocle_theme_opts()
    
    if(horizontal_orientation == TRUE){
    
        q <- q + coord_flip()
       
    }
    q
}
                         
auc_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)", 
    min_expr = NULL, df = 3, maxit = 300, relative_expr = TRUE, pseudocount = 0, cores=1) 
{
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)

        message("Getting expression matrix")
    
        cds_exprs <- Biobase::exprs(cds_subset)
        cds_exprs <- cds_exprs + pseudocount
        cds_exprs <- t(t(cds_exprs)/sizeFactors(cds_subset))

        cds_exprs <- melt(as.matrix(cds_exprs))
    
    colnames(cds_exprs) <- c("gene_id", "Cell", "expression")
    
        message("Merging data")
    
    cds_exprs <- merge(cds_exprs, cds_fData, by.x = "gene_id", 
        by.y = "row.names")
    cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    
    if (is.null(min_expr)) {
        min_expr <- cds_subset@lowerDetectionLimit
    }
    else{
        
        min_expr <- min_expr
    }

    cds_exprs$gene_label <- cds_exprs$gene_id

    cds_exprs$gene_label <- factor(cds_exprs$gene_label)
    new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    
    message("Fitting smooth curves")
    model_expectation <- genSmoothCurves(cds_subset, cores = cores, 
        trend_formula = trend_formula, relative_expr = T, 
        new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
        x[1]])
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr

    pseudotime_max = max(cds_exprs$Pseudotime)
    
    message("Calculating area under the curve")
        
    cds_exprs.auc = cds_exprs %>%
        dplyr::select(id, Pseudotime, expectation) %>%
        dplyr::mutate(pseudotime_quantile = cut(Pseudotime, breaks = seq(0, pseudotime_max, pseudotime_max/5))) %>%
        unique() %>%
        dplyr::group_by(id, pseudotime_quantile) %>%
        dplyr::summarize(AUC=MESS::auc(Pseudotime, expectation, type='spline')) %>%
        dplyr::group_by(id, pseudotime_quantile) %>% 
        dplyr::summarize(auc = abs(AUC)) %>%
        dplyr::arrange(desc(auc))

    cds_exprs.auc$pseudotime_quantile = as.numeric(factor(cds_exprs.auc$pseudotime_quantile, levels=sort(unique(cds_exprs.auc$pseudotime_quantile))))
    
    cds_exprs.auc_score = tidyr::spread(cds_exprs.auc, key=pseudotime_quantile, value=auc)
        
    return(as.data.frame(cds_exprs.auc_score))
        }
                                   

plot_pseudotime_heatmap_binSize <- function (cds_subset, cluster_rows = TRUE, bin_size = 100,  
                                             hclust_method = "ward.D2", 
    num_clusters = NA, hmcols = NULL, add_annotation_row = NULL, annotation_row = NA,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
    norm_method = c("vstExprs", "log"), scale_max = 3, scale_min = -3, 
    trend_formula = "~sm.ns(Pseudotime, df=3)", main = NULL, return_heatmap = FALSE, 
    cores = 1) 
{
    pseudocount <- NA
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
        max(pData(cds_subset)$Pseudotime), length.out = bin_size))
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
        relative_expr = T, new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
        FALSE) {
        
        ### Possibly offending line ###
        m = vstExprs(cds_subset, expr_matrix = m)
    }
    else if (norm_method == "log") {
        m = log10(m + pseudocount)
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- inferno(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE, 
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F, 
        clustering_distance_rows = row_dist, clustering_method = hclust_method, 
        cutree_rows = num_clusters, silent = TRUE, filename = NA, 
        breaks = bks, color = hmcols)

        annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
        num_clusters)))

    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), 
            ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
    }
    
    else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) <- row_ann_labels
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE, 
        cluster_rows = cluster_rows, show_rownames = show_rownames, 
        show_colnames = F, clustering_distance_rows = row_dist, 
        clustering_method = hclust_method, cutree_rows = num_clusters, 
        annotation_row = annotation_row, treeheight_row = 10, 
        breaks = bks, fontsize = 6, color = hmcols, silent = TRUE, 
        main = main, annotation_colors = ann_colors,
        filename = NA)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    
    print(names(cluster_rows))
    
    if (return_heatmap) {
        return(ph_res)
    }
}

collect_gsa_hyper_results_clusters <- function (genes_list, clusters, gsc) 
{
    gene_universe <- unique(as.character(genes_list))
    gsa_results <- list()
    cluster_ids <- unique(clusters)
    for (i in (1:length(cluster_ids))) {
        cluster_genes <- unique(names(clusters[clusters == i]))
        gsaRes <- runGSAhyper(cluster_genes, gsc = gsc, universe = gene_universe, 
            adjMethod = "BH")
        gsa_results[[length(gsa_results) + 1]] <- gsaRes
    }
    names(gsa_results) <- cluster_ids
    gsa_results
}
                                   
gsea_bar_plots <- function(GSAhyper_list, qval_cutoff, pattern, width, height, sample, gsc){
    
    for(cluster in names(GSAhyper_list)){
        
        print(cluster)
        
            GSAhyper_df <- as.data.frame(GSAhyper_list[[cluster]]$p.adj)
            GSAhyper_df$gene_set <- row.names(GSAhyper_df)
            colnames(GSAhyper_df) <- c("qval","gene_set")

            if(is.null(pattern) == FALSE){
                GSAhyper_df$gene_set <- str_replace(string = GSAhyper_df$gene_set, pattern = pattern, replace = "")
            }

            GSAhyper_df_cutoff <- GSAhyper_df %>% filter(qval < qval_cutoff) %>% arrange(desc(qval)) %>% 
            mutate(gene_set = factor(gene_set, levels = gene_set))

            plot_title <- paste0(sample,"_",as.character(cluster),"_",gsc,".png")
            print(plot_title)
            
            ggplot(GSAhyper_df_cutoff, aes(x = gene_set, y = -log10(qval))) + 
            geom_bar(stat = "identity") + 
            coord_flip() +
            theme_classic(base_size = 8) +
            ggsave(plot_title, width = width, height = height)
        
        
        
    }
    
}
                                   
getDTWcds <- function(query_cds, ref_cds, ref, query, expressed_genes, cores){

alignment_genes <- intersect(row.names(subset(fData(ref_cds), use_for_ordering)),
                         row.names(subset(fData(query_cds), use_for_ordering)))

ref_align_cds <- ref_cds[alignment_genes]
query_align_cds <- query_cds[alignment_genes]

### Set a consistent Pseudospace between both ordering sets
pData(ref_align_cds)$cell_id <- row.names(pData(ref_align_cds))
pData(ref_align_cds)$Pseudotime <- 100 * pData(ref_align_cds)$Pseudotime / max(pData(ref_align_cds)$Pseudotime)
ref_align_cds <- ref_align_cds[alignment_genes,as.character(arrange(pData(ref_align_cds), Pseudotime)$cell_id)]

pData(query_align_cds)$cell_id <- row.names(pData(query_align_cds))
pData(query_align_cds)$Pseudotime <- 100 * pData(query_align_cds)$Pseudotime / max(pData(query_align_cds)$Pseudotime)
query_align_cds <- query_align_cds[alignment_genes,as.character(arrange(pData(query_align_cds), Pseudotime)$cell_id)]

# Fits a smoothed curve to alignment genes accross Pseudotime
    #closeAllConnections()
smoothed_ref_exprs <- genSmoothCurves(ref_align_cds[alignment_genes], data.frame(Pseudotime=seq(0,100, by=1)), cores=cores)
smoothed_ref_exprs <- smoothed_ref_exprs[rowSums(is.na(smoothed_ref_exprs)) == 0,]

vst_smoothed_ref_exprs <- vstExprs(ref_cds, expr_matrix=smoothed_ref_exprs)

    #closeAllConnections()
smoothed_query_exprs <- genSmoothCurves(query_align_cds[alignment_genes], data.frame(Pseudotime=seq(0,100, by=1)), cores=cores)
smoothed_query_exprs <- smoothed_query_exprs[rowSums(is.na(smoothed_query_exprs)) == 0,]

vst_smoothed_query_exprs <- vstExprs(query_cds, expr_matrix=smoothed_query_exprs)

alignment_genes <- intersect(row.names(vst_smoothed_ref_exprs), row.names(vst_smoothed_query_exprs))

ref_matrix <- t(scale(t(vst_smoothed_ref_exprs[alignment_genes,])))
query_matrix <- t(scale(t(vst_smoothed_query_exprs[alignment_genes,])))

ref_query_dtw <- align_cells(ref_matrix, query_matrix, step_pattern=rabinerJuangStepPattern(3, "c"), open.begin=F, open.end=F)

align_res <- warp_pseudotime(ref_align_cds, query_align_cds, ref_query_dtw)
    
    
query_ref_aligned <- align_res$query_cds

pData(query_ref_aligned)$Pseudotime <- pData(query_ref_aligned)$Alignment_Pseudotime

ref_aligned_cell_ids <- setdiff(row.names(pData(ref_align_cds)), "duplicate_root")
query_aligned_cell_ids <- setdiff(row.names(pData(query_align_cds)), "duplicate_root")

combined_exprs <- cBind(exprs(query_cds[expressed_genes,query_aligned_cell_ids]), 
                        exprs(ref_cds[expressed_genes,ref_aligned_cell_ids]))

pData_ref <- pData(ref_align_cds)[,c("treatment_id", "spatial_id", "Pseudotime")]
pData_ref$Cell.Type <- ref


pData_query_aligned <- pData(query_ref_aligned)[,c("treatment_id", "spatial_id", "Pseudotime")]
pData_query_aligned$Cell.Type <- query


combined_pData <- rbind(pData_query_aligned, pData_ref)
combined_pData <- combined_pData[colnames(combined_exprs),]

combined_pd <- new("AnnotatedDataFrame", data = combined_pData)

fd <- new("AnnotatedDataFrame", data = fData(ref_cds)[row.names(combined_exprs),1:2])

ref_queryToRef_combined_cds <- newCellDataSet(combined_exprs, 
                 phenoData = combined_pd, 
                 featureData = fd,
                 expressionFamily= VGAM::negbinomial.size(), 
                 lowerDetectionLimit=1)

pData(ref_queryToRef_combined_cds)$cell_id <- row.names(pData(ref_queryToRef_combined_cds))

    return(ref_queryToRef_combined_cds)
    
    }

rank_auc_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
    min_expr = NULL, df = 3, maxit = 300, relative_expr = TRUE, pseudocount = 0, cores=1) 
{
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
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
        cds_exprs <- exprs(cds_subset)
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
        trend_formula = trend_formula, relative_expr = T, 
        new_data = new_data)
    colnames(model_expectation) <- colnames(cds_subset)
    cds_exprs$expectation <- apply(cds_exprs, 1, function(x) model_expectation[x[2], 
        x[1]])
    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
    cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr

    cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr

    pseudotime_max = max(cds_exprs$Pseudotime)

    cds_exprs.auc = cds_exprs %>%
        dplyr::select(id, Cell.Type, Pseudotime, expectation) %>%
        dplyr::mutate(pseudotime_quantile = cut(Pseudotime, breaks = seq(0, pseudotime_max, pseudotime_max/5))) %>%
        unique() %>%
        dplyr::group_by(id, Cell.Type, pseudotime_quantile) %>%
        dplyr::summarize(AUC=MESS::auc(Pseudotime, expectation, type='spline')) %>%
        dplyr::group_by(id, pseudotime_quantile) %>% 
        dplyr::summarize(auc_difference = abs(AUC[1] - AUC[2]) / (AUC[1] + AUC[2])) %>%
        dplyr::arrange(desc(auc_difference))

    cds_exprs.auc$pseudotime_quantile = as.numeric(factor(cds_exprs.auc$pseudotime_quantile, levels=sort(unique(cds_exprs.auc$pseudotime_quantile))))
    
    cds_exprs.auc_rank = tidyr::spread(cds_exprs.auc, key=pseudotime_quantile, value=auc_difference)
        
    return(as.data.frame(cds_exprs.auc_rank))
}                                   
        
                                   
compare_cell_types_in_pseudospace <- function (cds_subset, trend_formula = "~ sm.ns(Pseudotime, df=3)*Cell.Type", 
    min_expr = NULL, cell_size = 0.75, cell_alpha = 1, line_size = 1, nrow = NULL, 
    ncol = 1, panel_order = NULL, color_by = "Cell.Type", shade_by = NULL, 
    df = 3, maxit = 300, relative_expr = TRUE, pseudocount = 0) 
{
    cds_pData <- pData(cds_subset)
    cds_fData <- fData(cds_subset)
    if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
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
            cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
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
        trend_formula = trend_formula, relative_expr = T, 
        new_data = new_data)
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
    q <- q + geom_line(aes(Pseudotime, expectation, color = Cell.Type), size = line_size)
    q <- q + scale_y_log10() + facet_wrap(~gene_label, nrow = nrow, 
        ncol = ncol, scales = "free_y")
    q <- q + ylab("Expression") + xlab("Pseudospace")
    q <- q + monocle:::monocle_theme_opts()
    q
}


# qp <- compare_cell_types_in_pseudotime(combined_cds[row.names(subset(fData(combined_cds),
# 								 gene_short_name %in% missing_myo_genes
# 								 )),], color_by="Cell.Type", ncol=3, df=4) + scale_color_brewer(palette="Set1")
# pdf("missing_myogenesis_genes.pdf", width=6, height=8)
# print (qp)
# dev.off()
#

# We're gonna want to use this implementation of Mahalanobis distance for various clustering applications
fastPwMahal = function(x1) {
  require(corpcor)
  #SQRT = with(irlba(invCovMat, min(dc)-1, min(dc)-1), u %*% diag(d^0.5) %*% t(v))
  SQRT = with(fast.svd(ginv(cov(x1))), u %*% diag(d^0.5) %*% t(v))
  dist(x1 %*% SQRT)
}


loadGSCSafe <- function (file, type = "auto", addInfo, sep="\t", encoding="latin1") 
{
  if (missing(addInfo)) {
    addUserInfo <- "skip"
    addInfo <- "none"
  }
  else {
    addUserInfo <- "yes"
  }
  tmp <- try(type <- match.arg(type, c("auto", "gmt", "sbml", 
                                       "sif", "data.frame"), several.ok = FALSE), silent = TRUE)
  if (class(tmp) == "try-error") {
    stop("argument type set to unknown value")
  }
  if (type == "auto") {
    if (class(file) == "character") {
      tmp <- unlist(strsplit(file, "\\."))
      type <- tolower(tmp[length(tmp)])
      if (!type %in% c("gmt", "sif", "sbml", "xml")) 
        stop(paste("can not handle .", type, " file extension, read manually using e.g. read.delim() and load as data.frame", 
                   sep = ""))
    }
    else {
      type <- "data.frame"
    }
  }
  if (type == "gmt") {
    con <- file(file, encoding=encoding)
    tmp <- try(suppressWarnings(open(con)), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    if (addUserInfo == "skip") 
      addInfo <- vector()
    gscList <- list()
    i <- 1
    tmp <- try(suppressWarnings(while (length(l <- scan(con, 
                                                        nlines = 1, what = "character", quiet = T, sep=sep)) > 0) {
      if (addUserInfo == "skip") 
        addInfo <- rbind(addInfo, l[1:2])
      tmp <- l[3:length(l)]
      gscList[[l[1]]] <- unique(tmp[tmp != "" & tmp != 
                                      " " & !is.na(tmp)])
      i <- i + 1
    }), silent = TRUE)
    if (class(tmp) == "try-error") 
      stop("file could not be read")
    close(con)
    gsc <- gscList[!duplicated(names(gscList))]
    if (addUserInfo == "skip") 
      addInfo <- unique(addInfo)
  }
  else if (type %in% c("sbml", "xml")) {
    require(rsbml)
    tmp <- try(sbml <- rsbml_read(file))
    if (class(tmp) == "try-error") {
      stop("file could not be read by rsbml_read()")
    }
    gsc <- list()
    for (iReaction in 1:length(reactions(model(sbml)))) {
      metIDs <- names(c(reactants(reactions(model(sbml))[[iReaction]]), 
                        products(reactions(model(sbml))[[iReaction]])))
      geneIDs <- names(modifiers(reactions(model(sbml))[[iReaction]]))
      if (length(geneIDs) > 0) {
        geneNames <- rep(NA, length(geneIDs))
        for (iGene in 1:length(geneIDs)) {
          geneNames[iGene] <- name(species(model(sbml))[[geneIDs[iGene]]])
        }
        for (iMet in 1:length(metIDs)) {
          gsc[[metIDs[iMet]]] <- c(gsc[[metIDs[iMet]]], 
                                   geneNames)
        }
      }
    }
    if (length(gsc) == 0) {
      stop("no gene association found")
    }
    else {
      for (iMet in 1:length(gsc)) {
        tmp1 <- name(species(model(sbml))[[names(gsc)[iMet]]])
        tmp2 <- compartment(species(model(sbml))[[names(gsc)[iMet]]])
        names(gsc)[iMet] <- paste(tmp1, " (", tmp2, ")", 
                                  sep = "")
      }
    }
  }
  else if (type == "sif") {
    tmp <- try(gsc <- as.data.frame(read.delim(file, header = FALSE, 
                                               quote = "", as.is = TRUE), stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be read and converted into a data.frame")
    }
    if (ncol(gsc) != 3) {
      stop("sif file should contain three columns")
    }
    if (addUserInfo == "skip") 
      addInfo <- gsc[, c(1, 2)]
    gsc <- gsc[, c(3, 1)]
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  else if (type == "data.frame") {
    tmp <- try(gsc <- as.data.frame(file, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("argument file could not be converted into a data.frame")
    }
    for (i in 1:ncol(gsc)) {
      gsc[, i] <- as.character(gsc[, i])
    }
    if (ncol(gsc) != 2) {
      stop("argument file has to contain exactly two columns")
    }
    tmp <- nrow(gsc)
    gsc <- unique(gsc)
    geneSets <- unique(gsc[, 2])
    gscList <- list()
    for (iGeneSet in 1:length(geneSets)) {
      gscList[[iGeneSet]] <- gsc[gsc[, 2] == geneSets[iGeneSet], 
                                 1]
    }
    names(gscList) <- geneSets
    gsc <- gscList
  }
  if (addUserInfo == "yes") {
    tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors = FALSE), 
               silent = TRUE)
    if (class(tmp) == "try-error") {
      stop("failed to convert additional info in argument 'addInfo' into a data.frame")
    }
  }
  if (class(addInfo) == "data.frame") {
    if (ncol(addInfo) != 2) 
      stop("additional info in argument 'file' or 'addInfo' has to contain 2 columns")
    tmp <- nrow(addInfo)
    addInfo <- unique(addInfo[addInfo[, 1] %in% names(gsc), 
                              ])
  }
  else {
  }
  res <- list(gsc, addInfo)
  names(res) <- c("gsc", "addInfo")
  class(res) <- "GSC"
  return(res)
}



plot_gsa_hyper_heatmap <- function(cds, gsa_results, significance=0.05)
{
  hyper_df <- ldply(gsa_results, function(gsa_res)
  {
    data.frame(gene_set = names(gsa_res$pvalues), pval = gsa_res$pvalues, qval = gsa_res$p.adj)
  })
  colnames(hyper_df)[1] <- "cluster_id"
  
  #hyper_df 
  
  hyper_df <- subset(hyper_df, qval <= significance)
  print (head(hyper_df))
  hyper_df <- merge(hyper_df, ddply(hyper_df, .(gene_set), function(x) { nrow(x) }), by="gene_set")
  #print (hyper_df)
  hyper_df$gene_set <- factor(hyper_df$gene_set, levels=unique(arrange(hyper_df, V1, cluster_id)$gene_set))
  
  qplot(cluster_id, gene_set, fill=-log10(qval), geom="tile", data=hyper_df) + scale_fill_gradientn(colours=rainbow(7))
}



collect_gsa_hyper_results <- function(cds, gsc, clusters)
{
  gene_universe <- unique(as.character(fData(cds)$gene_short_name))
  gsa_results <- list()
  cluster_ids <- unique(clusters)
  for (i in (1:length(cluster_ids))) {
    cluster_genes <- unique(fData(cds[names(clusters[clusters == i]),])$gene_short_name)
    gsaRes <- runGSAhyper(cluster_genes, gsc=gsc, universe=gene_universe)
    gsa_results[[length(gsa_results) + 1]] <- gsaRes
  }
  names(gsa_results) <- cluster_ids
  gsa_results
}

filter_gsc <- function(gsc, cds, alias_table, allowed_genes, gsc_whitelist)
{
  gene_gsc_names_allowed <- sapply(names(gsc$gsc), function(gsc_name, cds, alias_table, allowed_genes) { 
    
    allowed_gene_symbols <- fData(cds[allowed_genes,])$gene_short_name
    
    gsc_matches_to_official_symbol <- sum(grepl(paste("^",gsc_name, sep=""), allowed_gene_symbols))
    gsc_matches_to_alias <- sum(grepl(paste("^", gsc_name, sep=""), subset(alias_table, symbol %in% allowed_gene_symbols)$alias_symbol))
    allowed_gsc_names <- gsc_matches_to_official_symbol + gsc_matches_to_alias > 0
    allowed_gsc_names
  },
  cds, alias_table, allowed_genes
  )
  gsc$gsc <- gsc$gsc[union(names(gene_gsc_names_allowed)[gene_gsc_names_allowed], gsc_whitelist)]
  gsc
}

clusterGenesRelative<-function(expr_matrix, k, method=function(x){as.dist((1 - cor(t(x)))/2)}, ...){
  n<-method(expr_matrix)
  
  clusters<-pam(n,k, ...)
  class(clusters)<-"list"
  clusters$exprs<-expr_matrix
  clusters
}

plot_clusters_relative<-function(cds, 
                                 clustering,
                                 drawSummary=TRUE, 
                                 sumFun=mean_cl_normal,
                                 ncol=NULL, 
                                 nrow=NULL, 
                                 row_samples=NULL, 
                                 callout_ids=NULL,
                                 conf_int=0.68,
                                 draw_individual_genes=FALSE){
  m <- as.data.frame(clustering$exprs)
  
  baseline_expr <- mean(rowMeans(m))
  
  m$ids <- rownames(clustering$exprs)
  if (is.null(clustering$labels) == FALSE)
  {
    m$cluster = factor(clustering$labels[clustering$clustering], levels = levels(clustering$labels))
  }else{
    m$cluster <- factor(clustering$clustering)
  }
  
  cluster_sizes <- as.data.frame(table(m$cluster))    
  
  cluster_sizes$Freq <- paste("(", cluster_sizes$Freq, ")")   
  facet_labels <- str_join(cluster_sizes$Var1, cluster_sizes$Freq, sep=" ")
  
  facet_wrap_labeller <- function(gg.plot,labels=NULL) {
    #works with R 3.0.1 and ggplot2 0.9.3.1
    require(gridExtra)
    
    g <- ggplotGrob(gg.plot)
    gg <- g$grobs      
    strips <- grep("strip_t", names(gg))
    
    for(ii in seq_along(labels))  {
      modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                         grep=TRUE, global=TRUE)
      gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
    }
    
    g$grobs <- gg
    class(g) = c("arrange", "ggplot",class(g)) 
    g
  }
  
  
  m.melt <- melt(m, id.vars = c("ids", "cluster"))
  
  m.melt <- merge(m.melt, pData(cds), by.x="variable", by.y="row.names")
  
  
  if (is.null(row_samples) == FALSE){
    m.melt <- m.melt[sample(nrow(m.melt), row_samples),]
  }
  
  c <- ggplot(m.melt) + facet_wrap("cluster", ncol=ncol, nrow=nrow, scales="free_y")
  #c <- c + stat_density2d(aes(x = Pseudotime, y = value), geom="polygon", fill="white", color="black", size=I(0.1)) + facet_wrap("cluster", ncol=ncol, nrow=nrow)
  
  if (draw_individual_genes){
    c <- c + geom_line(aes(x=Pseudotime, y=value, group=ids), alpha=I(0.1), size=I(0.5))
  }
  
  if (drawSummary) {
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1),
                          fun.data = sumFun, color = "red",
                          alpha = 0.2, size = 0.5, geom = "smooth")
    c <- c + stat_summary(aes(x = Pseudotime, y = value, group = 1), 
                          fun.data = "median_hilow", fill = "black", 
                          alpha = 0.2, size = 0.5, geom = "ribbon", conf.int=conf_int)
  }
  
  #cluster_medians <- subset(m.melt, ids %in% clustering$medoids)
  
  #c <- c + geom_line()
  #c <- c + geom_line(aes(x=Pseudotime, y=value), data=cluster_medians, color=I("red"))
  c <- c + scale_color_hue(l = 50, h.start = 200) + theme(axis.text.x = element_text(angle = 0, 
                                                                                     hjust = 0)) + xlab("Pseudo-time") + ylab("Expression")
  c <- c + theme(strip.background = element_rect(colour = 'white', fill = 'white')) + 
    theme(panel.border = element_blank()) +
    theme(legend.position="none") +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank())
  
  c <- c + geom_hline(yintercept=baseline_expr, color="steelblue")
  #   if (draw_cluster_size){
  #     cluster_sizes <- as.data.frame(table(m$cluster))
  #     colnames(cluster_sizes) <- c("cluster", "Freq")
  #     cluster_sizes <- cbind (cluster_sizes, Pseudotime = cluster_label_text_x, value = cluster_label_text_y)
  #     c <- c + geom_text(aes(x=Pseudotime, y=value, label=Freq), data=cluster_sizes, size=cluster_label_text_size)
  #   }
  
  if (is.null(callout_ids) == FALSE)
  {
    callout_melt <- subset(m.melt, ids %in% callout_ids)
    c <- c + geom_line(aes(x=Pseudotime, y=value), data=callout_melt, color=I("steelblue"))
  }
  #c <- c + monocle_theme_opts()
  #c <- facet_wrap_labeller(c, facet_labels)
  c
}

# selectHighDispersionGenes <- function(cds, 
#                                       detectionLimit=-Inf, 
#                                       relative_expr=TRUE){
#   disp_df<-esApply(cds,1,
#     function(f_expression) { 
#       if (relative_expr && cds@expressionFamily@vfamily == "negbinomial"){
#         f_expression <- f_expression / Size_Factor
#       }
#       f_expression <- f_expression[f_expression > detectionLimit]
#       expr_mean <- mean(f_expression[f_expression > 0])
#       if (is.null(expr_mean) == FALSE) {
#         disp_guess_fit <- cds@dispFitInfo[["blind"]]$disp_func(expr_mean)

#         # For NB: Var(Y)=mu*(1+mu/k)
#         f_expression_var <- var(f_expression)
#         f_expression_mean <- mean(f_expression)

#         disp_guess_meth_moments <- f_expression_var - f_expression_mean 
#         disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k 

#         return (data.frame(mean_exp=expr_mean, disp_fit=disp_guess_fit, disp_empirical=disp_guess_meth_moments))
#       }
#       return (NULL)
#   } )
#   do.call(rbind,disp_df)
# }

# NOTE: our definition of "reference" and "query" are reversed w.r.t what DTW uses.
align_cells <- function(query_cell_matrix, ref_cell_matrix, 
                        step_pattern=rabinerJuangStepPattern(3, "d", smoothed=TRUE),
                        open.end=FALSE,
                        open.begin=FALSE)
{
  # we want the query and reference to have the same number of cells, so downsample the larger
  # one to match the smaller one
  
  
  if (ncol(ref_cell_matrix) < ncol(query_cell_matrix)){
    ref_cell_matrix_sub <- ref_cell_matrix
    sampled_idxs <- c(1, sort(sample(seq(2,ncol(query_cell_matrix)-1), ncol(ref_cell_matrix)-2)), ncol(query_cell_matrix))
    query_cell_matrix_sub <- query_cell_matrix[,sampled_idxs]
  }else if (ncol(ref_cell_matrix) > ncol(query_cell_matrix)){
    query_cell_matrix_sub <- query_cell_matrix
    sampled_idxs <- c(1, sort(sample(seq(2,ncol(ref_cell_matrix)-1), ncol(query_cell_matrix)-2)), ncol(ref_cell_matrix))
    ref_cell_matrix_sub <- ref_cell_matrix[,sampled_idxs]
  }else{
    query_cell_matrix_sub <- query_cell_matrix
    ref_cell_matrix_sub <- ref_cell_matrix
  }
  
  ref_cell_matrix_sub <- t(scale(t(ref_cell_matrix_sub)))
  query_cell_matrix_sub <- t(scale(t(query_cell_matrix_sub)))
  
  valid_genes <- rowSums(is.na(ref_cell_matrix_sub)) == 0 & rowSums(is.na(query_cell_matrix_sub)) == 0
  ref_cell_matrix_sub <- ref_cell_matrix_sub[valid_genes,]
  query_cell_matrix_sub <- query_cell_matrix_sub[valid_genes,]
  
  print (dim(ref_cell_matrix_sub))
  #dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))^2
  dist_matrix <- (1 - cor(ref_cell_matrix_sub, query_cell_matrix_sub, use="pairwise.complete.obs"))
  
  cell_alignment_dtw <- dtw(dist_matrix,
                            step.pattern=step_pattern, 
                            open.end=open.end,
                            open.begin=open.begin,
                            keep=T,
                            keep.internals=T)
  cell_alignment_dtw
}


customDtwPlotDensity <- function (d, normalize = FALSE, xlab = "Query index", ylab = "Reference index", 
                                  ...) 
{
  cm <- d$costMatrix
  if (is.null(cm)) 
    stop("dtwPlotDensity requires dtw internals (set keep.internals=TRUE on dtw() call)")
  if (normalize) {
    norm <- attr(d$stepPattern, "norm")
    if (is.na(norm)) 
      stop("No normalization known for step pattern used")
    if (norm == "N") {
      cm <- cm/row(cm)
    }
    else if (norm == "N+M") {
      cm <- cm/(row(cm) + col(cm))
    }
    else if (norm == "M") {
      cm <- cm/col(cm)
    }
  }
  xd <- dim(cm)[1]
  yd <- dim(cm)[2]
  image(cm, col = terrain.colors(100), x = 1:xd, y = 1:yd, 
        xlab = xlab, ylab = ylab, useRaster=TRUE, ...)
  contour(cm, x = 1:xd, y = 1:yd, add = TRUE)
  lines(d$index1, d$index2, col = "cyan", lwd = 2)
}

warp_pseudotime <- function(ref_cds, query_cds, alignment_dtw)
{
 ref_pseudotimes <- pData(ref_cds)$Pseudotime
 names(ref_pseudotimes) <- row.names(pData(ref_cds))
 ref_pseudotimes <- sort(ref_pseudotimes)
 
 query_pseudotimes <- pData(query_cds)$Pseudotime
 names(query_pseudotimes) <- row.names(pData(query_cds))
 query_pseudotimes <- sort(query_pseudotimes)
 
 query_alignment_time <- warp(alignment_dtw, index.reference=T)
 
 #print (query_alignment_time)
 pData(ref_cds)$Alignment_Pseudotime <- pData(ref_cds)$Pseudotime
 pData(query_cds)$Alignment_Pseudotime <- NA
 
 #names(query_alignment_time) <- rownames(alignment_dtw$localCostMatrix)
 
 #pData(query_cds)[names(query_alignment_time),]$Alignment_Pseudotime <- query_alignment_time
 
 fc <- approxfun(seq(1, 101, by=1), query_alignment_time, rule=2)
 pData(query_cds)$Alignment_Pseudotime  <- fc(query_cds$Pseudotime)
 return(list(ref_cds=ref_cds, query_cds=query_cds))
 
 # BJ_MYO_aligned <- BJ_MYO_selected
 # pData(BJ_MYO_aligned)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
 # pData(BJ_MYO_selected)[colnames(smoothed_BJ_MYO_exprs),]$Pseudotime <- BJ_MYO_alignment_time
}

calABCs_from_alignment <- function(cds, 
                                   trend_formula = "~sm.ns(Pseudotime, df = 3)*Cell.Type",
                                   method = 'fitting', 
                                   ABC_method = 'integral', 
                                   points_num = 1000, 
                                   fc_limit = 3, 
                                   branchTest = FALSE, 
                                   cell_types = c("BJ-MYO", "HSMM"), 
                                   relative_expr = TRUE,
                                   stretch = TRUE,
                                   pseudocount=0,
                                   cores = 1, 
                                   weighted = TRUE, 
                                   min_expr = 0.5, 
                                   integer_expression = FALSE, 
                                   num = 100, 
                                   cell_type_labels = NULL, ...) {
  
  
  new_data = list()
  for (cell_type in unique(pData(cds)$Cell.Type)){
    new_data[[length(new_data) + 1]] = data.frame(Pseudotime=seq(1, 100, length.out=num),
                                                  Cell.Type=cell_type)
  }
  new_data = do.call(rbind, new_data)
  
  str_branchAB_expression_curve_matrix <- genSmoothCurves(cds, 
                                                          cores=cores, trend_formula = trend_formula,
                                                          relative_expr = relative_expr, 
                                                          pseudocount = pseudocount, 
                                                          new_data = new_data)
  
  str_branchA_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, 1:num]
  str_branchB_expression_curve_matrix <- str_branchAB_expression_curve_matrix[, (num + 1):(2 * num)]
  
  ABCs_res <- str_branchA_expression_curve_matrix - str_branchB_expression_curve_matrix
  ILR_res <- log2(str_branchA_expression_curve_matrix / (str_branchB_expression_curve_matrix + 0.1))
  
  ABCs_res <- apply(ABCs_res, 1, function(x, num, ABC_method) {
    avg_delta_x <- (x[1:(num - 1)] + x[2:(num)])/2
    step <- (100/(num - 1))
    
    if (ABC_method == "integral") {
      res <- round(sum(avg_delta_x * step), 3)
    }
    else if (ABC_method == "global_normalization") {
      max <- max(max(predictBranchOri), max(x))
      res <- round(sum(avg_delta_x/max * step), 3)
    }
    else if (ABC_method == "local_normalization") {
      pair_wise_max <- apply(data.frame(x = x, y = predictBranchOri),
                             1, max)
      res <- round(sum((((predictBranchOri - x)/pair_wise_max)[1:(num -
                                                                    1)] + ((predictBranchOri - x)/pair_wise_max)[2:(num)])/2 *
                         step), 3)
    }
    else if (ABC_method == "four_values") {
      ori_ABCs <- round(sum((x[1:(num - 1)] + x[2:(num)])/2 *
                              step), 3)
      other_ABCs <- round(sum((predictBranchOri[1:(num -
                                                     1)] + predictBranchOri[2:(num)])/2 * step),
                          3)
      ori_ABCs_H <- round(sum(avg_delta_x[avg_delta_x >
                                            0] * step), 3)
      other_ABCs_H <- round(sum(avg_delta_x[avg_delta_x <
                                              0] * step), 3)
      res <- c(ori_ABCs = ori_ABCs, other_ABCs = other_ABCs,
               ori_ABCs_H = ori_ABCs_H, other_ABCs_H = other_ABCs_H)
    }
    else if (ABC_method == "ILRs") {
      str_logfc_df <- log2((predictBranchOri + 1)/(x +
                                                     1))
      res <- sum(str_logfc_df)
    }
    return(res)}, num = num, ABC_method = ABC_method
  )
  
  ABCs_res <- cbind(ABCs_res, ILR_res[,ncol(ILR_res)])
  colnames(ABCs_res)<- c("ABCs", "Endpoint_ILR")
  
  ABCs_res <- merge(ABCs_res, fData(cds), by = "row.names")
  row.names(ABCs_res) <- ABCs_res[, 1]
  ABCs_res[, 1] <- NULL
  colnames(ABCs_res)[1] <- "ABCs"
  
  
  return(ABCs_res)
}


#The following code is swipped from colorRamps package which is used to make the pallette
table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
  x <- seq(0, 1, length.out = n)
  y <- rep(0, length(x))
  sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
  sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
  y[sill.min:sill.max] <- 1
  base.min <- round((n - 1) * (mid - base / 2)) + 1
  base.max <- round((n - 1) * (mid + base / 2)) + 1
  xi <- base.min:sill.min
  yi <- seq(0, 1, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  xi <- sill.max:base.max
  yi <- seq(1, 0, length.out = length(xi))
  i <- which(xi > 0 & xi <= n)
  y[xi[i]] <- yi[i]
  height * y
}


rgb.tables <- function(n,
                       red = c(0.75, 0.25, 1),
                       green = c(0.5, 0.25, 1),
                       blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

matlab.like <- function(n) rgb.tables(n)

matlab.like2 <- function(n)
  rgb.tables(n,
             red = c(0.8, 0.2, 1),
             green = c(0.5, 0.4, 0.8),
             blue = c(0.2, 0.2, 1))

blue2green2red <- matlab.like2

                                   
calculate_signature_scores <- function (cds, gene_list) 
{
    vst = as.data.frame(as.matrix(log2((Matrix::t(Matrix::t(Biobase:::exprs(cds))/pData(cds)$Size_Factor)) + 1)))

    vst = vst[gene_list$id, ]
    vst = merge(vst, gene_list[, names(gene_list) %in% c("id", 
        "signature")], by.x = "row.names", by.y = "id")
    vst <- melt(vst, id.vars = c("Row.names", "signature"), variable.name = "cell", 
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
                                   

                                   
                                   
                                   
                                   