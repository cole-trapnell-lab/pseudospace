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
library(tidyr)
library(matrixStats)
library(cluster)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(viridis)
library(ggrepel)})

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

preprocess_cds <- function(cds){
    
    cds <- detectGenes(cds, min_expr = 0.1)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    return(cds)
}

getPseudospaceTrajectory <- function(cds, sig_genes){
    
    cds <- setOrderingFilter(cds, sig_genes)
    cds <- reduceDimension(cds, max_components = 2, norm_method = "log")
    cds <- orderCells(cds, reverse = FALSE)
    
    return(cds)
}

## Need to update function in pseudospace_support_functions to specify which columns of pData to keep after alignment
getDTWcds <- function(query_cds, ref_cds, ref, query, expressed_genes, cores = 1){

alignment_genes <- intersect(row.names(subset(fData(ref_cds), use_for_ordering)),
                         row.names(subset(fData(query_cds), use_for_ordering)))

ref_align_cds <- ref_cds[alignment_genes]
query_align_cds <- query_cds[alignment_genes]

### Set a consistent Pseudospace between both ordering sets
    message("Normalizing pseudospace for each sample")
pData(ref_align_cds)$cell_id <- row.names(pData(ref_align_cds))
pData(ref_align_cds)$Pseudotime <- 100 * pData(ref_align_cds)$Pseudotime / max(pData(ref_align_cds)$Pseudotime)
ref_align_cds <- ref_align_cds[alignment_genes,as.character(arrange(pData(ref_align_cds), Pseudotime)$cell_id)]

pData(query_align_cds)$cell_id <- row.names(pData(query_align_cds))
pData(query_align_cds)$Pseudotime <- 100 * pData(query_align_cds)$Pseudotime / max(pData(query_align_cds)$Pseudotime)
query_align_cds <- query_align_cds[alignment_genes,as.character(arrange(pData(query_align_cds), Pseudotime)$cell_id)]

# Fits a smoothed curve to alignment genes accross Pseudotime
    message("Fitting smooth curves across pseudospace")
    #closeAllConnections()
smoothed_ref_exprs <- genSmoothCurves(ref_align_cds[alignment_genes], data.frame(Pseudotime=seq(0,100, by=1)), cores= cores)
smoothed_ref_exprs <- smoothed_ref_exprs[rowSums(is.na(smoothed_ref_exprs)) == 0,]

vst_smoothed_ref_exprs <- vstExprs(ref_cds, expr_matrix=smoothed_ref_exprs)

    #closeAllConnections()
smoothed_query_exprs <- genSmoothCurves(query_align_cds[alignment_genes], data.frame(Pseudotime=seq(0,100, by=1)), cores= cores)
smoothed_query_exprs <- smoothed_query_exprs[rowSums(is.na(smoothed_query_exprs)) == 0,]

vst_smoothed_query_exprs <- vstExprs(query_cds, expr_matrix=smoothed_query_exprs)

alignment_genes <- intersect(row.names(vst_smoothed_ref_exprs), row.names(vst_smoothed_query_exprs))

ref_matrix <- t(scale(t(vst_smoothed_ref_exprs[alignment_genes,])))
query_matrix <- t(scale(t(vst_smoothed_query_exprs[alignment_genes,])))

    message("Aligning pseudopsatial trajectories with dynamic time warping")
ref_query_dtw <- align_cells(ref_matrix, query_matrix, step_pattern=rabinerJuangStepPattern(3, "c"), open.begin=F, open.end=F)
    
    message("Warping pseudospace")
align_res <- warp_pseudotime(ref_align_cds, query_align_cds, ref_query_dtw)
    
    
query_ref_aligned <- align_res$query_cds

pData(query_ref_aligned)$Pseudotime <- pData(query_ref_aligned)$Alignment_Pseudotime

ref_aligned_cell_ids <- setdiff(row.names(pData(ref_align_cds)), "duplicate_root")
query_aligned_cell_ids <- setdiff(row.names(pData(query_align_cds)), "duplicate_root")

combined_exprs <- cBind(Biobase::exprs(query_cds[expressed_genes,query_aligned_cell_ids]), 
                        Biobase::exprs(ref_cds[expressed_genes,ref_aligned_cell_ids]))

pData_ref <- pData(ref_align_cds)[,c("gene","all_gene","barcode","proportion","guide_count","condition","treatment", "position", "Pseudotime")]
pData_ref$Cell.Type <- ref


pData_query_aligned <- pData(query_ref_aligned)[,c("gene","all_gene","barcode","proportion","guide_count","condition","treatment", "position", "Pseudotime")]
pData_query_aligned$Cell.Type <- query


combined_pData <- rbind(pData_query_aligned, pData_ref)
combined_pData <- combined_pData[colnames(combined_exprs),]

combined_pd <- new("AnnotatedDataFrame", data = combined_pData)

fd <- new("AnnotatedDataFrame", data = fData(ref_cds)[row.names(combined_exprs),1:2])

    message("Creating a new cds object with a common pseudospatial axes")
ref_queryToRef_combined_cds <- newCellDataSet(combined_exprs, 
                 phenoData = combined_pd, 
                 featureData = fd,
                 expressionFamily=negbinomial.size(), 
                 lowerDetectionLimit=1)

pData(ref_queryToRef_combined_cds)$cell_id <- row.names(pData(ref_queryToRef_combined_cds))

    return(ref_queryToRef_combined_cds)
    
}

# Expectation maximiation model t0 correct for different efficiencies across sgRNAs
get.guide.weights = function(mat, ntc.dist, n.iterations = 30) {
    n.guides = nrow(mat)
    n.cells = rowSums(mat)
    empirical.dist = sweep(mat, 1, n.cells, "/")
    
    lof.prop = rep(0.5, n.guides)
    expected.n.lof = n.cells * lof.prop
    
    for (i in 1:n.iterations) {
        lof.dist = sapply(1:n.guides, function(guide) {
            p = lof.prop[guide]
            (empirical.dist[guide,] - (1-p) * ntc.dist) / p
        })

        lof.dist = rowSums(sweep(lof.dist, 2, expected.n.lof / sum(expected.n.lof), "*"))
        lof.dist = ifelse(lof.dist < 0, 0, lof.dist)
        lof.dist = lof.dist / sum(lof.dist)

        lof.prop = sapply(1:n.guides, function(guide) {
            optimize(function(p) dmultinom(mat[guide,], prob = p * lof.dist + (1-p) * ntc.dist, log = T),
                c(0.0, 1.0), maximum = T)$maximum
        })

        expected.n.lof = n.cells * lof.prop
    }

    return(lof.prop)
}

calculate_ntc_empirical_fdr <- function(cds, iterations){

    chisq_qval.list <- list()
    
    median_NTC.list <- list()
    
        median_NTC.list[["Mock"]] <- median((pData(cds.aligned.list[["Mock"]]) %>% 
                                             group_by(gene) %>%
                                             summarize(n = n()))$n)
        median_NTC.list[["TGFB"]] <- median((pData(cds.aligned.list[["TGFB"]]) %>% 
                                             group_by(gene) %>% 
                                             summarize(n = n()))$n)
    
    NTC_cell_subset.list <- list()
    
        for(sample in names(cds.aligned.list)){
    
            NTC_cell_subset.list[[sample]] <- row.names(subset(pData(cds.aligned.list[[sample]]),
                                                               gene == "NONTARGETING"))

        }
    
    for(i in 1:iterations){
        
        if(i %in% c(1,10,100,250,500,750,100)){message(paste0("Iteration ", i," of ",as.character(iterations)))}
    
    cds_list <- cds
    
    random_NTC_subset.list <- list()
    
    for(sample in names(cds_list)){

        set.seed(i)
        random_NTC_subset.list[[sample]] <- sample(NTC_cell_subset.list[[sample]], 50, replace = FALSE)
        
    }
    
    new_gene_assignments.list <- list()
    
    for(sample in names(cds_list)){
    
        new_gene_assignments.list[[sample]] <- sapply(pData(cds_list[[sample]])$cell,function(x){
        
        if(x %in% random_NTC_subset.list[[sample]]) return("NTC_decoy")
        return(pData(cds_list[[sample]])[x,]$gene)
            
    })
    
    }
        pData(cds_list[["Mock"]])$gene <- new_gene_assignments.list[["Mock"]]
        pData(cds_list[["TGFB"]])$gene <- new_gene_assignments.list[["TGFB"]]
        
    analysis.targets = list()
    
    analysis.targets[["Mock"]] = as.data.frame(pData(cds_list[["Mock"]]) %>% filter(gene != "NONTARGETING") %>%
                                               group_by(gene) %>% summarize(
                                                   n.cells = n(),
                                                   n.guides = length(intersect(unique(barcode), analysis.guides[["Mock"]]))) %>%
                                               filter(n.cells >= 15, n.guides >= 1) %>% dplyr::select(gene))[,1]

    analysis.targets[["TGFB"]] = as.data.frame(pData(cds_list[["TGFB"]]) %>% filter(gene != "NONTARGETING") %>%
                                               group_by(gene) %>% summarize(
                                                   n.cells = n(),
                                                   n.guides = length(intersect(unique(barcode), analysis.guides[["TGFB"]]))) %>%
                                               filter(n.cells >= 15, n.guides >= 1) %>% dplyr::select(gene))[,1]
        
        
    target.to.guide.map <- list()
    
    for (target in analysis.targets[["Mock"]]) {
        target.to.guide.map[["Mock"]][[target]] = 
        sort(unique(as.data.frame(pData(cds_list[["Mock"]]) %>%
                                  filter(gene == target, barcode %in% analysis.guides[["Mock"]]) %>%
                                  dplyr::select(barcode))[, 1]))
    }

    for (target in analysis.targets[["TGFB"]]) {
        target.to.guide.map[["TGFB"]][[target]] = 
        sort(unique(as.data.frame(pData(cds_list[["TGFB"]]) %>%
                                  filter(gene == target, barcode %in% analysis.guides[["TGFB"]]) %>%
                                  dplyr::select(barcode))[, 1]))
    }
    
    guide.to.target.map = list()
    
    for(sample in names(cds_list)){
        guide.to.target.map[[sample]] = list()

        for (target in analysis.targets[[sample]]) {
            for (guide in target.to.guide.map[[sample]][[target]]) {
                guide.to.target.map[[sample]][[guide]] = target
            }
        }
    }
     
    target.region.mat = list()
    
    target.region.mat[["Mock"]] = acast(
        pData(cds_list[["Mock"]]) %>%
        filter(barcode %in% analysis.guides[["Mock"]] | gene == "NONTARGETING") %>%
        mutate(dummy = 1) %>% dplyr::select(gene, region, dummy),
        gene ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)

    target.region.mat[["TGFB"]] = acast(
        pData(cds_list[["TGFB"]]) %>%
        filter(barcode %in% analysis.guides[["TGFB"]] | gene == "NONTARGETING") %>%
        mutate(dummy = 1) %>% dplyr::select(gene, region, dummy),
        gene ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)
        
        NTC_decoy.region.mat <- list()
        NTC_decoy.region.mat[["Mock"]] <- matrix(target.region.mat[["Mock"]][row.names(target.region.mat[["Mock"]]) == "NTC_decoy",],
                           nrow = 1)
        row.names(NTC_decoy.region.mat[["Mock"]]) <- "NTC_decoy"
        NTC_decoy.region.mat[["TGFB"]] <- matrix(target.region.mat[["TGFB"]][row.names(target.region.mat[["TGFB"]]) == "NTC_decoy",],
                           nrow = 1)
        row.names(NTC_decoy.region.mat[["TGFB"]]) <- "NTC_decoy"
        
        weighted.target.region.mat[["Mock"]] <- rbind(weighted.target.region.mat[["Mock"]],
                                              NTC_decoy.region.mat[["Mock"]])
        weighted.target.region.mat[["TGFB"]] <- rbind(weighted.target.region.mat[["TGFB"]],
                                              NTC_decoy.region.mat[["TGFB"]])

    NTC.region.p = list()
    
    for(sample in names(cds_list)){
        
        pData(cds_list[[sample]])$gene <- as.factor(pData(cds_list[[sample]])$gene)
        pData(cds_list[[sample]])$region <- as.factor(pData(cds_list[[sample]])$region)
        
        NTC.region.p[[sample]] <- pData(cds_list[[sample]]) %>%
        group_by(gene, region) %>%
        summarize(n = n()) %>%
        tidyr::complete(region, fill = list(n = 0.1)) %>%
        filter(gene == "NONTARGETING")
        
        
    }
        
    
    ntc.distribution = list()

    for(sample in names(cds_list)){
        ntc.distribution[[sample]] = weighted.target.region.mat[[sample]]["NONTARGETING",]
        ntc.distribution[[sample]] = ntc.distribution[[sample]] / sum(ntc.distribution[[sample]])
    }
    
    initial.target.level.chisq.pval = list()
    
    for(sample in names(cds_list)){
        set.seed(42)
        initial.target.level.chisq.pval[[sample]] = sapply(
            analysis.targets[[sample]], function(target) {
                suppressWarnings({chisq.test(
                    weighted.target.region.mat[[sample]][target,],
                    p = NTC.region.p[[sample]]$n,
                    simulate.p.value = F, rescale.p = T, B = 1000)$p.value})
            })
    }
    
    initial.target.level.chisq.qval <- list()
    
    for(sample in names(cds_list)){
    
        initial.target.level.chisq.qval[[sample]] <- p.adjust(initial.target.level.chisq.pval[[sample]], method = "BH")
        initial.target.level.chisq.qval[[sample]] <- sapply(initial.target.level.chisq.qval[[sample]], 
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})
    }   
    
    chisq_qval.list[[i]] <- initial.target.level.chisq.qval

}
     return(chisq_qval.list)

}

#### Load data ####
Pseudospace_lof_cds <- readRDS("CROPseq_pseudospace_cds.rds")

# Create a cds subset for each stimulation condition that contains spatially isolated cells
cds.list <- list()

cds.list[["Mock"]] <- Pseudospace_lof_cds[,!is.na(pData(Pseudospace_lof_cds)$proportion) & 
                                     pData(Pseudospace_lof_cds)$guide_count == 1 &
                                     pData(Pseudospace_lof_cds)$treatment == "mock"]

cds.list[["TGFB"]] <- Pseudospace_lof_cds[,!is.na(pData(Pseudospace_lof_cds)$proportion) &
                                     pData(Pseudospace_lof_cds)$guide_count == 1 &
                                     pData(Pseudospace_lof_cds)$treatment == "tgfb"]

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

for(sample in names(cds.list)) {
    
    cds.list[[sample]] <- preprocess_cds(cds.list[[sample]])
}

# Identify genes that vary significantly between inner and outer CROPseq cell fractions

Spatial.DEG.test.list <- list()

for(sample in names(cds.list)){
    
    Spatial.DEG.test.list[[sample]] <- differentialGeneTest(cds.list[[sample]][expressed_genes.list[[sample]]], 
                                         fullModelFormulaStr = "~position",
                                         reducedModelFormulaStr = "~1",
                                         cores = 1)
}

# Calculate fold change in expression levels of significant genes between CROPseq cell fractions isolated by space
for(sample in names(Spatial.DEG.test.list)){
    
    diff_test_genes <- row.names(Spatial.DEG.test.list[[sample]])
    diff_cds <- cds.list[[sample]][diff_test_genes]
    diff_FC <- diff_foldChange(diff_cds, "position","inner")
    Spatial.DEG.test.list[[sample]]$log2_foldChange <- diff_FC$log2FC_outer
    
    rm(diff_test_genes,diff_cds,diff_FC)
}

Spatial_sig_genes.list <- list()

for(sample in names(Spatial.DEG.test.list)){
    
    Spatial_sig_genes.list[[sample]] <- row.names(subset(Spatial.DEG.test.list[[sample]], qval <= 1e-6 &
                                                        abs(log2_foldChange) >= 1))
    print(length(Spatial_sig_genes.list[[sample]]))
}

# Create pseudospatial trajectories and examine the distribution of inner and outer cells within them

for(sample in names(cds.list)){
    
    cds.list[[sample]] <- getPseudospaceTrajectory(cds.list[[sample]],
                                                   Spatial_sig_genes.list[[sample]])
}

cds.list[["Mock"]] <- orderCells(cds.list[["Mock"]], reverse = F)
cds.list[["TGFB"]] <- orderCells(cds.list[["TGFB"]], reverse = F)


plot_cell_trajectory(cds.list[["Mock"]], color_by = "position",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=20), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#0075F2", "#D62828"), 
                   name = "Spatial Context") +
ggsave(file = "MCF10A_Mock_loss_of_function_PseudospatialTrajectory.png", height = 6, width = 6)

ggplot(pData(cds.list[["Mock"]]), aes(x = Pseudotime, fill = position, color = position)) + 
  geom_density() + 
  facet_wrap(~position, ncol = 1) +
  theme_classic() +
  scale_color_manual("Spatial Context", labels = c("inner colony", "outer colony"),
                     values = c("#000000","#000000")) +
  scale_fill_manual("Spatial Context", labels = c("inner colony", "outer colony")
                    , values = c("#0075F2","#D62828")) +
  xlab("Pseudospace") +
  ylab("Cell density") +
  monocle:::monocle_theme_opts() +
  theme(legend.position = "top", legend.direction = "vertical",  text=element_text(size=20)) +
 ggsave("Mock_loss_of_function_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

plot_cell_trajectory(cds.list[["TGFB"]], color_by = "position",show_branch_points = FALSE) + 
theme(legend.position="top", text=element_text(size=20), legend.direction = "vertical") +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#70163C", "#38726C"), 
                   name = "Spatial Context") +
ggsave(file = "MCF10A_TGFB_loss_of_function_PseudospatialTrajectory.png", height = 6, width = 6)

ggplot(pData(cds.list[["TGFB"]]), aes(x = Pseudotime, fill = position, color = position)) + 
  geom_density() + 
  facet_wrap(~position, ncol = 1) +
  theme_classic() +
  scale_color_manual("Spatial Context", labels = c("inner colony", "outer colony"),
                     values = c("#000000","#000000")) +
  scale_fill_manual("Spatial Context", labels = c("inner colony", "outer colony")
                    , values = c("#70163C", "#38726C")) +
  xlab("Pseudospace") +
  ylab("Cell density") +
  monocle:::monocle_theme_opts() +
  theme(legend.position = "top", legend.direction = "vertical",  text=element_text(size=20)) +
  ggsave("TGFB_loss_of_function_cell_density_accross_pseudospace_geom_density.png", height = 6, width = 5)

plot_genes_in_pseudotime(cds.list[["Mock"]][fData(cds.list[["Mock"]])$gene_short_name == "CDH1",], 
                         color_by = "position", min_expr = 0.1) +
theme(text=element_text(size=20)) +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#0075F2", "#D62828"), 
                   name = "Spatial Context")

plot_genes_in_pseudotime(cds.list[["TGFB"]][fData(cds.list[["TGFB"]])$gene_short_name == "CDH1",], 
                         color_by = "position", min_expr = 0.1)+
theme(text=element_text(size=20)) +
scale_color_manual(labels = c("inner colony", "outer colony"), values = c("#38726C", "#70163C"), 
                   name = "Spatial Context")

expressed_genes <- unique(union(expressed_genes.list[["Mock"]],expressed_genes.list[["TGFB"]]))

# Plot the expression of known EMT markers across pseudospace
Mock_Figure1_Mar <- cds.list[["Mock"]][row.names(subset(fData(cds.list[["Mock"]]), gene_short_name %in% 
  c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]
plot_genes_in_pseudotime(Mock_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#0075F2","outer"="#D62828")) + 
ggsave("MCF10A_Mock_CFG_Figure1Markers_byPseudospace.png", width = 6, height =5)

TGFB_Figure1_Mar <- cds.list[["TGFB"]][row.names(subset(fData(cds.list[["TGFB"]]), gene_short_name %in% 
  c("CDH1","CRB3","DSP", "CDH2","FN1","VIM"))),]
plot_genes_in_pseudotime(TGFB_Figure1_Mar, color_by = "spatial_id", ncol = 2, min_expr = 0.1,
                        panel_order = c("CDH1","CDH2","CRB3","FN1","DSP","VIM")) + 
xlab("Pseudospace") +                
theme(legend.position = "none",text=element_text(size=20)) + 
scale_color_manual(values = c("inner" = "#70163C","outer"="#38726C")) + 
ggsave("MCF10A_TGFB_CFG_Figure1Markers_byPseudospace.png", width = 6, height =5)

# Use dynamic time warping to align Mock and TGFB pseudospatial trajectories and create a cds object of aligned trajectories
TGFB.to.Mock.CFG.aligned.cds <- getDTWcds(cds.list[["TGFB"]],cds.list[["Mock"]], 
                                      ref = "Mock", query = "TGFB", 
                                      expressed_genes = expressed_genes, cores = 1)

TGFB.to.Mock.CFG.aligned.cds <- estimateSizeFactors(TGFB.to.Mock.CFG.aligned.cds)

# Divide the aligned cds by treatment to test accumulation of knockouts along pseudospace independently
cds.aligned.list <- list()

cds.aligned.list[["Mock"]] <- TGFB.to.Mock.CFG.aligned.cds[,pData(TGFB.to.Mock.CFG.aligned.cds)$Cell.Type == "Mock"]
cds.aligned.list[["TGFB"]] <- TGFB.to.Mock.CFG.aligned.cds[,pData(TGFB.to.Mock.CFG.aligned.cds)$Cell.Type == "TGFB"]

for(sample in names(cds.aligned.list)){
    
    cds.aligned.list[[sample]] <- preprocess_cds(cds.aligned.list[[sample]])
}

for(sample in names(cds.aligned.list)){
    
    cds.aligned.list[[sample]]@reducedDimA <- t(as.matrix(pData(cds.aligned.list[[sample]])$Pseudotime))
    colnames(cds.aligned.list[[sample]]@reducedDimA) <- row.names(pData(cds.aligned.list[[sample]]))
    
}

for(sample in names(cds.aligned.list)){
    
    cds.aligned.list[[sample]] <- clusterCells(cds.aligned.list[[sample]], method = "densityPeak")
    
}

plot_rho_delta(cds.aligned.list[["Mock"]],rho_threshold = 50, delta_threshold = 5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
plot_rho_delta(cds.aligned.list[["TGFB"]], rho_threshold = 10, delta_threshold = 5) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 50)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

rho_delta.list <- list()

rho_delta.list[["Mock"]] <- c(50,5)
rho_delta.list[["TGFB"]] <- c(10,5)

for(sample in names(cds.aligned.list)){
    
    cds.aligned.list[[sample]] <- clusterCells(cds.aligned.list[[sample]], method = "densityPeak",
                                              verbose = T,
                                               rho_threshold = rho_delta.list[[sample]][1], 
                                               delta_threshold = rho_delta.list[[sample]][2],
                                               skip_rho_sigma = T)
    
}

ggplot(pData(cds.aligned.list[["Mock"]]), aes(x = Pseudotime, fill = Cluster)) + 
geom_density() + monocle:::monocle_theme_opts()

ggplot(pData(cds.aligned.list[["TGFB"]]), aes(x = Pseudotime, fill = Cluster)) + 
geom_density() + monocle:::monocle_theme_opts()

# Re-order regions to be in order from low to high pseudospace
region.list <- list()

region.list[["Mock"]] <- sapply(pData(cds.aligned.list[["Mock"]])$Cluster, function(x){
    
    if(x == "3")return("1")
    if(x == "6")return("2")
    if(x == "7")return("3")
    if(x == "1")return("4")
    if(x == "4")return("5")
    if(x == "5")return("6")
    if(x == "2")return("7")
    
})

region.list[["TGFB"]] <- sapply(pData(cds.aligned.list[["TGFB"]])$Cluster, function(x){
    
    if(x == "3")return("1")
    if(x == "7")return("2")
    if(x == "8")return("3")
    if(x == "6")return("4")
    if(x == "4")return("5")
    if(x == "1")return("6")
    if(x == "5")return("7")
    if(x == "2")return("8")
    
})

for(sample in names(cds.aligned.list)){
    
    pData(cds.aligned.list[[sample]])$region <- region.list[[sample]]
    
}

ggplot(pData(cds.aligned.list[["Mock"]]), aes(x = Pseudotime, fill = region)) + 
geom_density() + monocle:::monocle_theme_opts()

ggplot(pData(cds.aligned.list[["TGFB"]]), aes(x = Pseudotime, fill = region)) + 
geom_density() + monocle:::monocle_theme_opts()

mock_regions <- sapply(row.names(pData(cds.list[["Mock"]])), function(x){
    
    return(pData(cds.aligned.list[["Mock"]])[x,]$region)
    
})

tgfb_regions <- sapply(row.names(pData(cds.list[["TGFB"]])), function(x){
    
    return(pData(cds.aligned.list[["TGFB"]])[x,]$region)
    
})

pData(cds.list[["Mock"]])$region <- mock_regions
pData(cds.list[["TGFB"]])$region <- tgfb_regions

analysis.guides = list()

    for(sample in names(cds.aligned.list)){
        pData(cds.aligned.list[[sample]]) %>% filter(guide_count == 1) %>% group_by(gene, barcode) %>%
        summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
        filter(n.guide.cells >= 10) %>% ungroup() %>%
        arrange(-n.target.cells, -n.guide.cells) %>% head(10)

        analysis.guides[[sample]] = 
        (pData(cds.aligned.list[[sample]]) %>% filter(guide_count == 1) %>% group_by(gene, barcode) %>%
        summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
        filter(n.guide.cells >= 10) %>% ungroup())$barcode

    }

analysis.targets = list()

analysis.targets[["Mock"]] = as.data.frame(pData(cds.aligned.list[["Mock"]]) %>% filter(gene != "NONTARGETING") %>%
                                               group_by(gene) %>% summarize(
                                                   n.cells = n(),
                                                   n.guides = length(intersect(unique(barcode), analysis.guides[["Mock"]]))) %>%
                                               filter(n.cells >= 15, n.guides >= 1) %>% dplyr::select(gene))[,1]

analysis.targets[["TGFB"]] = as.data.frame(pData(cds.aligned.list[["TGFB"]]) %>% filter(gene != "NONTARGETING") %>%
                                               group_by(gene) %>% summarize(
                                                   n.cells = n(),
                                                   n.guides = length(intersect(unique(barcode), analysis.guides[["TGFB"]]))) %>%
                                               filter(n.cells >= 15, n.guides >= 1) %>% dplyr::select(gene))[,1]

target.region.mat = list()
    
    target.region.mat[["Mock"]] = acast(
        pData(cds.aligned.list[["Mock"]]) %>%
        filter(barcode %in% analysis.guides[["Mock"]] | gene == "NONTARGETING") %>%
        mutate(dummy = 1) %>% dplyr::select(gene, region, dummy),
        gene ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)

    target.region.mat[["TGFB"]] = acast(
        pData(cds.aligned.list[["TGFB"]]) %>%
        filter(barcode %in% analysis.guides[["TGFB"]] | gene == "NONTARGETING") %>%
        mutate(dummy = 1) %>% dplyr::select(gene, region, dummy),
        gene ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)

target.to.guide.map <- list()
    
    for (target in analysis.targets[["Mock"]]) {
        target.to.guide.map[["Mock"]][[target]] = 
        sort(unique(as.data.frame(pData(cds.aligned.list[["Mock"]]) %>%
                                  filter(gene == target, barcode %in% analysis.guides[["Mock"]]) %>%
                                  dplyr::select(barcode))[, 1]))
    }

    for (target in analysis.targets[["TGFB"]]) {
        target.to.guide.map[["TGFB"]][[target]] = 
        sort(unique(as.data.frame(pData(cds.aligned.list[["TGFB"]]) %>%
                                  filter(gene == target, barcode %in% analysis.guides[["TGFB"]]) %>%
                                  dplyr::select(barcode))[, 1]))
    }

NTC.guides <- unique(pData(cds.aligned.list[["Mock"]])[pData(cds.aligned.list[["Mock"]])$gene == "NONTARGETING",]$barcode)


guide.region.mat = list()

guide.region.mat[["Mock"]] = acast(
    pData(cds.aligned.list[["Mock"]]) %>% filter(barcode %in% analysis.guides[["Mock"]][!(analysis.guides[["Mock"]] %in% NTC.guides)]) %>%
    mutate(dummy = 1) %>% dplyr::select(barcode, region, dummy),
    barcode ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)

guide.region.mat[["TGFB"]] = acast(
    pData(cds.aligned.list[["TGFB"]]) %>% filter(barcode %in% analysis.guides[["TGFB"]][!(analysis.guides[["TGFB"]] %in% NTC.guides)]) %>%
    mutate(dummy = 1) %>% dplyr::select(barcode, region, dummy),
    barcode ~ region, value.var = "dummy", fun.aggregate = sum, fill = 0)

ntc.distribution = list()

    for(sample in names(cds.aligned.list)){
        
        ntc.distribution[[sample]] = target.region.mat[[sample]]["NONTARGETING",]
        ntc.distribution[[sample]] = ntc.distribution[[sample]] / sum(ntc.distribution[[sample]])
    }

weighted.target.region.mat <- list()

for (condition in c("Mock", "TGFB")) {
    weighted.target.region.mat[[condition]] = t(sapply(analysis.targets[[condition]],
            function(target) {
        guides = target.to.guide.map[[condition]][[target]]
        if (length(guides) == 1) {
            return(target.region.mat[[condition]][target,])
        } else {
            mat = guide.region.mat[[condition]][guides,]
            guide.weights = get.guide.weights(mat, ntc.distribution[[condition]])
            guide.weights = guide.weights / max(guide.weights)
            print(condition)
            print(target)
        
            print(round(guide.weights, 3))
            #return(target.cluster.mat[[condition]][target,])
            return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
        }
    }))
}

NTC.region.mat <- list()
NTC.region.mat[["Mock"]] <- matrix(target.region.mat[["Mock"]][row.names(target.region.mat[["Mock"]]) == "NONTARGETING",],
                           nrow = 1)
row.names(NTC.region.mat[["Mock"]]) <- "NONTARGETING"
NTC.region.mat[["TGFB"]] <- matrix(target.region.mat[["TGFB"]][row.names(target.region.mat[["TGFB"]]) == "NONTARGETING",],
                           nrow = 1)
row.names(NTC.region.mat[["TGFB"]]) <- "NONTARGETING"

weighted.target.region.mat[["Mock"]] <- rbind(weighted.target.region.mat[["Mock"]],
                                              NTC.region.mat[["Mock"]])
weighted.target.region.mat[["TGFB"]] <- rbind(weighted.target.region.mat[["TGFB"]],
                                              NTC.region.mat[["TGFB"]])

median((pData(cds.aligned.list[["TGFB"]]) %>% 
                                         group_by(gene) %>% 
                                         summarize(n = n()))$n)

row.names(weighted.target.region.mat[["Mock"]])

# Calculate empirical FDR
chisq_qval.list <- calculate_ntc_empirical_fdr(cds.aligned.list, 1000)

mock_chisq_pval.list <- lapply(chisq_qval.list, `[[`, 1)
tgfb_chisq_pval.list <- lapply(chisq_qval.list, `[[`, 2)

mock_chisq_pval_df <- do.call("rbind",mock_chisq_pval.list)
mock_chisq_pval_df <- t(mock_chisq_pval_df)
tgfb_chisq_pval_df <- do.call("rbind",tgfb_chisq_pval.list)
tgfb_chisq_pval_df <- t(tgfb_chisq_pval_df)

mock_chisq_pval_df_test <- as.data.frame(mock_chisq_pval_df)
tgfb_chisq_pval_df_test <- as.data.frame(tgfb_chisq_pval_df)

mock_chisq_pval_NTC_df <- t(mock_chisq_pval_df_test["NTC_decoy",])
colnames(mock_chisq_pval_NTC_df) <- "mock_NTC_decoy"
tgfb_chisq_pval_NTC_df <- t(tgfb_chisq_pval_df_test["NTC_decoy",])
colnames(tgfb_chisq_pval_NTC_df) <- "tgfb_NTC_decoy"

mock_chisq_pval_df_test["NTC_decoy",]

chisq_pval_NTC_df <- merge(mock_chisq_pval_NTC_df,tgfb_chisq_pval_NTC_df, by = "row.names")

met <- tgfb_chisq_pval_df[22,]

length(met[met < 0.005])

ggplot(chisq_pval_NTC_df, aes(x = mock_NTC_decoy)) + 
geom_histogram() +
xlim(0,1)

ggplot(chisq_pval_NTC_df, aes(x = tgfb_NTC_decoy)) + 
geom_histogram() +
xlim(0,1)

mock_empirical_FDR_df <- apply(mock_chisq_pval_df,2,function(x){
  x >= x["NTC_decoy"]
})

tgfb_empirical_FDR_df <- apply(tgfb_chisq_pval_df,2,function(x){
  x >= x["NTC_decoy"]
})

mock_empirical_FDR <- rowSums(mock_empirical_FDR_df)/1000
tgfb_empirical_FDR <- rowSums(tgfb_empirical_FDR_df)/1000

empirical_FDR_df <- as.data.frame(cbind(names(mock_empirical_FDR),mock_empirical_FDR,tgfb_empirical_FDR))
empirical_FDR_df <- melt(empirical_FDR_df)
colnames(empirical_FDR_df) <- c("target", "mock_FDR","TGFB_FDR")

empirical_FDR_df$target <- as.character(empirical_FDR_df$target)
empirical_FDR_df$mock_FDR <- as.numeric(as.character(empirical_FDR_df$mock_FDR))
empirical_FDR_df$TGFB_FDR <- as.numeric(as.character(empirical_FDR_df$TGFB_FDR))

# Plot FDR by KO at the target and individual guide levels
empirical_FDR_df <- empirical_FDR_df[order(empirical_FDR_df$mock_FDR, decreasing = FALSE),]
empirical_FDR_df$target <- factor(empirical_FDR_df$target, levels = empirical_FDR_df$target)

ggplot(empirical_FDR_df,aes( x = as.factor(target), y = mock_FDR, fill = mock_FDR < 0.1)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "dimgrey") +
  xlab("Target") +
  ylab("FDR\n(empirically determined)") +
  ggtitle("Spontaneous EMT") +
  scale_fill_manual("FDR < 0.1", values = c("TRUE" = "red","FALSE" = "black")) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  monocle:::monocle_theme_opts() +
  ggsave("Mock_empirical_target_level_FDR.png", height = 2, width = 4, units = "in")

empirical_FDR_df <- empirical_FDR_df[order(empirical_FDR_df$TGFB_FDR, decreasing = FALSE),]
empirical_FDR_df$target <- factor(empirical_FDR_df$target, levels = empirical_FDR_df$target)

ggplot(empirical_FDR_df,aes( x = as.factor(target), y = TGFB_FDR, fill = TGFB_FDR < 0.1)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "dimgrey") +
  xlab("Target") +
  ylab("FDR\n(empirically determined)") +
  ggtitle("TGF-B-driven EMT") +
  scale_fill_manual("FDR < 0.1",values = c("TRUE" = "red","FALSE" = "black")) +
  theme(text = element_text(size = 6), 
        axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  monocle:::monocle_theme_opts() +
  ggsave("TGFB_empirical_target_level_FDR.png", height = 2, width = 4, units = "in")

NTC.region.p = list()

for(sample in names(cds.aligned.list)){
NTC.region.p[[sample]] <- pData(cds.aligned.list[[sample]]) %>% 
    group_by(region) %>% 
    filter(gene == "NONTARGETING") %>%
    summarize(n = n()) %>% 
    complete(region, fill = list(n = 0.1))
}

pass.target.level.screen = list()

pass.target.level.screen[["Mock"]] <- as.character(empirical_FDR_df[empirical_FDR_df$mock_FDR < 0.1,]$target)
pass.target.level.screen[["TGFB"]] <- as.character(empirical_FDR_df[empirical_FDR_df$TGFB_FDR < 0.1,]$target)

print(pass.target.level.screen[["Mock"]])
print(pass.target.level.screen[["TGFB"]])

targets.passing.initial.screen = list()

for(sample in names(cds.aligned.list)){
targets.passing.initial.screen[[sample]] = pass.target.level.screen[[sample]]
}

length(targets.passing.initial.screen[["Mock"]])
length(targets.passing.initial.screen[["TGFB"]])

length(intersect(targets.passing.initial.screen[["Mock"]],
                 targets.passing.initial.screen[["TGFB"]]))

length(unique(union(targets.passing.initial.screen[["Mock"]],
                 targets.passing.initial.screen[["TGFB"]])))

# Re-weigh guides for calculating log2 odds ratios and for plotting enrichemnt heatmaps
weighted.target.region.mat = list()

for (condition in c("Mock", "TGFB")) {
    weighted.target.region.mat[[condition]] = t(sapply(targets.passing.initial.screen[[condition]],
            function(target) {
        guides = target.to.guide.map[[condition]][[target]]
        if (length(guides) == 1) {
            return(target.region.mat[[condition]][target,])
        } else {
            mat = guide.region.mat[[condition]][guides,]
            guide.weights = get.guide.weights(mat, ntc.distribution[[condition]])
            guide.weights = guide.weights / max(guide.weights)
            print(condition)
            print(target)
        
            print(round(guide.weights, 3))
            #return(target.cluster.mat[[condition]][target,])
            return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
        }
    }))
}

region.enrichment.df = list()

for (condition in c("Mock", "TGFB")) {
    weighted.mat = weighted.target.region.mat[[condition]]
    ntc.counts = target.region.mat[[condition]]["NONTARGETING",]

    region.enrichment.df[[condition]] = do.call(rbind, lapply(rownames(weighted.mat), function(target) {
        do.call(rbind, lapply(1:ncol(weighted.mat), function(region) {
            test = fisher.test(cbind(
                c(weighted.mat[target, region], sum(weighted.mat[target, -region])),
                c(ntc.counts[region], sum(ntc.counts[-region]))))

            data.frame(
                target = target,
                region = region,
                odds.ratio = unname(test$estimate),
                p.value = test$p.value)
        }))
    }))

   region.enrichment.df[[condition]]$q.value = p.adjust(region.enrichment.df[[condition]]$p.value, "BH")

    region.enrichment.df[[condition]]$log2.odds = with(region.enrichment.df[[condition]],
        ifelse(odds.ratio == 0, -5, round(log2(odds.ratio),2)))
}

region.enrichment.df[["Mock"]]$target <- as.character(region.enrichment.df[["Mock"]]$target)
region.enrichment.df[["TGFB"]]$target <- as.character(region.enrichment.df[["TGFB"]]$target)

region.enrichment.heatmap.df <- region.enrichment.df

for(sample in names(region.enrichment.heatmap.df)){

    region.enrichment.heatmap.df[[sample]]$log2.odds[region.enrichment.heatmap.df[[sample]]$log2.odds < -2] <- -2
    region.enrichment.heatmap.df[[sample]]$log2.odds[region.enrichment.heatmap.df[[sample]]$log2.odds > 2] <- 2
    
}

region.enrichment.heatmap.matrix <- list()

for(sample in names(region.enrichment.heatmap.df)){
    
    region.enrichment.heatmap.matrix[[sample]] <- recast(region.enrichment.heatmap.df[[sample]], 
                                                            target ~ region, measure.var = "log2.odds")
    row.names(region.enrichment.heatmap.matrix[[sample]]) <- region.enrichment.heatmap.matrix[[sample]]$target
    region.enrichment.heatmap.matrix[[sample]] <- region.enrichment.heatmap.matrix[[sample]][,-1]
    
}

# Plot enrichment heatmaps
pheatmap(region.enrichment.heatmap.matrix[["Mock"]],
         clustering_method = "ward.D2",
         show_rownames = T,
         show_colanmes = T,
         na_col = "grey90",
         col = colorspace::diverge_hsv(30),
         cluster_col= FALSE,
         cluster_row = TRUE,
         useRaster = TRUE,
         width = 4,
         height = 6,
         file = "Mock_region_enrichment_heatmap.png")

pheatmap(region.enrichment.heatmap.matrix[["TGFB"]],
         clustering_method = "ward.D2",
         show_rownames = T,
         show_colanmes = T,
         na_col = "grey90",
         col = colorspace::diverge_hsv(30),
         cluster_col= FALSE,
         cluster_row = TRUE,
         useRaster = TRUE,
         width = 4,
         height = 6,
        file = "TGFB_region_enrichment_heatmap.png")

# Highlight examples of accumulation across pseudospace for EGFR and MET in spontaneous EMT and TGFBRs in TGF-B driven EMT
region_3_min_pseudospace_value <- min(pData(cds.aligned.list[["Mock"]])[pData(cds.aligned.list[["Mock"]])$region == "3",]$Pseudotime)
region_3_max_pseudospace_value <- max(pData(cds.aligned.list[["Mock"]])[pData(cds.aligned.list[["Mock"]])$region == "3",]$Pseudotime)


ggplot(pData(cds.aligned.list[["Mock"]])[pData(cds.aligned.list[["Mock"]])$gene %in% c("NONTARGETING", "EGFR", "MET"),], 
       aes(x = Pseudotime, fill = gene)) + 
geom_rect(xmin = region_3_min_pseudospace_value, xmax = region_3_max_pseudospace_value, ymin = 0, ymax = Inf, 
          fill = "slategray1", alpha = 0.01) +
geom_density() +
facet_wrap(~factor(gene, levels = c("EGFR", "MET","NONTARGETING")), scales = "free_y", ncol = 1) +
theme(legend.position = "none", strip.text.x=element_text(size=18), 
      axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
      axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) +
scale_fill_manual(values = c("EGFR" = "firebrick3", "MET" = "brown4", "TGFBR2" = "navy", "NONTARGETING" = "dimgrey")) +
monocle:::monocle_theme_opts() +
xlab("Pseudospace") +
ggsave("Density_of_EGFR_MET_NTC_across_spontaneous_EMT.png", width = 5, height = 10)


region_4_max_pseudospace_value <- max(pData(cds.aligned.list[["TGFB"]])[pData(cds.aligned.list[["TGFB"]])$region == "4",]$Pseudotime)


ggplot(pData(cds.aligned.list[["TGFB"]])[pData(cds.aligned.list[["TGFB"]])$gene %in% c("NONTARGETING", "TGFBR1", "TGFBR2"),], 
       aes(x = Pseudotime, fill = gene)) + 
geom_rect(xmin = 0, xmax = region_4_max_pseudospace_value, ymin = 0, ymax = Inf, 
          fill = "slategray1", alpha = 0.01) +
geom_density() +
facet_wrap(~factor(gene, levels = c("TGFBR2", "TGFBR1", "NONTARGETING")), scales = "free_y", ncol = 1) +
theme(legend.position = "none", strip.text.x=element_text(size=18), 
      axis.text.x = element_text(size=10), axis.text.y = element_text(size=10),
      axis.title.x = element_text(size=18), axis.title.y = element_text(size=18)) +
scale_fill_manual(values = c("TGFBR2" = "firebrick3", "TGFBR1" = "brown4", "ITGAV" = "navy", "NONTARGETING" = "dimgrey")) +
monocle:::monocle_theme_opts() +
xlab("Pseudospace") +
ggsave("Density_of_TGFBR1_TGFBR2_NTC_across_TGFB_driven_EMT.png", width = 5, height = 10)

# FOr supplemental
# Determine the fraction of CDH1 single positive, CDH1/VIM double positive and VIM single positive cells upon confluence dependent EMT
mock_CDH1_VIM_cds_subset <- cds.aligned.list[["Mock"]][fData(cds.aligned.list[["Mock"]])$gene_short_name %in% c("CDH1","VIM"),
                                                  !is.na(pData(cds.aligned.list[["Mock"]])$region)]
mock_CDH1_VIM_cds_exprs <- Biobase::exprs(mock_CDH1_VIM_cds_subset)
mock_CDH1_VIM_cds_exprs <- Matrix::t(Matrix::t(mock_CDH1_VIM_cds_exprs)/sizeFactors(mock_CDH1_VIM_cds_subset))

tgfb_CDH1_VIM_cds_subset <- cds.aligned.list[["TGFB"]][fData(cds.aligned.list[["TGFB"]])$gene_short_name %in% c("CDH1","VIM"),
                                                  !is.na(pData(cds.aligned.list[["TGFB"]])$region)]
tgfb_CDH1_VIM_cds_exprs <- Biobase::exprs(tgfb_CDH1_VIM_cds_subset)
tgfb_CDH1_VIM_cds_exprs <- Matrix::t(Matrix::t(tgfb_CDH1_VIM_cds_exprs)/sizeFactors(tgfb_CDH1_VIM_cds_subset))

CDH1_expression_cutoff <- mean(mock_CDH1_VIM_cds_exprs[2,])
VIM_expression_cutoff <- mean(mock_CDH1_VIM_cds_exprs[1,])
CDH1_expression_cutoff
VIM_expression_cutoff

mock_CDH1_VIM_double_positive_cells <- colnames(mock_CDH1_VIM_cds_exprs[,mock_CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & 
                                                              mock_CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff])
mock_CDH1_positive_cells <- colnames(mock_CDH1_VIM_cds_exprs[,mock_CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff & 
                                                   !(colnames(mock_CDH1_VIM_cds_exprs) %in% mock_CDH1_VIM_double_positive_cells)])
mock_VIM_positive_cells <- colnames(mock_CDH1_VIM_cds_exprs[,mock_CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & 
                                                  !(colnames(mock_CDH1_VIM_cds_exprs) %in% mock_CDH1_VIM_double_positive_cells)])

tgfb_CDH1_VIM_double_positive_cells <- colnames(tgfb_CDH1_VIM_cds_exprs[,tgfb_CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & 
                                                              tgfb_CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff])
tgfb_CDH1_positive_cells <- colnames(tgfb_CDH1_VIM_cds_exprs[,tgfb_CDH1_VIM_cds_exprs[2,] > CDH1_expression_cutoff & 
                                                   !(colnames(tgfb_CDH1_VIM_cds_exprs) %in% tgfb_CDH1_VIM_double_positive_cells)])
tgfb_VIM_positive_cells <- colnames(tgfb_CDH1_VIM_cds_exprs[,tgfb_CDH1_VIM_cds_exprs[1,] > VIM_expression_cutoff & 
                                                  !(colnames(tgfb_CDH1_VIM_cds_exprs) %in% tgfb_CDH1_VIM_double_positive_cells)])

mock_pData <- pData(cds.aligned.list[["Mock"]][,!is.na(pData(cds.aligned.list[["Mock"]])$region)])
mock_pData$positive_marker <- sapply(mock_pData$cell, function(x){
    
    if(x %in% mock_CDH1_VIM_double_positive_cells){
        return("CDH1/VIM double positive")
    }
    
    if(x %in% mock_CDH1_positive_cells){
        return("CDH1 single positive")
    }
    
    if(x %in% mock_VIM_positive_cells){
        return("VIM single positive")
    }
    
    return(NA)
    
})

tgfb_pData <- pData(cds.aligned.list[["TGFB"]][,!is.na(pData(cds.aligned.list[["TGFB"]])$region)])
tgfb_pData$positive_marker <- sapply(tgfb_pData$cell, function(x){
    
    if(x %in% tgfb_CDH1_VIM_double_positive_cells){
        return("CDH1/VIM double positive")
    }
    
    if(x %in% tgfb_CDH1_positive_cells){
        return("CDH1 single positive")
    }
    
    if(x %in% tgfb_VIM_positive_cells){
        return("VIM single positive")
    }
    
    return(NA)
    
})

ggplot(pData(cds.aligned.list[["Mock"]])[mock_CDH1_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1 single positive cells") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("Mock_loss_of_function_cell_density_accross_pseudospace_CDH1_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.aligned.list[["Mock"]])[mock_CDH1_VIM_double_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1/VIM double positive cells") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("Mock_loss_of_function_cell_density_accross_pseudospace_CDH1_VIM_double_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.aligned.list[["Mock"]])[mock_VIM_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("VIM single positive cells") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("Mock_loss_of_function_cell_density_accross_pseudospace_VIM_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.aligned.list[["TGFB"]])[tgfb_CDH1_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1 single positive cells") +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("TGFB_loss_of_function_cell_density_accross_pseudospace_CDH1_positive_cells_geom_density.png", height = 3, width = 6)


ggplot(pData(cds.aligned.list[["TGFB"]])[tgfb_CDH1_VIM_double_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("CDH1/VIM double positive cells") +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("TGFB_loss_of_function_cell_density_accross_pseudospace_CDH1_VIM_double_positive_cells_geom_density.png", height = 3, width = 6)

ggplot(pData(cds.aligned.list[["TGFB"]])[tgfb_VIM_positive_cells,], aes(x = Pseudotime)) + 
  geom_density(fill = "gray70", color = "black") + 
  theme_classic() +
  monocle:::monocle_theme_opts() +
  xlab("Pseudospace") +
  ylab("Cell density") +
  ggtitle("VIM single positive cells") +
  theme(legend.position = "none",  
        axis.title.y = element_text(size = 14), axis.title.x = element_text(size = 24),
       axis.text.y = element_text(size = 14), plot.title = element_text(size = 24, hjust = 0.5)) +
  ggsave("TGFB_loss_of_function_cell_density_accross_pseudospace_VIM_positive_cells_geom_density.png", height = 3, width = 6)


# Identify differentially expressed genes across Pseudospace for every ko vs NTC combination
mock_target_NTC_diff_test.list <- list()

for(target in analysis.targets[["Mock"]]){
    
    message("Obtaining differentially expressed genes between ",target, " and NTC")
    subset_cds <- cds.aligned.list[["Mock"]][,pData(cds.aligned.list[["Mock"]])$gene == target | 
                                               pData(cds.aligned.list[["Mock"]])$gene == "NONTARGETING"]
    subset_cds <- estimateDispersions(subset_cds)
    mock_target_NTC_diff_test.list[[target]] <- myDifferentialGeneTest(subset_cds[expressed_genes.list[["Mock"]]],
                                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)+gene",
                                                                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                                                                       cores = 1)
    rm(subset_cds)
    message("Done")
    
}

tgfb_target_NTC_diff_test.list <- list()


for(target in analysis.targets[["TGFB"]]){
    
    message("Obtaining differentially expressed genes between ",target, " and NTC")
    subset_cds <- cds.aligned.list[["TGFB"]][,pData(cds.aligned.list[["TGFB"]])$gene == target |
                                            pData(cds.aligned.list[["TGFB"]])$gene == "NONTARGETING"]
    subset_cds <- estimateDispersions(subset_cds)
    tgfb_target_NTC_diff_test.list[[target]] <- myDifferentialGeneTest(subset_cds[expressed_genes.list[["TGFB"]]],
                                                                      fullModelFormulaStr = "~sm.ns(Pseudotime, df=3)+gene",
                                                                       reducedModelFormulaStr = "~sm.ns(Pseudotime, df=3)",
                                                                       cores = 1)
    rm(subset_cds)
    message("Done")
    
}


for(target in names(mock_target_NTC_diff_test.list)){
    
    mock_target_NTC_diff_test.list[[target]]$target <- rep(target, 
                                                           length(row.names(mock_target_NTC_diff_test.list[[target]])))
}

for(target in names(tgfb_target_NTC_diff_test.list)){
    
    tgfb_target_NTC_diff_test.list[[target]]$target <- rep(target, 
                                                           length(row.names(tgfb_target_NTC_diff_test.list[[target]])))
}

mock_target_NTC_diff_test <- do.call("rbind", mock_target_NTC_diff_test.list)

tgfb_target_NTC_diff_test <- do.call("rbind", tgfb_target_NTC_diff_test.list)

mock_target_NTC_diff_test$qval <- p.adjust(mock_target_NTC_diff_test$pval, method = "BH")
tgfb_target_NTC_diff_test$qval <- p.adjust(tgfb_target_NTC_diff_test$pval, method = "BH")

mock_target_NTC_sig_genes <- unique(subset(mock_target_NTC_diff_test, qval < 0.05)$id)
length(mock_target_NTC_sig_genes)
tgfb_target_NTC_sig_genes <- unique(subset(tgfb_target_NTC_diff_test, qval < 0.05)$id)
length(tgfb_target_NTC_sig_genes)

mock_target_NTC_dif_test_sig_subset <- mock_target_NTC_diff_test[mock_target_NTC_diff_test$id %in% mock_target_NTC_sig_genes,]
tgfb_target_NTC_dif_test_sig_subset <- tgfb_target_NTC_diff_test[tgfb_target_NTC_diff_test$id %in% tgfb_target_NTC_sig_genes,]


mock_diff_test_summary <- mock_target_NTC_diff_test %>% filter(qval < 0.05) %>% group_by(target) %>% 
summarize(n = n()) %>% arrange(desc(n))
colnames(mock_diff_test_summary) <- c("target","total_degs")

mock_cell_number_summary <- pData(cds.aligned.list[["Mock"]]) %>% filter(gene != "NONTARGETING") %>% 
group_by(gene) %>% summarize(n = n()) %>% arrange(desc(n))

tgfb_diff_test_summary <- tgfb_target_NTC_diff_test %>% filter(qval < 0.05) %>% group_by(target) %>% 
summarize(n = n()) %>% arrange(desc(n))
colnames(tgfb_diff_test_summary) <- c("target","total_degs")

tgfb_cell_number_summary <- pData(cds.aligned.list[["TGFB"]]) %>% filter(gene != "NONTARGETING") %>% 
group_by(gene) %>% summarize(n = n()) %>% arrange(desc(n))


mock_summary_df <- merge(mock_diff_test_summary, mock_cell_number_summary, by.x = "target", by.y = "gene")
tgfb_summary_df <- merge(tgfb_diff_test_summary, tgfb_cell_number_summary, by.x = "target", by.y = "gene")

# Plot number of DEGs for every KO vs number of KO cells in the experiment
ggplot(mock_summary_df, aes(x = n, y = total_degs, label = target)) + 
geom_point() +
monocle:::monocle_theme_opts() +
theme(text = element_text(size = 24)) +
xlab("Total number of target cells") +
ylab("Total number of DEGs") +
ggsave("Mock_CROPseq_TargetvsNTC_degs.png", width = 6, height = 6)

ggplot(tgfb_summary_df, aes(x = n, y = total_degs, label = target)) + 
geom_point() +
monocle:::monocle_theme_opts() +
theme(text = element_text(size = 24)) +
xlab("Total number of target cells") +
ylab("Total number of DEGs") +
ggsave("TGFB_CROPseq_TargetvsNTC_degs.png", width = 6, height = 6)

ggplot(tgfb_summary_df[!(tgfb_summary_df$target %in% c("TGFBR1","TGFBR2")),], aes(x = n, y = total_degs, label = target)) + 
geom_point() +
monocle:::monocle_theme_opts() +
theme(text = element_text(size = 24)) +
xlab("Total number of target cells") +
ylab("Total number of DEGs") +
ggsave("TGFB_CROPseq_TargetvsNTC_degs_woTGFRBs.png", width = 6, height = 6)


