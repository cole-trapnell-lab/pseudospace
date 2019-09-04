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
library(FNN)
library(monocle)})

sessionInfo()

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

knnProject <- function (cds, max_components = 2,  norm_method = c("log", 
    "vstExprs", "none"), verbose = FALSE, scaling = TRUE, knn_k = 10, 
                        pseudo_expr = 1, num_dim = 25, ...) 
{
    
    message("Scaling expression data")
    extra_arguments <- list(...)
    set.seed(2016)
    FM <- monocle:::normalize_expr_data(cds, norm_method, pseudo_expr)
    xm <- Matrix::rowMeans(FM)
    xsd <- sqrt(Matrix::rowMeans((FM - xm)^2))
    FM <- FM[xsd > 0, ]
    
    if (scaling) {
        FM <- as.matrix(Matrix::t(scale(Matrix::t(FM))))
        FM <- FM[!is.na(row.names(FM)), ]
    }
    else FM <- as.matrix(FM)
    if (nrow(FM) == 0) {
        stop("Error: all rows have standard deviation zero")
    }
    FM <- FM[apply(FM, 1, function(x) all(is.finite(x))), ]
                   
    message("Performing dimensionality reduction with PCA")
            
            FM <- (FM)
            irlba_res <- prcomp_irlba(t(FM), n = min(num_dim, 
                min(dim(FM)) - 1), center = TRUE, scale. = TRUE)
            irlba_pca_res <- irlba_res$x
            
    row.names(irlba_pca_res) <- colnames(FM)
                   
    MCF10A_cells <- row.names(pData(cds)[pData(cds)$Cell.Type == "MCF10A",])
    MCF10A.mtx <- irlba_pca_res[MCF10A_cells,]
                   
    knn.list <- list()
                   
    knn.list[["MCF10A_matrix"]] <- MCF10A.mtx
                   
    for(patient in unique(pData(cds)[pData(cds)$Cell.Type != "MCF10A",]$Cell.Type)){ 
    
    HNSCC_cells <- row.names(pData(cds)[pData(cds)$Cell.Type == patient,])
    HNSCC.mtx <- irlba_pca_res[HNSCC_cells,]
                   
    message(paste0("Identifying k nearest neighbors for tumor sample ",patient))
                
    knn.list[[patient]] = FNN::get.knnx(MCF10A.mtx, HNSCC.mtx, k=knn_k, algorithm="kd_tree")
    row.names(knn.list[[patient]]$nn.index) <- row.names(HNSCC.mtx) 
    row.names(knn.list[[patient]]$nn.dist) <- row.names(HNSCC.mtx) 

        
        }
           
    return(knn.list)
}

# Load Mock and TGFB cds objects created in Figure1 code
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

# Load cds object containing HNSCC single-cell data from Puram et al 
HNSCC_cds <- readRDS("HSNCC_cds.rds")

# Filter HNSCC data for samples with more than 40, non-lymph node, cancer cells processesd with Maxima enzyme
HNSCC_patient_list <- as.data.frame(pData(HNSCC_cds) %>% 
filter(Maxima_enzyme == "0" & clasified_as_cancer_cell == "1" &
      lymph_node == "0") %>%
group_by(patient_id) %>% summarize(n = n()) %>% arrange(desc(n)))

HNSCC_patients_forAnalysis <- HNSCC_patient_list[HNSCC_patient_list$n > 40,]$patient_id

HNSCC_patient_list

HNSCC_cds <- HNSCC_cds[,pData(HNSCC_cds)$patient_id %in% HNSCC_patients_forAnalysis & 
                       pData(HNSCC_cds)$Maxima_enzyme == "0" & 
                       pData(HNSCC_cds)$clasified_as_cancer_cell == "1" &
                      pData(HNSCC_cds)$lymph_node == "0"]

# Generate a new cds objects containing MCF10A mock or tgfb exposed cells and HNSCC tumor data
HNSCC_fData <- fData(HNSCC_cds)
MCF10A_mock_fData <- fData(cds.list[["Mock"]])[,c("id","gene_short_name")]
MCF10A_tgfb_fData <- fData(cds.list[["TGFB"]])[,c("id","gene_short_name")]

pData(HNSCC_cds)$cell <- row.names(pData(HNSCC_cds))

HNSCC_pData <- pData(HNSCC_cds)[,c("cell","patient_id")]
MCF10A_mock_pData <- pData(cds.list[["Mock"]])[,c("cell","sample","Pseudotime")]
MCF10A_tgfb_pData <- pData(cds.list[["TGFB"]])[,c("cell","sample","Pseudotime")]

colnames(HNSCC_pData) <- c("cell","sample")

HNSCC_pData$Pseudotime <- rep(NA, length(row.names(HNSCC_pData)))

HNSCC_exprs <- as(round(Biobase::exprs(HNSCC_cds)),"dgTMatrix")
MCF10A_mock_exprs <- Biobase::exprs(cds.list[["Mock"]])
MCF10A_tgfb_exprs <- Biobase::exprs(cds.list[["TGFB"]])

MCF10A_mock_exprs <- MCF10A_mock_exprs[row.names(HNSCC_exprs),]
MCF10A_tgfb_exprs <- MCF10A_tgfb_exprs[row.names(HNSCC_exprs),]

MCF10A_mock_fData <- MCF10A_mock_fData[row.names(HNSCC_fData),]
MCF10A_tgfb_fData <- MCF10A_tgfb_fData[row.names(HNSCC_fData),]

mock_exprs <- cbind(MCF10A_mock_exprs,HNSCC_exprs)
tgfb_exprs <- cbind(MCF10A_tgfb_exprs,HNSCC_exprs)

mock_pData <- rbind(MCF10A_mock_pData, HNSCC_pData)
tgfb_pData <- rbind(MCF10A_tgfb_pData, HNSCC_pData)

mock_pData <- new("AnnotatedDataFrame", data = mock_pData)
tgfb_pData <- new("AnnotatedDataFrame", data = tgfb_pData)

mock_fData <- MCF10A_mock_fData
tgfb_fData <- MCF10A_tgfb_fData

mock_fData <- new("AnnotatedDataFrame", data = mock_fData)
tgfb_fData <- new("AnnotatedDataFrame", data = tgfb_fData)

HNSCC_MCF10A_mock_cds <- newCellDataSet(as(as.matrix(mock_exprs), "sparseMatrix"),
                       phenoData = mock_pData,
                       featureData = mock_fData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

HNSCC_MCF10A_tgfb_cds <- newCellDataSet(as(as.matrix(tgfb_exprs), "sparseMatrix"),
                       phenoData = tgfb_pData,
                       featureData = tgfb_fData,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

# Add a pData column that specifies whether cells are from MCF10A or HNSCC tumor samples
mock_cell_type <- sapply(pData(HNSCC_MCF10A_mock_cds)$sample,function(x){
    
    if(x %in% c('Mock_inner','Mock_outer'))return("MCF10A")
    return(x)
    
})

tgfb_cell_type <- sapply(pData(HNSCC_MCF10A_tgfb_cds)$sample,function(x){
    
    if(x %in% c('TGFB_inner','TGFB_outer'))return("MCF10A")
    return(x)
    
})

pData(HNSCC_MCF10A_mock_cds)$Cell.Type <- mock_cell_type
pData(HNSCC_MCF10A_tgfb_cds)$Cell.Type <- tgfb_cell_type

# Estimate size factors for normalization
HNSCC_MCF10A_mock_cds <- estimateSizeFactors(HNSCC_MCF10A_mock_cds)
HNSCC_MCF10A_tgfb_cds <- estimateSizeFactors(HNSCC_MCF10A_tgfb_cds)

# Identify genes expressed in at least 50 cells in each dataset
mock_expressed_genes <- row.names(fData(HNSCC_MCF10A_mock_cds)[Matrix::rowSums(Biobase::exprs(HNSCC_MCF10A_mock_cds) > 0) > 50 ,])
length(mock_expressed_genes)

tgfb_expressed_genes <- row.names(fData(HNSCC_MCF10A_tgfb_cds)[Matrix::rowSums(Biobase::exprs(HNSCC_MCF10A_tgfb_cds) > 0) > 50 ,])
length(tgfb_expressed_genes)


# Identify the 20 nearest neighbors of every HNSCC tumors cell to MCF10A  
mock_MCF10AtoHNSCC_knn <- knnProject(HNSCC_MCF10A_mock_cds[mock_expressed_genes,], 
                                     knn_k = 20, num_dim = 20)


tgfb_MCF10AtoHNSCC_knn <- knnProject(HNSCC_MCF10A_tgfb_cds[tgfb_expressed_genes,], 
                                     knn_k = 20, num_dim = 20)


# Identify the average pseudospace of the 20 MCF10A nearest neighbors of every HNSCC cell
knn_mock_pseudospace.list <- list()

for(patient in HNSCC_patients_forAnalysis){
    
     knn_pseudospace <- apply(mock_MCF10AtoHNSCC_knn[[patient]]$nn.index,1,function(x){
    
        nearest_cells <- row.names(mock_MCF10AtoHNSCC_knn[["MCF10A_matrix"]][as.vector(x),])
        average_knn_pseudotime <- mean(pData(HNSCC_MCF10A_mock_cds)[nearest_cells,]$Pseudotime)
        
        })
    print(patient)
    print(length(knn_pseudospace))
    patient_id <- rep(patient, length(knn_pseudospace))
    knn_mock_pseudospace.list[[patient]] <- data.frame(Pseudospace = knn_pseudospace,
                                                      patient_id = patient_id)
    }                               

knn_tgfb_pseudospace.list <- list()

for(patient in HNSCC_patients_forAnalysis){
    
     knn_pseudospace <- apply(tgfb_MCF10AtoHNSCC_knn[[patient]]$nn.index,1,function(x){
    
        nearest_cells <- row.names(tgfb_MCF10AtoHNSCC_knn[["MCF10A_matrix"]][as.vector(x),])
        average_knn_pseudotime <- mean(pData(HNSCC_MCF10A_tgfb_cds)[nearest_cells,]$Pseudotime)
        
        })
    print(patient)
    print(length(knn_pseudospace))
    patient_id <- rep(patient, length(knn_pseudospace))
    knn_tgfb_pseudospace.list[[patient]] <- data.frame(Pseudospace = knn_pseudospace,
                                                      patient_id = patient_id)
    }                               

# Collapse into a single data.frame
knn_mock_pseudopace <- do.call("rbind", knn_mock_pseudospace.list)
knn_tgfb_pseudopace <- do.call("rbind", knn_tgfb_pseudospace.list)

# Order HNSCC tumors according to the partial EMT rank indentified in Puram et al. 
knn_mock_pseudopace$patient_id <- factor(knn_mock_pseudopace$patient_id, 
                                         levels = c("HNSCC6","HNSCC20","HNSCC5",
                                                    "HNSCC18","HNSCC22","HNSCC25","HNSCC17",
                                                    "HNSCC16"))

knn_tgfb_pseudopace$patient_id <- factor(knn_tgfb_pseudopace$patient_id, 
                                         levels = c("HNSCC6","HNSCC20","HNSCC5",
                                                    "HNSCC18","HNSCC22","HNSCC25","HNSCC17",
                                                    "HNSCC16"))

# Generate plots of the density of projected HSNCC cells across pseudospace
ggplot(knn_mock_pseudopace, aes(x = Pseudospace, fill = ..density..)) + 
geom_density(fill = "#0075F2") +
facet_wrap(~patient_id, ncol = 1, scales = "free_y") +
xlab("Spontaneous EMT") +
ylab("Nearest Neighbor density across Pseudospace") +
theme(text = element_text(size = 11),
     axis.text.y = element_text(size = 3)) +
monocle:::monocle_theme_opts() +
ggsave("Figure2i_KNN_projection_HNSCC_to_MCF10A_mock_treated_cells.png", width = 2, height = 4)

ggplot(knn_tgfb_pseudopace, aes(x = Pseudospace, fill = ..density..)) + 
geom_density(fill = "#70163C") +
facet_wrap(~patient_id, ncol = 1, scales = "free_y") +
xlab("TGF-B-driven EMT") +
ylab("Nearest Neighbor density across Pseudospace") +
theme(text = element_text(size = 11),
     axis.text.y = element_text(size = 3)) +
monocle:::monocle_theme_opts() +
ggsave("Figure2i_KNN_projection_HNSCC_to_MCF10A_tgfb_treated_cells.png", width = 2, height = 4)


