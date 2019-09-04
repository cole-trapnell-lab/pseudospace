# Load necessary packages 
library(reshape2)
library(dplyr)
library(monocle)

sessionInfo()

# Load single-cell RNA-seq datat from Puram et al. and create an expression matrix and metadata data.frame to create a cds object
puram_all_HSNCC_data <- read.table("GSE103322_HNSCC_all_data.txt", 
                                   header = TRUE, sep = "\t")

row.names(puram_all_HSNCC_data) <- puram_all_HSNCC_data[,1]
puram_all_HSNCC_data <- puram_all_HSNCC_data[,-1]

# Generate exprs data
puram_all_HSNCC_expression_data <- as.matrix(puram_all_HSNCC_data[6:dim(puram_all_HSNCC_data)[1],])
puram_all_HSNCC_expression_data <- apply(puram_all_HSNCC_expression_data,c(1,2),as.numeric)
puram_all_HSNCC_TPM_expression_data <- puram_all_HSNCC_expression_data

# Generate metadata
puram_all_HSNCC_metadata <- t(puram_all_HSNCC_data[1:5,])
puram_all_HSNCC_metadata <- as.data.frame(apply(puram_all_HSNCC_metadata,c(1,2),as.character))
patient_id <-  gsub( "_.*$", "", row.names(puram_all_HSNCC_metadata))
puram_all_HSNCC_metadata$patient_id <- patient_id

# Load Mock and TGFB cds objects created in Figure1 code to swap gene_short_names for ensembl gene ids
cds.list <- readRDS("pseudospace_processed_trajectories_cds.list.rds")

fData_subset <- fData(cds.list[["Mock"]])[fData(cds.list[["Mock"]])$gene_short_name %in% row.names(puram_all_HSNCC_expression_data),]

puram_all_HSNCC_TPM_expression_data_ENSid <- puram_all_HSNCC_TPM_expression_data[fData_subset$gene_short_name,]

# Generate feature data
puram_all_HSNCC_ENSid_genedata <- data.frame(row.names = as.character(fData_subset$id),
                                       id = as.character(fData_subset$id),
                                       gene_short_name = as.character(row.names(puram_all_HSNCC_TPM_expression_data_ENSid)),
                                      test_shortname = as.character(fData_subset$gene_short_name))

puram_all_HSNCC_TPM_expression_data_ENSid <- puram_all_HSNCC_TPM_expression_data_ENSid[order(puram_all_HSNCC_ENSid_genedata$gene_short_name),]
puram_all_HSNCC_ENSid_genedata <- puram_all_HSNCC_ENSid_genedata[order(puram_all_HSNCC_ENSid_genedata$gene_short_name),]

# Re-name exprs data row.names to match fData
row.names(puram_all_HSNCC_TPM_expression_data_ENSid) <- row.names(puram_all_HSNCC_ENSid_genedata)

puram_all_HSNCC_ENSid_genedata <- puram_all_HSNCC_ENSid_genedata[,-3]

# Generate a cds object
pd <- new("AnnotatedDataFrame", data = puram_all_HSNCC_metadata)
fd <- new("AnnotatedDataFrame", data = puram_all_HSNCC_ENSid_genedata)

# First create a CellDataSet from the relative expression levels
cds <- newCellDataSet(as.matrix(puram_all_HSNCC_TPM_expression_data_ENSid),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.1,
                       expressionFamily = tobit(Lower = 0.1))

# Use relative to absolute to approximate counts to be able to compare to MCF10A data
rpc_matrix <- relative2abs(cds, method = "num_genes")

HSNCC_cds <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd,
                       featureData = fd,
                       lowerDetectionLimit = 0.5,
                       expressionFamily = negbinomial.size())

# Re-work metadata column names
colnames(pData(HSNCC_cds)) <- c("Maxima_enzyme","lymph_node","clasified_as_cancer_cell","clasified_as_non_cancer_cell",
                                "non_cancer_cell_type","patient_id","Size_Factor")

saveRDS(HSNCC_cds,"HSNCC_cds.rds")

