suppressPackageStartupMessages({
    library(devtools)
    library(RColorBrewer)
    library(plyr)
    library(dplyr)
    library(tidyr)
    library(piano)
    library(ggplot2)
    library(glmnet)
    library(parallel)
    library(pheatmap)
    library(reshape2)
    library(scales)
    library(viridis)
    library(stringr)
    library(gridExtra)
    library(ggrastr)
    #### Note that monocle 3 alpha (monocle_2.99.3) is used in this script as monocle2 
    #### does not have the functionality to cluster in PCA space
    library(monocle)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

##### Load and define necessary functions #####
source("Pseudospace_support_functions.R")

plot_marker_tSNE <- function(cds, marker_genes, ncol){
    
    plot.list <- list()
    
    for(marker in marker_genes){
     
        cds_subset <- cds[fData(cds)$gene_short_name == marker,]
        marker_expr <- Matrix::t(Biobase::exprs(cds_subset))/pData(cds_subset)$Size_Factor
        
        pData(cds_subset)$marker <- log10(as.numeric(marker_expr[,1]) + 0.1)
        
        color_str <- paste0("log10(",marker," + 0.1)")
        
        plot.list[[marker]] <- plot_cell_clusters(cds_subset, color_by = "marker", cell_size = 0.5) +
        scale_color_viridis(option = "magma", name = color_str) +
        ggtitle(marker) +
        theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size = 6),
              legend.key.width = unit(0.25, "cm"),
              legend.key.height = unit(0.25, "cm"))
        
    }
    
    plot.list[["NTC"]] <- plot_cell_clusters(cds.list[["TGFB"]], color_by = "NONTARGETING" , cell_size = 0.5) +
        scale_color_manual("Expressed sgRNA", labels = c("2_TRUE" = "Nontargeting", "1_FALSE" = "Other"), 
                          values = c("2_TRUE" = "#D00000", "1_FALSE" = "#B0B5B3")) +
        theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size = 10),
              legend.key.width = unit(0.25, "cm"),
              legend.key.height = unit(0.3, "cm")) +
    guides(color = guide_legend(override.aes = list(size=2)))

    
    plot.list[["TGFBR2"]] <- plot_cell_clusters(cds.list[["TGFB"]], color_by = "TGFBR2", cell_size = 0.5) +
        scale_color_manual("Expressed sgRNA", labels = c("2_TRUE" = "TGFBR2", "1_FALSE" = "Other"), 
                           values = c("2_TRUE" = "#D00000", "1_FALSE" = "#B0B5B3")) +
        theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size = 10),
              legend.key.width = unit(0.25, "cm"),
              legend.key.height = unit(0.3, "cm")) +
    guides(color = guide_legend(override.aes = list(size=2)))
    
    
    print(names(plot.list))
    pdf("Marker_array.pdf")
    do.call("grid.arrange", c(plot.list, ncol = ncol))
    dev.off()

}

# Expecation maximitation to infer guide efficiency
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

#### Load data ####
pseudospace_cds <- readRDS("CROPseq_pseudospace_cds.rds")

pseudospace_cds <- pseudospace_cds[,with(pData(pseudospace_cds), !is.na(proportion) & 
                                         guide_count == 1 &
                                        proportion > 0.8)]

pseudospace_cds <- updateCDS(pseudospace_cds)

# Create a cds subset for each stimulation condition that contains spatially isolated CROPseq cells
cds.list <- list()

cds.list[["Mock"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment == "mock"]

cds.list[["TGFB"]] <- pseudospace_cds[,pData(pseudospace_cds)$treatment == "tgfb"]

# Identify genes that are expressed in at least 50 of cells
expressed_genes.list <- list()

for(sample in names(cds.list)){
    
    expressed_genes.list[[sample]] <- row.names(fData(cds.list[[sample]])[Matrix::rowSums(Biobase::exprs(cds.list[[sample]]) > 0) > 50 ,])

}

for(sample in names(cds.list)) {
    
    cds.list[[sample]] <- detectGenes(cds.list[[sample]], min_expr = 0.5)
    cds.list[[sample]] <- estimateSizeFactors(cds.list[[sample]])
    
}

for(sample in names(cds.list)){
    
    cds.list[[sample]] <- preprocessCDS(cds.list[[sample]][expressed_genes.list[[sample]]],
                                        method = "PCA",
                                        num_dim = 25,
                                        norm_method = "log", 
                                        verbose = T,
                                        cores = 1)
    
}

for(sample in names(cds.list)){
    
    cds.list[[sample]] <- reduceDimension(cds.list[[sample]][expressed_genes.list[[sample]]], 
                                          reduction_method = "tSNE",
                                          max_components = 2, 
                                          norm_method = "log",
                                          num_dim = 25,
                                          verbose = T,
                                          cores = 1)
    
}

for(sample in names(cds.list)){
    
    cds.list[[sample]] <- clusterCells(cds.list[[sample]],
                                       use_pca = TRUE,
                                       method = 'louvain',
                                       cores = 1)
    
}

plot_cell_clusters(cds.list[["TGFB"]], cell_size = 0.5) +
theme(legend.position = "right", 
      plot.title = element_text(hjust = 0.5),
      legend.title=element_text(size=10),
      legend.text=element_text(size = 10),
      legend.key.width = unit(0.25, "cm"),
      legend.key.height = unit(0.3, "cm")) +
scale_color_manual("Cluster", 
                   values = c("1"="firebrick1", 
                              "2"="lemonchiffon4", 
                              "3"="deepskyblue2", 
                              "4"="slategray4", 
                              "5"="navy", 
                              "6"="brown4", 
                              "7"="darkgreen", 
                              "8"="gold2", 
                              "9" = "orangered3", 
                              "10" = "darkolivegreen",
                              "11" = "chartreuse4",
                              "12" = "plum4", 
                              "13" = "darkorchid2",
                              "14" = "darkcyan",
                             "15" = "aquamarine3",
                             "16" = "darkgoldenrod4",
                             "17" = "dimgrey",
                             "18" = "firebrick3")) +
guides(color = guide_legend(override.aes = list(size=2))) +
ggsave("Crop-seq_TGFB_driven_EMT_tSNE_by_PCA_louvain_cluster.png",
      width = 4, height = 3.5)

plot_genes_violin(cds.list[["TGFB"]][fData(cds.list[["TGFB"]])$gene_short_name == 
                                     "VIM",], 
                  grouping = "Cluster",
                  fill_by = "Cluster",
                  plot_trend = TRUE,
                  log_scale = TRUE) +
scale_fill_manual("Cluster", 
                   values = c("1"="firebrick1", 
                              "2"="lemonchiffon4", 
                              "3"="deepskyblue2", 
                              "4"="slategray4", 
                              "5"="navy", 
                              "6"="brown4", 
                              "7"="darkgreen", 
                              "8"="gold2", 
                              "9" = "orangered3", 
                              "10" = "darkolivegreen",
                              "11" = "chartreuse4",
                              "12" = "plum4", 
                              "13" = "darkorchid2",
                              "14" = "darkcyan",
                             "15" = "aquamarine3",
                             "16" = "darkgoldenrod4",
                             "17" = "dimgrey",
                             "18" = "firebrick3")) +
theme(legend.position = "right", 
      plot.title = element_text(hjust = 0.5),
      legend.title=element_text(size=10),
      legend.text=element_text(size = 10),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.3, "cm")) +
ggsave("Crop-seq_TGFB_driven_EMT_VIM_levels_by_PCA_louvain_cluster.png",
      width = 6, height = 3.5)
    

plot_marker_tSNE(cds.list[["TGFB"]], 
                 c("CDH1","CRB3","FN1","VIM"), 
                 ncol = 2)

plot_cell_clusters(cds.list[["TGFB"]], color_by = "TGFBR1", cell_size = 0.5) +
scale_color_manual("Expressed sgRNA", labels = c("2_TRUE" = "TGFBR1", "1_FALSE" = "Other"), 
                   values = c("2_TRUE" = "#D00000", "1_FALSE" = "#B0B5B3")) +
 theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size = 10),
              legend.key.width = unit(0.25, "cm"),
              legend.key.height = unit(0.3, "cm")) +
guides(color = guide_legend(override.aes = list(size=2))) +
ggsave("Crop-seq_TGFB_driven_EMT_tSNE_by_TGFBR1.png",
      width = 4.8, height = 3.5)

plot_cell_clusters(cds.list[["TGFB"]], color_by = "TGFBR2", cell_size = 0.5) +
scale_color_manual("Expressed sgRNA", labels = c("2_TRUE" = "TGFBR2", "1_FALSE" = "Other"), 
                   values = c("2_TRUE" = "#D00000", "1_FALSE" = "#B0B5B3")) +
 theme(legend.position = "right", 
              plot.title = element_text(hjust = 0.5),
              legend.title=element_text(size=10),
              legend.text=element_text(size = 10),
              legend.key.width = unit(0.25, "cm"),
              legend.key.height = unit(0.3, "cm")) +
guides(color = guide_legend(override.aes = list(size=2))) +
ggsave("Crop-seq_TGFB_driven_EMT_tSNE_by_TGFBR2.png",
      width = 4.8, height = 3.5)

overall_TGFBR <- pData(cds.list[["TGFB"]]) %>% 
                 filter(gene %in% c("NONTARGETING","TGFBR1","TGFBR2")) %>% 
                 group_by(gene) %>% summarize(n = n()) 

VIM_low_TGFBR <- pData(cds.list[["TGFB"]])%>% 
                 filter(gene %in% c("NONTARGETING","TGFBR1","TGFBR2")) %>% 
                 filter(Cluster %in% c(1,3,6,7,14)) %>% 
                 group_by(gene) %>% 
                 summarize(n = n()) %>%  
                 arrange(n)

ggplot(overall_TGFBR, aes(x = gene, y = n)) + 
geom_bar(stat = "identity") +
ylab("Total cell number\n") +
xlab("Expressed sgRNA") +
theme_classic(base_size = 20) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_x_discrete(labels = c("NONTARGETING"="NTC", "TGFBR1" = "TGFBR1", "TGFBR2" = "TGFBR2")) +
ggsave("Pseudospace_CFG_TGFB_NTC-TGFBRs_cell_number.png", width = 4, height = 5)

VIM_low_freq <- data.frame(as.numeric(VIM_low_TGFBR$n/overall_TGFBR$n))
row.names(VIM_low_freq) <- c("NONTARGETING","TGFBR1","TGFBR2")
VIM_low_freq$gene <- row.names(VIM_low_freq)
colnames(VIM_low_freq) <- c("fraction", "gene")

ggplot(VIM_low_freq, aes(x = gene, y = fraction)) + 
geom_bar(stat = "identity") +
ylab("Fraction in\nVIM low clusers") +
xlab("Expressed sgRNA") +
theme_classic(base_size = 20) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_x_discrete(labels = c("NONTARGETING"="NTC", "TGFBR1" = "TGFBR1", "TGFBR2" = "TGFBR2")) +
ggsave("Pseudospace_CFG_TGFB_NTC-TGFBRs_VIM_low_frequency.png", width = 4, height = 5)

# Examine the expression levels of meenchymal markers across TGFB exposed cells expressing non-targeting control sgRNAs or sgRNA against TGFBRs

VIM_FN1_cds <- cds.list[["TGFB"]][fData(cds.list[["TGFB"]])$gene_short_name %in% c("FN1","VIM"),
                                  pData(cds.list[["TGFB"]])$gene %in% c("NONTARGETING","TGFBR1","TGFBR2")]

plot_genes_violin(VIM_FN1_cds, grouping = "gene",fill = "gene", min_expr = 0.1, 
                  plot_trend = TRUE, log_scale = TRUE, ncol = 2) +
theme_classic() + 
monocle:::monocle_theme_opts() +
scale_fill_manual(values = c("NONTARGETING" = "#9B9696","TGFBR1"="#458F96","TGFBR2" = "#7DAA92")) +
scale_color_manual(values = c("NONTARGETING" = "#000000","TGFBR1"="#000000","TGFBR2" = "#000000")) +
scale_x_discrete(labels = c("NTC","TGFBR1","TGFBR2"))+
xlab("Expressed sgRNA") +
theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none",text=element_text(size=20, family="Arial")) + 
ggsave("TGFB_treated_pseudospace_CFG_FN1_VIM_violin.png", width = 6, height = 5)

# Bin cluster with low vimenting together and test for enrichment of KOs across this group of clusters
new_cluster <- sapply(pData(cds.list[["TGFB"]])$Cluster, function(x){

    if(x %in% c(1,3,6,7,14))return("vim_low")
    return(x)

})

pData(cds.list[["TGFB"]])$Cluster <- factor(new_cluster, levels = unique(new_cluster))

#### Determine enrichment across tSNE clusters as implemented in Hill, A.J., McFaline-Figueroa, J.L. et al. (2018) Nat Methods doi:10.1038/nmeth.4604
analysis.guides = list()

for(sample in names(cds.list)){
pData(cds.list[[sample]]) %>% filter(guide_count == 1, gene != "NONTARGETING") %>% group_by(gene, barcode) %>%
    summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
    filter(n.guide.cells >= 10) %>% ungroup() %>%
    arrange(-n.target.cells, -n.guide.cells) %>% head(10)

analysis.guides[[sample]] = 
    (pData(cds.list[[sample]]) %>% filter(guide_count == 1, gene != "NONTARGETING") %>% group_by(gene, barcode) %>%
    summarize(n.guide.cells = n()) %>% group_by(gene) %>% mutate(n.target.cells = sum(n.guide.cells)) %>%
    filter(n.guide.cells >= 10) %>% ungroup())$barcode

length(analysis.guides[[sample]])
    }

analysis.targets = list()

analysis.targets[["Mock"]] = as.data.frame(pData(cds.list[["Mock"]]) %>%
    group_by(gene) %>% summarize(
        n.cells = n(),
        n.guides = length(intersect(unique(barcode), analysis.guides[["Mock"]]))) %>%
    filter(n.cells >= 15, n.guides >= 1) %>% select(gene))[, 1]

analysis.targets[["TGFB"]] = as.data.frame(pData(cds.list[["TGFB"]]) %>%
    group_by(gene) %>% summarize(
        n.cells = n(),
        n.guides = length(intersect(unique(barcode), analysis.guides[["TGFB"]]))) %>%
    filter(n.cells >= 15, n.guides >= 1) %>% select(gene))[, 1]

analysis.targets[["Mock"]]
analysis.targets[["TGFB"]]

target.to.guide.map = list()

target.to.guide.map[["Mock"]] = list()

target.to.guide.map[["TGFB"]] = list()

for (target in analysis.targets[["Mock"]]) {
    target.to.guide.map[["Mock"]][[target]] = 
        sort(unique(as.data.frame(pData(cds.list[["Mock"]]) %>%
            filter(gene == target, barcode %in% analysis.guides[["Mock"]]) %>%
            select(barcode))[, 1]))
}

for (target in analysis.targets[["TGFB"]]) {
    target.to.guide.map[["TGFB"]][[target]] = 
        sort(unique(as.data.frame(pData(cds.list[["TGFB"]]) %>%
            filter(gene == target, barcode %in% analysis.guides[["TGFB"]]) %>%
            select(barcode))[, 1]))
}

guide.to.target.map = list()

for(sample in names(cds.list)){
guide.to.target.map[[sample]] = list()

for (target in analysis.targets[[sample]]) {
    for (guide in target.to.guide.map[[sample]][[target]]) {
        guide.to.target.map[[sample]][[guide]] = target
    }
}

}

target.cluster.mat = list()

target.cluster.mat[["Mock"]] = acast(
    pData(cds.list[["Mock"]]) %>%
    filter(barcode %in% analysis.guides[["Mock"]] | gene == "NONTARGETING") %>%
    mutate(dummy = 1) %>% select(gene, Cluster, dummy),
    gene ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

target.cluster.mat[["TGFB"]] = acast(
    pData(cds.list[["TGFB"]]) %>%
    filter(barcode %in% analysis.guides[["TGFB"]] | gene == "NONTARGETING") %>%
    mutate(dummy = 1) %>% select(gene, Cluster, dummy),
    gene ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

NTC.cluster.p = list()

for(sample in names(cds.list)){
NTC.cluster.p[[sample]] <- pData(cds.list[[sample]])[pData(cds.list[[sample]])$gene == "NONTARGETING",] %>% 
    group_by(Cluster) %>% 
    summarize(n = n()) %>% 
    complete(Cluster, fill = list(n = 0.1))
}

guide.cluster.mat = list()

guide.cluster.mat[["Mock"]] = acast(
    pData(cds.list[["Mock"]]) %>% filter(barcode %in% analysis.guides[["Mock"]]) %>%
    mutate(dummy = 1) %>% select(barcode, Cluster, dummy),
    barcode ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

guide.cluster.mat[["TGFB"]] = acast(
    pData(cds.list[["TGFB"]]) %>% filter(barcode %in% analysis.guides[["TGFB"]]) %>%
    mutate(dummy = 1) %>% select(barcode, Cluster, dummy),
    barcode ~ Cluster, value.var = "dummy", fun.aggregate = sum, fill = 0)

ntc.distribution = list()

for(sample in names(cds.list)){
ntc.distribution[[sample]] = target.cluster.mat[[sample]]["NONTARGETING",]
ntc.distribution[[sample]] = ntc.distribution[[sample]] / sum(ntc.distribution[[sample]])
}

initial.target.level.chisq.pval = list()

for(sample in names(cds.list)){
set.seed(42)
initial.target.level.chisq.pval[[sample]] = sapply(
    analysis.targets[[sample]], function(target) {
    
    message(target)
    chisq.test(
        target.cluster.mat[[sample]][target,],
        p = NTC.cluster.p[[sample]]$n,
        simulate.p.value = F, rescale.p = T, B = 20000)$p.value  
})
}

initial.guide.level.chisq.pval <- list()

for(sample in names(cds.list)){
set.seed(42)
initial.guide.level.chisq.pval[[sample]] = sapply(
    analysis.guides[[sample]], function(guide) {
    
    message(guide)
    chisq.test(
        guide.cluster.mat[[sample]][guide,],
        p = NTC.cluster.p[[sample]]$n,
        simulate.p.value = F, rescale.p = T, B = 20000)$p.value  
})
}

initial.target.level.chisq.qval <- list()

for(sample in names(cds.list)){
    
    initial.target.level.chisq.qval[[sample]] <- p.adjust(initial.target.level.chisq.pval[[sample]], method = "BH")
  
}

for(sample in names(cds.list)){
    
    initial.target.level.chisq.qval[[sample]] <- sapply(initial.target.level.chisq.qval[[sample]], 
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})
}

initial.guide.level.chisq.qval <- list()

for(sample in names(cds.list)){
    
    initial.guide.level.chisq.qval[[sample]] <- p.adjust(initial.guide.level.chisq.pval[[sample]], method = "BH")

}

for(sample in names(cds.list)){
    
    initial.guide.level.chisq.qval[[sample]] <- sapply(initial.guide.level.chisq.qval[[sample]], 
                                                       function(x){if(x < 1e-50){return(1e-50)}else{return(x)}})
}

pass.target.level.screen = list()

for(sample in names(cds.list)){
    
    pass.target.level.screen[[sample]] = 
    sort(names(which(initial.target.level.chisq.qval[[sample]] < 0.05)))
    print(pass.target.level.screen[[sample]])
}

pass.guide.level.screen = list()

for(sample in names(cds.list)){
    
    pass.guide.level.screen[[sample]] = sort(unlist(unique(sapply(
        names(which(initial.guide.level.chisq.qval[[sample]] < 0.05)), function(guide) {
        guide.to.target.map[[sample]][[guide]]
    }))))
}

targets.passing.initial.screen = list()

for(sample in names(cds.list)){
targets.passing.initial.screen[[sample]] = sort(union(
    pass.target.level.screen[[sample]], pass.guide.level.screen[[sample]]))
}


targets.passing.initial.screen[["Mock"]]
targets.passing.initial.screen[["TGFB"]]


weighted.target.cluster.mat = list()

for (condition in c("Mock", "TGFB")) {
    weighted.target.cluster.mat[[condition]] = t(sapply(targets.passing.initial.screen[[condition]],
            function(target) {
        guides = target.to.guide.map[[condition]][[target]]
        if (length(guides) == 1) {
            return(target.cluster.mat[[condition]][target,])
        } else {
            mat = guide.cluster.mat[[condition]][guides,]
            guide.weights = get.guide.weights(mat, ntc.distribution[[condition]])
            guide.weights = guide.weights / max(guide.weights)
            print(condition)
            print(target)
        
            print(round(guide.weights, 3))
            return(round(colSums(sweep(mat, 1, guide.weights, "*"))))
        }
    }))
}

cluster.enrichment.df = list()

for (condition in c("Mock", "TGFB")) {
    weighted.mat = weighted.target.cluster.mat[[condition]]
    ntc.counts = target.cluster.mat[[condition]]["NONTARGETING",]

    cluster.enrichment.df[[condition]] = do.call(rbind, lapply(rownames(weighted.mat), function(target) {
        do.call(rbind, lapply(1:ncol(weighted.mat), function(cluster) {
            test = fisher.test(cbind(
                c(weighted.mat[target, cluster], sum(weighted.mat[target, -cluster])),
                c(ntc.counts[cluster], sum(ntc.counts[-cluster]))))

            data.frame(
                target = target,
                cluster = colnames(weighted.mat)[cluster],
                odds.ratio = unname(test$estimate),
                p.value = test$p.value)
        }))
    }))

    cluster.enrichment.df[[condition]]$q.value = p.adjust(cluster.enrichment.df[[condition]]$p.value, "fdr")

    cluster.enrichment.df[[condition]]$log2.odds = with(cluster.enrichment.df[[condition]],
        ifelse(odds.ratio == 0, -5, log2(odds.ratio)))
}

#### Print the enrichment of TGFBR KOs across vimentin low clusters 
cluster.enrichment.df[["TGFB"]] %>% 
filter(target %in% c("TGFBR1","TGFBR2")) %>% 
arrange(target, cluster, desc(q.value))

cluster.enrichment.df[["TGFB"]] %>% 
filter(cluster == "vim_low", q.value < 0.05) %>% 
arrange(target, cluster, desc(q.value))


