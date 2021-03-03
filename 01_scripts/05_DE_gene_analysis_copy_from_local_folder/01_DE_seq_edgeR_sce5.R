#sce4: convert gene names to gene ID in the row name of cd11b.integrated object

# This script is to determine DE genes affected by Strain, Genotype, and Strain:Genotype interaction.
# Use EdgeR QFT method adapted from: Tutorial https://osca.bioconductor.org/multi-sample-comparisons.html#motivation-8
# Bench mark DE gene analysis identified EdgeR QFT as one of the top methods https://www.nature.com/articles/nmeth.4612 (Figure 5)


# load Seurat project 
library(tidyverse)
library(Seurat)
library(scater)
library(edgeR)


cd11b.integrated <- readRDS("../emase_2/output/Integration_mg/mg_int.rds")  # i=17 (PCA dim) # j=0.6 (resolution)

DefaultAssay(cd11b.integrated) <- "RNA"
cd11b.integrated@meta.data %>% group_by(Strain, Genotype) %>% summarise(N=n())

# bin homeostatic clusters 0,1,2,3,4,5 into H - homeostatic cluster
cd11b.integrated$clusters2 <- ifelse(cd11b.integrated$seurat_clusters %in% 0:5, "H", cd11b.integrated$seurat_clusters %>% as.character())
cd11b.integrated$clusters2 <- factor(cd11b.integrated$clusters2, levels = c("H", "6", "7", "8", "9", "10", "11", "12"))

# convert to single cell experiment
cd11b.integrated <- as.SingleCellExperiment(cd11b.integrated)

EnsembleID <- read_tsv("../documents/Ensembl2Symbol.tsv")

x <- rownames(cd11b.integrated)==EnsembleID$Symbol
### find the difference
rownames(cd11b.integrated)[!x]
EnsembleID$Symbol[!x]

all(gsub(".[0-9]$", "", rownames(cd11b.integrated)[!x]) == EnsembleID$Symbol[!x])

### replace the gene name to gene ensembl ID
orig.symbol <- rownames(cd11b.integrated)
rownames(cd11b.integrated) <- EnsembleID$`#Ensemble_ID` ## this command works here but doesn't work on cloud. 

# explore the Single Cell Experiment object
#cd11b.integrated

head(rownames(cd11b.integrated))

head(colnames(cd11b.integrated))

colData(cd11b.integrated) %>% head()

table(cd11b.integrated$seurat_clusters, cd11b.integrated$Group)

table(cd11b.integrated$seurat_clusters, cd11b.integrated$Strain)

# gridExtra::grid.arrange(
#   plotUMAP(cd11b.integrated, colour_by="Strain", text_by="seurat_clusters"),
#   plotUMAP(cd11b.integrated, colour_by="Genotype"),
#   ncol=2
# )


# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(cd11b.integrated, 
                               id=DataFrame(
                                 label=cd11b.integrated$clusters2,
                                 sample=cd11b.integrated$rep)
)

summed
save(summed, file = "01_DE_gene/summed_H_cluster.rda") 
save(orig.symbol, file = "01_DE_gene/orig_symbol.rda")

### Below are the testing to make DE gene list comparing all other strains to B6 or all other strains to WSB. Use "01_DE_seq_edgeR_sce6.R" to run DE genes across all strain background. The below code are the base single strain code, for record.

# unique(summed$label)
# label="H"
# current <- summed[,label==summed$label]
# 
# # Creating up a DGEList object for use in edgeR:
# 
# y <- DGEList(counts(current), samples=colData(current))
# y
# 
# ### add gene names
# library(org.Mm.eg.db)
# 
# y$genes$Orig_Symbol <- orig.symbol
# 
# y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
#                          keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist
# y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
#                            keytype="ENSEMBL", column="ENTREZID")
# y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
#                            keytype="ENSEMBL", column="GENENAME")
# y$genes <- as.data.frame(y$genes)
# 
# 
# head(y$genes$Genename)
# 
# # remove samples with low library size
# discarded <- isOutlier(y$samples$lib.size, log=TRUE, type="lower")
# y <- y[,!discarded]
# summary(discarded)
# 
# # remove genes that are lowly expressed
# keep <- filterByExpr(y, group=current$Group) # check group argument when filtering
# y <- y[keep,]
# #save(y, file = "01_DE_gene/y_H_cluster.rda") 
# summary(keep)
# 
# # Trimmed means of M-values of methods for normalization
# y <- calcNormFactors(y)
# y$samples
# 
# # statistical modeling
# 
# y$samples$Strain= factor(y$samples$Strain, levels=c("B6J", "CAST","PWK", "WSB"))
# y$samples$Strain2= factor(y$samples$Strain, levels=c("WSB", "CAST","PWK","B6J"))
# y$samples$Genotype= factor(y$samples$Genotype, levels = c("WT", "APP/PS1"))
# y$samples$batch= factor(y$samples$batch)
# 
# 
# str(y$samples)
# 
# design <- model.matrix(~Strain+Genotype+Strain:Genotype+batch, y$samples)
# design
# 
# y <- estimateDisp(y, design)
# summary(y$trended.dispersion)
# 
# plotBCV(y)
# 
# fit <- glmQLFit(y, design, robust=TRUE)
# summary(fit$var.prior)
# 
# summary(fit$df.prior)
# 
# plotQLDisp(fit)
# 
# colnames(coef(fit))
# res <- glmQLFTest(fit, coef=4)
# 
# 
# summary(decideTests(res))
# 
# x <- topTags(res, n=NULL, p.value=1)
# 
# 
# ## change model for everything comparing WSB
# design <- model.matrix(~Strain2+Genotype+Strain2:Genotype+batch, y$samples)
# design
# 
# y <- estimateDisp(y, design)
# summary(y$trended.dispersion)
# 
# plotBCV(y)
# 
# fit <- glmQLFit(y, design, robust=TRUE)
# summary(fit$var.prior)
# 
# summary(fit$df.prior)
# 
# plotQLDisp(fit)
# 
# colnames(coef(fit))
# res <- glmQLFTest(fit, coef=2)
# 
# 
# summary(decideTests(res))
# 
# x <- topTags(res, n=NULL, p.value=1)
# 
# 
# ## for each cluster perform test for all relevant coefficient, loop it for all clusters
# colnames(design)
# # select coeff 2 (StrainCAST), 3 (StrainPWK), 4 (StrainWSB), 5 (GenotypeAPP/PS1), 9 (StrainCAST:GenotypeAPP/PS1),10 (StrainPWK:GenotypeAPP/PS1),11 (StrainWSB:GenotypeAPP/PS1)
# 
# summed$batch= factor(summed$batch)
# 
# 
# labels=unique(summed$label)
# 
# coef_num <- c(2:5, 9:11)
# coef_name <- colnames(design)[coef_num] %>% str_replace_all("\\:|\\/", "-")
# 
# ## generate DE gene list compared to B6
# for(label in labels){
#   current <- summed[,label==summed$label]
#   
#   y <- DGEList(counts(current), samples=colData(current))
#   
#   y$genes$Orig_Symbol <- orig.symbol
#   
#   y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
#                            keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
#   y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
#                              keytype="ENSEMBL", column="ENTREZID")
#   y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
#                              keytype="ENSEMBL", column="GENENAME")
#   y$genes <- as.data.frame(y$genes)
#   
#   discarded <- isOutlier(y$samples$lib.size, log=TRUE, type="lower")
#   y <- y[,!discarded]
#   summary(discarded)
#   
#   keep <- filterByExpr(y, group=current$Group) 
#   y <- y[keep,]
#   
#   y <- calcNormFactors(y)
#   
#   y$samples$Strain= factor(y$samples$Strain, levels=c("B6J", "CAST","PWK", "WSB"))
#   y$samples$Genotype= factor(y$samples$Genotype, levels = c("WT", "APP/PS1"))
#   y$samples$batch= factor(y$samples$batch)
#   design <- model.matrix(~Strain+Genotype+Strain:Genotype+batch, y$samples) # 
#   
#   y <- estimateDisp(y, design)
#   fit <- glmQLFit(y, design, robust=TRUE)
#   
#   my_glmQLF <-function(coef_num, fit){
#     ### glmQLFTest + output
#     res <- glmQLFTest(glmfit=fit, coef = coef_num)
#     res <- topTags(res, n=NULL, p.value=1)
#     res_sigTags <- res[[1]]
#     res_sigTags$coef <- res[[3]] %>% str_replace_all("\\:|\\/", "-")
#     res_sigTags$cluster <- label
#     res_sigTags <- rownames_to_column(res_sigTags, var = "ENSMUSG")
#     return(res_sigTags)
#   }
#   
#   res_table <- coef_num %>% map(my_glmQLF, fit=fit)
#   
#   for(i in seq_along(res_table)){
#      write_delim(res_table[[i]], paste("./01_DE_gene/DE_gene_list_vsB6/all_DE_gene_cluster_", label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
#     write_delim(res_table[[i]] %>% filter(FDR<0.05), paste("./01_DE_gene/DE_gene_list_vsB6/filtered_DE_gene_cluster_", label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
#   }
# 
# }
# 
# 
# 
# ## generate DE gene list compared to WSB
# for(label in labels){
#   current <- summed[,label==summed$label]
#   
#   y <- DGEList(counts(current), samples=colData(current))
#   
#   y$genes$Orig_Symbol <- orig.symbol
#   
#   y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
#                            keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
#   y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
#                              keytype="ENSEMBL", column="ENTREZID")
#   y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
#                              keytype="ENSEMBL", column="GENENAME")
#   y$genes <- as.data.frame(y$genes)
#   
#   discarded <- isOutlier(y$samples$lib.size, log=TRUE, type="lower")
#   y <- y[,!discarded]
#   summary(discarded)
#   
#   keep <- filterByExpr(y, group=current$Group) 
#   y <- y[keep,]
#   
#   y <- calcNormFactors(y)
#   
#   y$samples$Strain2= factor(y$samples$Strain, levels=c("WSB", "B6J", "CAST","PWK"))
#   y$samples$Genotype= factor(y$samples$Genotype, levels = c("WT", "APP/PS1"))
#   y$samples$batch= factor(y$samples$batch)
#   design <- model.matrix(~Strain2+Genotype+Strain2:Genotype+batch, y$samples) # 
#   
#   y <- estimateDisp(y, design)
#   fit <- glmQLFit(y, design, robust=TRUE)
#   
#   my_glmQLF <-function(coef_num, fit){
#     ### glmQLFTest + output
#     res <- glmQLFTest(glmfit=fit, coef = coef_num)
#     res <- topTags(res, n=NULL, p.value=1)
#     res_sigTags <- res[[1]]
#     res_sigTags$coef <- res[[3]] %>% str_replace_all("\\:|\\/", "-")
#     res_sigTags$cluster <- label
#     res_sigTags <- rownames_to_column(res_sigTags, var = "ENSMUSG")
#     return(res_sigTags)
#   }
#   
#   res_table <- coef_num %>% map(my_glmQLF, fit=fit)
#   
#   for(i in seq_along(res_table)){
#     write_delim(res_table[[i]], paste("./01_DE_gene/DE_gene_list_vsWSB/all_DE_gene_cluster_", label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
#     write_delim(res_table[[i]] %>% filter(FDR<0.05), paste("./01_DE_gene/DE_gene_list_vsWSB/filtered_DE_gene_cluster_", label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
#   }
#   
# }
# 

