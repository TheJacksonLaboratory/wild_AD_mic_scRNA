#sce4: convert gene names to gene ID in the row name of cd11b.integrated object

# This script is to determine DE genes affected by Strain, Genotype, and Strain:Genotype interaction.
# Use EdgeR QFT method adapted from: Tutorial https://osca.bioconductor.org/multi-sample-comparisons.html#motivation-8
# Bench mark DE gene analysis identified EdgeR QFT as one of the top methods https://www.nature.com/articles/nmeth.4612 (Figure 5)


# load Seurat project 
library(tidyverse)
library(scater)
library(edgeR)

load("01_DE_gene/summed_H_cluster.rda") 
load("01_DE_gene/orig_symbol.rda")

unique(summed$label)


## generate DE gene list compared to one backgroup. 
## Set the backgroup strain as the first strain in strain_levels 
strain_levels_list <- list(c("B6", "CAST", "PWK", "WSB") , 
                           c("CAST", "PWK", "WSB", "B6"), 
                           c("PWK", "WSB", "B6", "CAST"), 
                           c("WSB", "B6", "CAST", "PWK"))

DE_list <- function(strain_levels, summed=summed){
  file_path_filtered <- paste("./01_DE_gene/DE_gene_list_vs", strain_levels[1], "/all_DE_gene_cluster_", sep = "")
  file_path_all <- paste("./01_DE_gene/DE_gene_list_vs", strain_levels[1], "/filtered_DE_gene_cluster_", sep = "")
  
  labels=unique(summed$label)
  coef_num <- c(2:5, 9:11)
  
  for(label in labels){
    current <- summed[,label==summed$label]
    
    y <- DGEList(counts(current), samples=colData(current))
    
    y$genes$Orig_Symbol <- orig.symbol
    
    y$genes$Symbol <- mapIds(org.Mm.eg.db, rownames(y),
                             keytype="ENSEMBL", column="SYMBOL")  # keytype="ENSEMBL", attach gene Symbol from database to DGElist y 
    y$genes$Entrezid <- mapIds(org.Mm.eg.db, rownames(y),
                               keytype="ENSEMBL", column="ENTREZID")
    y$genes$Genename <- mapIds(org.Mm.eg.db, rownames(y),
                               keytype="ENSEMBL", column="GENENAME")
    y$genes <- as.data.frame(y$genes)
    
    discarded <- isOutlier(y$samples$lib.size, log=TRUE, type="lower")
    y <- y[,!discarded]
    summary(discarded)
    
    keep <- filterByExpr(y, group=current$Group) 
    y <- y[keep,]
    
    y <- calcNormFactors(y)
    
    y$samples$Strain=gsub("B6J", "B6", y$samples$Strain)
    y$samples$Strain= factor(y$samples$Strain, levels=strain_levels)
    y$samples$Genotype= factor(y$samples$Genotype, levels = c("WT", "APP/PS1"))
    y$samples$batch= factor(y$samples$batch)
    design <- model.matrix(~Strain+Genotype+Strain:Genotype+batch, y$samples) # 
    
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    
    my_glmQLF <-function(coef_num, fit){
      ### glmQLFTest + output
      res <- glmQLFTest(glmfit=fit, coef = coef_num)
      res <- topTags(res, n=NULL, p.value=1)
      res_sigTags <- res[[1]]
      res_sigTags$coef <- res[[3]] %>% str_replace_all("\\:|\\/", "-")
      res_sigTags$cluster <- label
      res_sigTags <- rownames_to_column(res_sigTags, var = "ENSMUSG")
      return(res_sigTags)
    }
   
    res_table <- coef_num %>% map(my_glmQLF, fit=fit)
    coef_name <- colnames(design)[coef_num] %>% str_replace_all("\\:|\\/", "-")
    
    for(i in seq_along(res_table)){
      write_delim(res_table[[i]], paste(file_path_filtered, label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
      write_delim(res_table[[i]] %>% filter(FDR<0.05), paste(file_path_all, label, "_", coef_name[[i]], ".txt", sep=""), delim = "\t")
    }
    
  }
  
}

walk(strain_levels_list, DE_list, summed=summed)

