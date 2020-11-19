## get DE genes based on FDR and FC. 
## The input file is DE gene analysis from Seurat package


cluster_DE <- function(logFDR_cut, logFC_cut, df){
  # input:
  # logFDR_cut, threthold of -log10FDR
  # logFC_cut, threthold of log2FC
  # df: DE gene table for each cluster
  # output: 
  # a list containing marker genes for a given cluster for all strains 
  # strain: global environment
  DE_list <- vector(mode = "list", length = length(strain))
  
  if(logFC_cut>0){
    for(i in seq_along(strain)){
      DE_list[[i]] <- df %>% select(Symbol, contains(strain[i])) %>% 
        filter(str_detect(Symbol, "^Gm", negate = TRUE)) %>% 
        filter_at(vars(contains("p_val_adj")), any_vars(-log10(.)>logFDR_cut)) %>% 
        filter_at(vars(contains("logFC")), any_vars(.>logFC_cut)) %>% 
        select(Symbol) %>% unlist()
    }
  }else{
    for(i in seq_along(strain)){
      DE_list[[i]] <- df %>% select(Symbol, contains(strain[i])) %>% 
        filter(str_detect(Symbol, "^Gm", negate = TRUE)) %>% 
        filter_at(vars(contains("p_val_adj")), any_vars(-log10(.)>logFDR_cut)) %>% 
        filter_at(vars(contains("logFC")), any_vars(.<logFC_cut)) %>% 
        select(Symbol) %>% unlist()
    }
  }
  
  names(DE_list) <- strain
  return(DE_list)
}
