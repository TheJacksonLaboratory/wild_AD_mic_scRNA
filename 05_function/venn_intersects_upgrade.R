## generate intersection table compatible with Vennerable library

library(Vennerable)
library(tidyverse)

venn_intersects_upgrade <- function(x_list){
  tmp <- Venn(x_list)
  intersect_name <- tmp@IndicatorWeight %>% rownames()
  Weight <- tmp@IndicatorWeight %>% as.tibble()
  names(Weight) <- str_remove(names(Weight), "\\.")
  Weight$intersect_name <- intersect_name
  Set_name_weight <- Weight %>% filter(Weight!=0) %>% select(Weight, intersect_name) # extract the the intersection terms with "non-zero" elements
  Sets <- tmp@IntersectionSets
  Gene_ID <- character()
  for (i in seq_along(Sets)){
    for(j in seq_along(Sets[[i]])){
      Gene_ID <- c(Gene_ID, Sets[[i]][j])
    }
  }
  intersect_label <- map2(Set_name_weight$intersect_name, Set_name_weight$Weight, rep)
  intersect_label_vecter <- character()
  for (i in seq_along(intersect_label)){
    for(j in seq_along(intersect_label[[i]])){
      intersect_label_vecter <- c(intersect_label_vecter, intersect_label[[i]][j])
    }
  }
  
  gene_intersect <- data.frame(Gene_ID, intersect_label_vecter) 
  colnames(gene_intersect) <- c("Orig_Symbol", "Intersections")
  return(gene_intersect)
}