
library(gplots)
# To get the venn diagram intersects
venn_intersects <- function(x_list){
  #input: 
  #x_list: list to be intersected, the names of the list is pre-defined
  x <- venn(x_list, names=names(x_list))
  intersections <- attr(x, "intersections")
  intersections_length <- map(intersections, length)
  intersections_sum <- map_df(intersections, length) %>% t() %>% as.data.frame()
  intersections_sum <- rownames_to_column(intersections_sum)
  colnames(intersections_sum) <- c("Intersection" ,"number_of_DE_genes")
  
  intersections_name <- names(intersections) %>% as.list()
  ## convert the list into a character vector
  Gene_ID <- character()
  for (i in seq_along(intersections)){
    for(j in seq_along(intersections[[i]])){
      Gene_ID <- c(Gene_ID, intersections[[i]][j])
    }
  }
  
  intersect_label <- map2(intersections_name, intersections_length, rep)
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