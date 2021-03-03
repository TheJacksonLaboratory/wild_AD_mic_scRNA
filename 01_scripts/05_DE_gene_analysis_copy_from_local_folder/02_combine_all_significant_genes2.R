library(tidyverse)

# combine genes from all significant genes of all coef terms into one table: 
## combine next for strain and cluster seperated...

strains <- c("B6", "CAST", "PWK", "WSB")

write_all_DEvsBgStrain <- function(strain){
  prefix <- paste("./01_DE_gene/DE_gene_list_vs", strain, sep = "")
  output_file <- paste("./02_combine/DE_vs", strain, "/DE_vs", strain, "_coef_cluster.txt", sep = "") 
  
  filtered_terms <- list.files(prefix, pattern = "filtered")
  filtered_files <- file.path(prefix, filtered_terms)
  
  df_filtered <- filtered_files%>% 
    map_df(read_delim, delim="\t", 
           col_types = list(
             col_character(), #ENSMUSG
             col_character(), #Orig_Symbol
             col_character(), #Symbol = 
             col_character(), #Entrezid = 
             col_character(), #Genename = 
             col_double(),# logFC = 
             col_double(),#logCPM = 
             col_double(),#F = 
             col_double(), #PValue = 
             col_double(), #FDR = 
             col_character(),#coef = 
             col_character()#cluster = 
           ))
  
  write_delim(df_filtered, output_file, delim = "\t")
}

walk(strains, write_all_DEvsBgStrain)

# combine DE_vs[BgStrain]_coef_cluster.txt into one

# combine genes from the same coef category for all strains

### grab all the q and FC from all comparisons
### select only the files with strain effects for each cluster ## vs B6

strain_levels_list <- list(c("B6", "CAST", "PWK", "WSB") , 
                           c("CAST", "PWK", "WSB", "B6"), 
                           c("PWK", "WSB", "B6", "CAST"), 
                           c("WSB", "B6", "CAST", "PWK"))

prefix <- paste("./01_DE_gene/DE_gene_list_vs", strain_levels_list[[1]][1], sep = "") # the loop item is strain_levels_list[[i]]
out

files <- list.files(prefix, pattern = "all")
df_files <- data_frame(file_name=files)
df_files <- df_files %>% separate(col= file_name, into=c("name"), sep = ".txt", remove=FALSE)
df_files <- df_files %>% separate(col= name, into=c("all","DE", "Gene", "cluster", "cluster_n","Coef"), sep = "_", remove=FALSE)
df_files <- df_files %>% select(file_name, cluster_n, Coef)
df_files <- df_files %>% filter(str_detect(Coef, pattern = "Strain"), str_detect(Coef, pattern = "-", negate = TRUE))

cluster <-  c("H", "6", "7", "8", "9", "10", "11", "12")
df_files <- df_files %>% mutate(cluster_n=factor(cluster_n, levels = cluster))

df_files_cluster <-df_files %>% group_split(cluster_n)

names(df_files_cluster) <- cluster

### combine logFC and FDR for each cluster: 
### take cluster H as an example: 

files <- file.path(prefix, df_files_cluster[[1]]$file_name) # i=1 
all_wilds_unfilter <- files %>% map(~read_delim(.,delim="\t"))
names(all_wilds_unfilter) <-strain_name <- strain_levels_list[[1]][-1] # # the loop item is strain_levels_list[[i]]

for (j in seq_along(all_wilds_unfilter)){
  colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], "vs", strain_levels_list[[1]][1], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
}

df_FCq <- all_wilds_unfilter[[strain_levels_list[[1]][2]]] %>% dplyr::select(1:5, 6,10) %>% 
  left_join(all_wilds_unfilter[[strain_levels_list[[1]][3]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
  left_join(all_wilds_unfilter[[strain_levels_list[[1]][4]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")

output_file <- paste("./02_combine/DE_vsB6/strain_vs", strain_levels_list[[1]][1], "_",cluster[[1]],".txt",sep = "") # cluster i=1
write_delim(df_FCq, output_file, delim = "\t") 


### now loop cluster, generate strain effect DE gene files for each cluster for one specific bg strain backgroud (vs B6)
#### also store the df_FCq object 
for(i in seq_along(cluster)){
  files <- file.path(prefix, df_files_cluster[[i]]$file_name) # i=1 
  all_wilds_unfilter <- files %>% map(~read_delim(.,delim="\t"))
  names(all_wilds_unfilter) <-strain_name <- strain_levels_list[[1]][-1] # # the loop item is strain_levels_list[[i]]
  
  for (j in seq_along(all_wilds_unfilter)){
    colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], "vs", strain_levels_list[[1]][1], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
  }
  
  df_FCq <- all_wilds_unfilter[[strain_levels_list[[1]][2]]] %>% dplyr::select(1:5, 6,10) %>% 
    left_join(all_wilds_unfilter[[strain_levels_list[[1]][3]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
    left_join(all_wilds_unfilter[[strain_levels_list[[1]][4]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")
  
  output_file <- paste("./02_combine/DE_vs",strain_levels_list[[1]][1], "/strain_vs", strain_levels_list[[1]][1], "_",cluster[[i]],".txt",sep = "") # cluster i=1
  write_delim(df_FCq, output_file, delim = "\t") 
}


### Now loop for all the strains

DE_vsBgStrain <- function(strain_levels){
  prefix <- paste("./01_DE_gene/DE_gene_list_vs", strain_levels[1], sep = "") # the loop item is strain_levels_list[[i]]
  files <- list.files(prefix, pattern = "all")
  df_files <- data_frame(file_name=files)
  df_files <- df_files %>% separate(col= file_name, into=c("name"), sep = ".txt", remove=FALSE)
  df_files <- df_files %>% separate(col= name, into=c("all","DE", "Gene", "cluster", "cluster_n","Coef"), sep = "_", remove=FALSE)
  df_files <- df_files %>% select(file_name, cluster_n, Coef)
  df_files <- df_files %>% filter(str_detect(Coef, pattern = "Strain"), str_detect(Coef, pattern = "-", negate = TRUE))
  
  cluster <-  c("H", "6", "7", "8", "9", "10", "11", "12")
  df_files <- df_files %>% mutate(cluster_n=factor(cluster_n, levels = cluster))
  
  df_files_cluster <-df_files %>% group_split(cluster_n)
  
  names(df_files_cluster) <- cluster
  
  
  for(i in seq_along(cluster)){
    files <- file.path(prefix, df_files_cluster[[i]]$file_name) # i=1 
    all_wilds_unfilter <- files %>% map(~read_delim(.,delim="\t"))
    names(all_wilds_unfilter) <-strain_name <- strain_levels[-1] # # the loop item is strain_levels_list[[i]]
    
    for (j in seq_along(all_wilds_unfilter)){
      colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], "vs", strain_levels[1], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
    }
    
    df_FCq <- all_wilds_unfilter[[strain_levels[2]]] %>% dplyr::select(1:5, 6,10) %>% 
      left_join(all_wilds_unfilter[[strain_levels[3]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
      left_join(all_wilds_unfilter[[strain_levels[4]]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")
    
    output_file <- paste("./02_combine/DE_vs",strain_levels[1], "/strain_vs", strain_levels[1], "_",cluster[[i]],".txt",sep = "") # cluster i=1
    write_delim(df_FCq, output_file, delim = "\t") 
  }
  
}

strain_levels_list <- list(c("B6", "CAST", "PWK", "WSB") , 
                           c("CAST", "PWK", "WSB", "B6"), 
                           c("PWK", "WSB", "B6", "CAST"), 
                           c("WSB", "B6", "CAST", "PWK"))

walk(strain_levels_list[[2]], DE_vsBgStrain)
            