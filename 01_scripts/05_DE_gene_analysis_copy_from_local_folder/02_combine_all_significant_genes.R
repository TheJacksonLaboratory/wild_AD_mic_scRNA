library(tidyverse)

# combine genes from all significant genes of all coef terms into one table: 
prefix <- "./01_DE_gene/DE_gene_list_vsB6"
filtered_terms <- list.files(prefix, pattern = "filtered")
filtered_files <- file.path(prefix, filtered_terms)

df_filtered <- filtered_files %>% 
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

write_delim(df_filtered, "./02_combine/DE_vsB6/DE_vsB6_coef_cluster.txt", delim = "\t")

df_filtered %>% filter(coef == "StrainCAST-GenotypeAPP-PS1") %>% dim()
df_filtered %>% filter(coef == "StrainPWK-GenotypeAPP-PS1") %>% dim()
df_filtered %>% filter(coef == "StrainWSB-GenotypeAPP-PS1") %>% dim()

## make the same file for WSB
prefix <- "./01_DE_gene/DE_gene_list_vsWSB"
filtered_terms <- list.files(prefix, pattern = "filtered")
filtered_files <- file.path(prefix, filtered_terms)

df_filtered <- filtered_files %>% 
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

write_delim(df_filtered, "./02_combine/DE_vsWSB/DE_vsWSB_coef_cluster.txt", delim = "\t")

df_filtered %>% filter(coef == "Strain2CAST-GenotypeAPP-PS1") %>% dim()
df_filtered %>% filter(coef == "Strain2PWK-GenotypeAPP-PS1") %>% dim()
df_filtered %>% filter(coef == "Strain2B6J-GenotypeAPP-PS1") %>% dim()

# combine genes from the same coef category
### grab all the q and FC from all comparisons
### select only the files with strain effects for each cluster ## vs B6
prefix <- "./01_DE_gene/DE_gene_list_vsB6"
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
names(all_wilds_unfilter) <-strain_name <- c("CAST", "PWK", "WSB")

for (j in seq_along(all_wilds_unfilter)){
  colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
}

df_FCq <- all_wilds_unfilter[["CAST"]] %>% dplyr::select(1:5, 6,10) %>% 
  left_join(all_wilds_unfilter[["PWK"]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
  left_join(all_wilds_unfilter[["WSB"]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")

write_delim(df_FCq, paste("./02_combine/DE_vsB6/strain_vsB6_",cluster[[1]],".txt",sep = ""), delim = "\t") # i=1


### now loop cluster, generate strain effect DE gene files for each cluster
#### also store the df_FCq object 
for (i in seq_along(cluster)){
  files <- file.path(prefix, df_files_cluster[[i]]$file_name) # i=1
  all_wilds_unfilter <- files %>% map(~read_delim(.,delim="\t"))
  names(all_wilds_unfilter) <-strain_name <- c("CAST", "PWK", "WSB")
  
  for (j in seq_along(all_wilds_unfilter)){
    colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
  }
  
  df_FCq <- all_wilds_unfilter[["CAST"]] %>% dplyr::select(1:5, 6,10) %>% 
    left_join(all_wilds_unfilter[["PWK"]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
    left_join(all_wilds_unfilter[["WSB"]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")
  
  write_delim(df_FCq, paste("./02_combine/DE_vsB6/strain_vsB6_",cluster[[i]],".txt", sep = ""), delim = "\t")
}


### for vs WSB

# combine genes from the same coef category
### grab all the q and FC from all comparisons
### select only the files with strain effects for each cluster ## vs B6
prefix <- "./01_DE_gene/DE_gene_list_vsWSB"
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
all_wilds_unfilter %>% map(~.$coef %>% unique)
names(all_wilds_unfilter) <-strain_name <- c("PWK", "B6J", "CAST")

for (j in seq_along(all_wilds_unfilter)){
  colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
}

df_FCq <- all_wilds_unfilter[[1]] %>% dplyr::select(1:5, 6,10) %>% 
  left_join(all_wilds_unfilter[[2]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
  left_join(all_wilds_unfilter[[3]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")

write_delim(df_FCq, paste("./02_combine/DE_vsWSB/strain_vsWSB_",cluster[[1]],".txt",sep = ""), delim = "\t") # i=1


### now loop cluster, generate strain effect DE gene files for each cluster
#### also store the df_FCq object 
for (i in seq_along(cluster)){
  files <- file.path(prefix, df_files_cluster[[i]]$file_name) # i=1
  all_wilds_unfilter <- files %>% map(~read_delim(.,delim="\t"))
  names(all_wilds_unfilter) <-strain_name <- c("PWK", "B6J", "CAST")
  
  for (j in seq_along(all_wilds_unfilter)){
    colnames(all_wilds_unfilter[[j]])[6:12] <- paste(strain_name[j], colnames(all_wilds_unfilter[[j]])[6:12], sep="_")
  }
  
  df_FCq <- all_wilds_unfilter[[1]] %>% dplyr::select(1:5, 6,10) %>% 
    left_join(all_wilds_unfilter[[2]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG") %>% 
    left_join(all_wilds_unfilter[[3]] %>% dplyr::select(ENSMUSG, 6, 10), by="ENSMUSG")
  
  write_delim(df_FCq, paste("./02_combine/DE_vsWSB/strain_vsWSB_",cluster[[i]],".txt", sep = ""), delim = "\t")
}

            