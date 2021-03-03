# 

library(tidyverse)
strain_color <- c("#00AA00", "#FF0000", "#9900EE")
#names(strain_color) <- c("CAST", "PWK", "WSB")


## plot strain effect genes overlapping GWAS all in one plot
input_path <- "./09_venn_strain_vsB6_all_cluster_loop_newGWAS/"

output <- "10_GWAS_heatmap/"

cluster <- c("H", "6", "7", "8", "9", "10", "11", "12")

DE_GWAS_files <- paste(input_path, "cluster_", cluster, "_vennlist_vsB6_GWAS", ".txt", sep="")


df_list <- DE_GWAS_files %>% map(read_delim, delim="\t")
names(df_list) <- cluster

for (i in cluster){
  df_list[[i]]$Cluster <- i
}

df <- do.call(rbind, df_list)
df$Cluster <- factor(df$Cluster, levels = cluster)
gene_order <- table(df$Orig_Symbol) %>% sort() %>% names()
df$Orig_Symbol <- factor(df$Orig_Symbol, levels = gene_order %>% rev())

ggplot(df, aes(x=Cluster, y=Orig_Symbol)) +
  geom_point(aes(color=Intersections), size=3, alpha=0.7)+
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.text= element_text(family = "Arial", face = "bold"), 
        axis.text.y = element_text(face = "bold.italic"))+
  coord_fixed()
ggsave("10_GWAS_heatmap/gene_intersection_scatter.png", dpi = 300, width = 6, height = 4.85)


#save df object to plot gene expression in single cell. 
save(df, file = "10_GWAS_heatmap/gene_table_intersection_cluster_combined.rda")


## plot strain effect genes overlapping GWAS all in one plot seperated by Intersections
## all the genes need to show up in each plot, colored by strain effect
df <- select(df, Orig_Symbol:WSB_FDR, Cluster)

df <- separate(df, Intersections, into = c("Intersections"), sep = ":GWAS")

df <- df %>% mutate(Intersection_count = str_count(Intersections, "vs"))

df_1 <- df %>% filter(Intersection_count==1)

df_2 <- df %>% 
  filter(Intersection_count==2) %>% 
  separate(Intersections,into = c("Intersection1", "Intersection2"), sep = ":" ) %>% 
  gather(key = "key", value = "Intersections", Intersection1:Intersection2) %>% 
  select(1, Intersections, everything(), -key)
all(names(df_1)==names(df_2))

df_3 <- df %>% 
  filter(Intersection_count==3) %>% 
  separate(Intersections,into = c("Intersection1", "Intersection2", "Intersection3"), sep = ":" ) %>% 
  gather(key = "key", value = "Intersections", Intersection1:Intersection3) %>% 
  select(1, Intersections, everything(), -key)
all(names(df_1)==names(df_3))

df_tidy <- rbind(df_1, df_2, df_3)

df_tidy$Intersections <- factor(df_tidy$Intersections, levels = c("CAST vs B6","PWK vs B6","WSB vs B6" ))

ggplot(df_tidy, aes(x=Cluster, y=Orig_Symbol)) +
  geom_point(aes(color=Intersections), size=3, alpha=0.7)+
  theme_bw()+
  facet_wrap(.~Intersections)+
  scale_color_manual(values=strain_color)+
  theme(axis.title= element_blank(),
        axis.text= element_text(family = "Arial", face = "bold"), 
        axis.text.y = element_text(face = "bold.italic"),
        strip.text = element_text(family = "Arial", face = "bold"),
        legend.position = "none")+
  coord_fixed()
ggsave("10_GWAS_heatmap/gene_intersection_scatter_sep.png", dpi = 300, width = 5.9, height = 4.5)

# then plot heatmap

# plot one cluster
# strain-specific microglia cluster genes with GWAS 
gene_to_plot <- filter(df, Cluster==cluster[[2]]) %>% .$Orig_Symbol %>% unique()

df_FC <-df %>% 
  filter(Cluster==cluster[[2]]) %>% 
  select(Orig_Symbol, contains("FC")) %>% 
  mutate(Orig_Symbol = factor(Orig_Symbol, levels = gene_to_plot)) %>% 
  gather(key = "Strain", value = "log2FC", -Orig_Symbol) %>% 
  mutate(Strain=factor(Strain, 
                       levels = c("CAST_logFC", "PWK_logFC", "WSB_logFC"),
                       labels = c("CAST", "PWK", "WSB")))

df_FDR <- df %>% 
  filter(Cluster==cluster[[2]]) %>% 
  select(Orig_Symbol, contains("FDR")) %>% 
  mutate(Orig_Symbol = factor(Orig_Symbol, levels = gene_to_plot)) %>% 
  gather(key = "Strain", value = "FDR", -Orig_Symbol) %>% 
  mutate(Strain=factor(Strain, 
                       levels = c("CAST_FDR", "PWK_FDR", "WSB_FDR"),
                       labels = c("CAST", "PWK", "WSB")))

all(df_FDR$Orig_Symbol==df_FC$Orig_Symbol)

df_FCq <- df_FC
df_FCq$FDR <- df_FDR$FDR

ggplot(df_FCq, aes(x=Strain, y=Orig_Symbol, fill=log2FC)) + 
  geom_tile(color="grey50")+
  scale_fill_gradient2(name="log2(FC)", low = scales::muted("blue"), mid = "white",
                       high = scales::muted("red"))+
  theme_bw()+
  #scale_x_discrete(label=c("WT", "Ddit3-/-", "Jun-/-", "Jun-/- Ddit3-/-"))+
  theme(text=element_text(family = "Arial"),
        axis.text.y = element_text(size = 8, face = "bold.italic"),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=8, 
        #                           face=c("bold"),
        #                           angle = 45, 
        #                           vjust = 1.2,
        #                           hjust = 1),
        axis.title = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "right",
        legend.direction = "vertical" )+ 
  coord_fixed()+
  geom_point(data=df_FCq %>% filter(FDR>=0.05), size=0.3) # overlay the q value on top of the matrix
ggsave(paste(output, "hmap_GWAS_cluster_", cluster[[2]], ".png", sep = ""), 
       width = 2.5, height = 2.4/10*length(gene_to_plot), units = "in", dpi = 300)

##########
# plot all the clusters


for (i in seq_along(cluster)){
  # plot one cluster
  # strain-specific microglia cluster genes with GWAS 
  gene_to_plot <- filter(df, Cluster==cluster[[i]]) %>% .$Orig_Symbol %>% unique()
  
  df_FC <-df %>% 
    filter(Cluster==cluster[[i]]) %>% 
    select(Orig_Symbol, contains("FC")) %>% 
    mutate(Orig_Symbol = factor(Orig_Symbol, levels = gene_to_plot)) %>% 
    gather(key = "Strain", value = "log2FC", -Orig_Symbol) %>% 
    mutate(Strain=factor(Strain, 
                         levels = c("CAST_logFC", "PWK_logFC", "WSB_logFC"),
                         labels = c("CAST", "PWK", "WSB")))
  
  df_FDR <- df %>% 
    filter(Cluster==cluster[[i]]) %>% 
    select(Orig_Symbol, contains("FDR")) %>% 
    mutate(Orig_Symbol = factor(Orig_Symbol, levels = gene_to_plot)) %>% 
    gather(key = "Strain", value = "FDR", -Orig_Symbol) %>% 
    mutate(Strain=factor(Strain, 
                         levels = c("CAST_FDR", "PWK_FDR", "WSB_FDR"),
                         labels = c("CAST", "PWK", "WSB")))
  
  all(df_FDR$Orig_Symbol==df_FC$Orig_Symbol)
  
  df_FCq <- df_FC
  df_FCq$FDR <- df_FDR$FDR
  
  #####
  
  ggplot(df_FCq, aes(x=Strain, y=Orig_Symbol, fill=log2FC)) + 
    geom_tile(color="grey50")+
    scale_fill_gradient2(name="log2(FC)", low = scales::muted("blue"), mid = "white",
                         high = scales::muted("red"))+
    theme_bw()+
    #scale_x_discrete(label=c("WT", "Ddit3-/-", "Jun-/-", "Jun-/- Ddit3-/-"))+
    theme(text=element_text(family = "Arial"),
          axis.text.y = element_text(size = 8, face = "bold.italic"),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(size=8, 
          #                           face=c("bold"),
          #                           angle = 45, 
          #                           vjust = 1.2,
          #                           hjust = 1),
          axis.title = element_blank(),
          axis.ticks = element_blank(), 
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "right",
          legend.direction = "vertical" )+ 
    coord_fixed()+
    geom_point(data=df_FCq %>% filter(FDR>=0.05), size=0.3) # overlay the q value on top of the matrix
  ggsave(paste(output, "hmap_GWAS_cluster_", cluster[[i]], ".png", sep = ""), 
         width = 2.5, height = 2.4/10*length(gene_to_plot), units = "in", dpi = 300)
}

