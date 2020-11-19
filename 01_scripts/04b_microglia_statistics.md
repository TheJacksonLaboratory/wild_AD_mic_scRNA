---
title: "Microglia  Statistics"
output: 
  html_document:
    keep_md: true
---




```r
library(tidyverse)
library(cowplot)
library(Seurat)

path = "../03_results/04b_microglia_statistics/"
```

### Load object

```r
cd11b.integrated <- readRDS("../02_data/intermediate_rds/mg_int.rds")  # i=17 (PCA dim) # j=0.6 (resolution)

DefaultAssay(cd11b.integrated) <- "RNA"
cd11b.integrated@meta.data %>% group_by(Strain, Genotype) %>% summarise(N=n())
```

```
## # A tibble: 8 x 3
## # Groups:   Strain [4]
##   Strain Genotype     N
##   <chr>  <chr>    <int>
## 1 B6J    APP/PS1  11531
## 2 B6J    WT        9201
## 3 CAST   APP/PS1   9232
## 4 CAST   WT       14892
## 5 PWK    APP/PS1   7906
## 6 PWK    WT       11796
## 7 WSB    APP/PS1  14425
## 8 WSB    WT        8763
```

```r
cd11b.integrated$Strain <- str_replace_all(cd11b.integrated$Strain, pattern = "B6J", replacement = "B6")
cd11b.integrated$Group <- str_replace_all(cd11b.integrated$Group, pattern = "B6J", replacement = "B6")
cd11b.integrated$Group <- factor(cd11b.integrated$Group, levels = c("B6_WT","B6_APP/PS1","CAST_WT", "CAST_APP/PS1",
                                                                    "PWK_WT",  "PWK_APP/PS1", "WSB_WT", "WSB_APP/PS1"))
cd11b.integrated$batch <- factor(cd11b.integrated$batch, levels = c("D", "C", "B", "A"))


sum_table <- cd11b.integrated@meta.data %>% group_by(seurat_clusters) %>% summarise(N=n(), ave_nCount_RNA=median(nCount_RNA), ave_nFeature_RNA=median(nFeature_RNA), ave_percent.mt=median(percent.mt))
prop.table(table(Idents(cd11b.integrated), cd11b.integrated$Group), margin = 2)
```

```
##     
##            B6_WT  B6_APP/PS1     CAST_WT CAST_APP/PS1      PWK_WT
##   0  0.230301054 0.194605845 0.234286865  0.193349220 0.152085453
##   1  0.138571894 0.113780245 0.266317486  0.177426343 0.327399118
##   2  0.223453972 0.124186974 0.158340048  0.102794627 0.131230926
##   3  0.177263341 0.149770185 0.095151759  0.104961005 0.061291963
##   4  0.087273122 0.068164079 0.047542305  0.018847487 0.103001017
##   5  0.006086295 0.003555633 0.087362342  0.042352686 0.120803662
##   6  0.007064450 0.128436389 0.005170561  0.172010399 0.001695490
##   7  0.052711662 0.055502558 0.014773033  0.030545927 0.016022380
##   8  0.017172047 0.018298500 0.026792909  0.017872617 0.024499830
##   9  0.021519400 0.026276992 0.016183186  0.015706239 0.033316378
##   10 0.019671775 0.036163386 0.022696750  0.025671577 0.019582909
##   11 0.016193892 0.031133466 0.023502552  0.035745234 0.008053577
##   12 0.002717096 0.050125748 0.001880204  0.062716638 0.001017294
##     
##      PWK_APP/PS1      WSB_WT WSB_APP/PS1
##   0  0.180369340 0.193769257 0.259896014
##   1  0.144194283 0.226292366 0.212478336
##   2  0.099418163 0.095971699 0.066828423
##   3  0.033392360 0.130320666 0.103570191
##   4  0.111307867 0.070295561 0.066620451
##   5  0.071211738 0.068013237 0.080623917
##   6  0.078800911 0.014835102 0.047209705
##   7  0.084745763 0.017459774 0.029254766
##   8  0.049203137 0.112746776 0.023986135
##   9  0.051479889 0.027958462 0.037781629
##   10 0.024791298 0.022252653 0.030641248
##   11 0.017961042 0.015748031 0.025303293
##   12 0.053124209 0.004336414 0.015805893
```

### plot both genotypes in all strains (all replicates combined)

```r
# generate meta data, combine cluster 0-5 into H for fraction plot (all replicates combined)
cd11b.meta <- cd11b.integrated@meta.data %>% 
  mutate(Genotype=factor(Genotype, levels=c("WT", "APP/PS1")),
         Strain=factor(Strain, levels = c("B6", "CAST", "PWK", "WSB")),
         new_clusters=ifelse(seurat_clusters %in% 0:5, "H", as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("H", "6", "7", "8", "9", "10", "11", "12"))) %>% 
  group_by(Strain, Genotype,new_clusters) %>% 
  summarise(N=n())

ggplot(cd11b.meta, aes(y=N, x=Genotype, fill= new_clusters)) + 
  geom_bar(stat="identity", position="fill", color="black") + 
  labs(y="Fraction", fill = "Clusters") +
  facet_grid(~ Strain) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = c("bold", "bold.italic")),
        axis.title.x = element_blank(), 
        strip.text.x = element_text(face = "bold"), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())
```

![](04b_microglia_statistics_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ggsave(paste(path, "fraction_WT_APPPS1.png", sep=""), width = 4.5, height = 4.5, units = "in")
```

### plot both genotypes in all strains (all replicates separated)

```r
### generate meta data, combine cluster 0-5 into H for fraction plot (all replicates separated)
cd11b.meta <- cd11b.integrated@meta.data %>% 
  mutate(Genotype=factor(Genotype, levels=c("WT", "APP/PS1")),
         Strain=factor(Strain, levels = c("B6", "CAST", "PWK", "WSB")),
         new_clusters=ifelse(seurat_clusters %in% 0:5, "H", as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("H", "6", "7", "8", "9", "10", "11", "12"))) %>% 
  group_by(rep, Strain, Genotype, Group, batch, new_clusters) %>% 
  arrange(Group) %>%
  summarise(N=n())

p <- ggplot(cd11b.meta, aes(y=N, x=batch, fill= new_clusters)) + 
  geom_bar(stat="identity", position=position_fill(), color="black") + 
  labs(y="Fraction", fill = "Clusters") +
  facet_grid(Group~., switch="y")+
  coord_flip()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(), 
        strip.text = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank(), 
        legend.position = "none")
ggsave(paste(path, "fraction_replicates_seperated.png", sep=""), p, width = 3.5, height = 5, units = "in")

# top to bottome: "B6_WT","B6_APP/PS1","CAST_WT", "CAST_APP/PS1", "PWK_WT",  "PWK_APP/PS1", "WSB_WT", "WSB_APP/PS1"
```

### Box plot for all microglia
#### generate meta data, combine cluster 0-5 into H for statistical testing and box plot

```r
cd11b.meta.stat <- cd11b.integrated@meta.data %>% 
  mutate(Genotype=factor(Genotype, levels=c("WT", "APP/PS1")),
         Strain=factor(Strain, levels = c("B6", "CAST", "PWK", "WSB")),
         new_clusters=ifelse(seurat_clusters %in% 0:5, "H", as.character(seurat_clusters)),
         new_clusters=factor(new_clusters, levels = c("H", "6", "7", "8", "9", "10", "11", "12"))) %>% 
  group_by(Strain, Genotype,rep, new_clusters) %>% 
  summarise(Med_nFeature=median(nFeature_RNA), 
            Med_percent_mt=median(percent.mt), 
            Med_percent_ribo=median(percent.ribo),
            N=n()) %>% 
  group_by(rep,Percent=N/sum(N)*100)
```



```r
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cd11b.meta.stat %>%
  ggplot(aes(y=Percent, x=Genotype, color=Genotype)) +
  geom_boxplot(outlier.size = 0, alpha=0.5) +
  geom_point(aes(color=Genotype), position=position_jitterdodge(), alpha=0.8) + 
  scale_colour_manual(values=cbPalette) + 
  theme_bw() +
  facet_grid(new_clusters ~ Strain, scales= "free_y")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 12), 
        axis.ticks.x = element_blank(), 
        axis.title= element_blank(), 
        #legend.text = element_text(face = c("plain", "italic")),
        legend.position = "bottom")
```

![](04b_microglia_statistics_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

```r
ggsave(paste(path, "cluster_box_all.png", sep=""), width = 4, height = 10, units = "in")
```

### Perform two-way ANOVA to determine the effect of strain and genotype on the percent of micrglia subclusters 

```r
####single function for stain and GT interaction debugging

clusters <- unique(cd11b.meta.stat$new_clusters) %>% as.list()

data = cd11b.meta.stat %>% filter(new_clusters %in% clusters[[1]])
aov_object = aov(Percent ~ Strain*Genotype , data=data)
aov.pvals = summary(aov_object)
aov.pvals= aov.pvals[[1]][5] %>% t() %>% as.data.frame()
names(aov.pvals) <- c("Strain", "Genotype", "Strain_Genotype", "Residuals")
aov.pvals <- aov.pvals %>% 
  select(-Residuals) %>% 
  mutate(Cluster = clusters[1] %>% as.character())

aov_StrainGT <- function(cluster, data){
  data = data %>% filter(new_clusters %in% cluster)
  aov_object = aov(Percent ~ Strain*Genotype, data=data)
  aov.pvals = summary(aov_object)
  aov.pvals= aov.pvals[[1]][5] %>% t() %>% as.data.frame()
  names(aov.pvals) <- c("Strain", "Genotype", "Strain_Genotype", "Residuals")
  aov.pvals <- aov.pvals %>% 
    select(-Residuals) %>% 
    mutate(Cluster = cluster %>% as.character())
  return(aov.pvals)
}

aov_StrainGT_object <- function(cluster, data){
  data = data %>% filter(new_clusters %in% cluster)
  aov_object = aov(Percent ~ Strain*Genotype, data=data)
  return(aov_object)
}

aov_StrainGT_table <- clusters %>% map_df(aov_StrainGT, data=cd11b.meta.stat)
aov_StrainGT_table <- aov_StrainGT_table %>% mutate_if(is.double, p.adjust)

aov_StrainGT_table$Cluster[aov_StrainGT_table$Strain_Genotype<0.05]
```

```
## [1] "H"  "6"  "7"  "12"
```

```r
# [1] "H"  "6"  "7"  "12"

aov_StrainGT_table$Cluster[aov_StrainGT_table$Strain<0.05]
```

```
## [1] "6"  "7"  "9"  "11" "12"
```

```r
# [1] "6"  "7"  "9" "11" "12" 

aov_StrainGT_table$Cluster[aov_StrainGT_table$Genotype<0.05]
```

```
## [1] "H"  "6"  "7"  "10" "11" "12"
```

```r
# [1] "6"  "7"  "H" "10" "11" "12"

# keep the annova object for  
aov_object_list <- clusters %>% map(aov_StrainGT_object, data=cd11b.meta.stat)
names(aov_object_list) <- clusters %>% unlist()
TukeyHSD(aov_object_list[["H"]]) %>% .$`Strain:Genotype` %>% data.frame(.,cluster="H")
```

```
##                                  diff        lwr         upr        p.adj
## CAST:WT-B6:WT              2.50564568  -8.247699  13.2589902 9.923753e-01
## PWK:WT-B6:WT               2.80327785  -8.811662  14.4182175 9.906220e-01
## WSB:WT-B6:WT              -6.85282267 -17.606167   3.9005219 4.246406e-01
## B6:APP/PS1-B6:WT         -20.94042199 -31.693767 -10.1870774 4.244607e-05
## CAST:APP/PS1-B6:WT       -22.62731099 -34.242251 -11.0123714 4.220632e-05
## PWK:APP/PS1-B6:WT        -22.67043094 -34.285371 -11.0554913 4.108533e-05
## WSB:APP/PS1-B6:WT         -7.57961006 -18.332955   3.1737345 3.075124e-01
## PWK:WT-CAST:WT             0.29763217 -11.317307  11.9125718 1.000000e+00
## WSB:WT-CAST:WT            -9.35846835 -20.111813   1.3948762 1.185954e-01
## B6:APP/PS1-CAST:WT       -23.44606767 -34.199412 -12.6927231 8.110079e-06
## CAST:APP/PS1-CAST:WT     -25.13295668 -36.747896 -13.5180171 9.096957e-06
## PWK:APP/PS1-CAST:WT      -25.17607662 -36.791016 -13.5611370 8.864693e-06
## WSB:APP/PS1-CAST:WT      -10.08525574 -20.838600   0.6680888 7.637308e-02
## WSB:WT-PWK:WT             -9.65610051 -21.271040   1.9588391 1.510316e-01
## B6:APP/PS1-PWK:WT        -23.74369984 -35.358639 -12.1287602 2.114176e-05
## CAST:APP/PS1-PWK:WT      -25.43058884 -37.847482 -13.0136961 2.057299e-05
## PWK:APP/PS1-PWK:WT       -25.47370878 -37.890602 -13.0568160 2.007033e-05
## WSB:APP/PS1-PWK:WT       -10.38288790 -21.997828   1.2320517 1.019553e-01
## B6:APP/PS1-WSB:WT        -14.08759932 -24.840944  -3.3342548 5.127222e-03
## CAST:APP/PS1-WSB:WT      -15.77448833 -27.389428  -4.1595487 3.557642e-03
## PWK:APP/PS1-WSB:WT       -15.81760827 -27.432548  -4.2026687 3.458418e-03
## WSB:APP/PS1-WSB:WT        -0.72678739 -11.480132  10.0265572 9.999978e-01
## CAST:APP/PS1-B6:APP/PS1   -1.68688900 -13.301829   9.9280506 9.996087e-01
## PWK:APP/PS1-B6:APP/PS1    -1.73000895 -13.344949   9.8849307 9.995384e-01
## WSB:APP/PS1-B6:APP/PS1    13.36081193   2.607467  24.1141565 8.546085e-03
## PWK:APP/PS1-CAST:APP/PS1  -0.04311994 -12.460013  12.3737728 1.000000e+00
## WSB:APP/PS1-CAST:APP/PS1  15.04770094   3.432761  26.6626406 5.724088e-03
## WSB:APP/PS1-PWK:APP/PS1   15.09082088   3.475881  26.7057605 5.565260e-03
##                          cluster
## CAST:WT-B6:WT                  H
## PWK:WT-B6:WT                   H
## WSB:WT-B6:WT                   H
## B6:APP/PS1-B6:WT               H
## CAST:APP/PS1-B6:WT             H
## PWK:APP/PS1-B6:WT              H
## WSB:APP/PS1-B6:WT              H
## PWK:WT-CAST:WT                 H
## WSB:WT-CAST:WT                 H
## B6:APP/PS1-CAST:WT             H
## CAST:APP/PS1-CAST:WT           H
## PWK:APP/PS1-CAST:WT            H
## WSB:APP/PS1-CAST:WT            H
## WSB:WT-PWK:WT                  H
## B6:APP/PS1-PWK:WT              H
## CAST:APP/PS1-PWK:WT            H
## PWK:APP/PS1-PWK:WT             H
## WSB:APP/PS1-PWK:WT             H
## B6:APP/PS1-WSB:WT              H
## CAST:APP/PS1-WSB:WT            H
## PWK:APP/PS1-WSB:WT             H
## WSB:APP/PS1-WSB:WT             H
## CAST:APP/PS1-B6:APP/PS1        H
## PWK:APP/PS1-B6:APP/PS1         H
## WSB:APP/PS1-B6:APP/PS1         H
## PWK:APP/PS1-CAST:APP/PS1       H
## WSB:APP/PS1-CAST:APP/PS1       H
## WSB:APP/PS1-PWK:APP/PS1        H
```

```r
# to export the statistic result: 
stat_list <- vector(mode = "list", length = length(clusters %>% unlist()))
names(stat_list) <- clusters %>% unlist()
for (i in clusters %>% unlist()){
  stat_list[[i]] <-  TukeyHSD(aov_object_list[[i]]) %>% 
    .$`Strain:Genotype` %>% 
    data.frame(.,cluster=i) %>% 
    rownames_to_column(var = "comparison")
}

stat_all <- do.call(rbind, stat_list)

# find genetype difference within strain comparison: 
stat_APP_cluster <- stat_all %>% 
  filter(comparison %in% c("B6:APP/PS1-B6:WT", "CAST:APP/PS1-CAST:WT", "PWK:APP/PS1-PWK:WT", "WSB:APP/PS1-WSB:WT")) %>% 
  mutate(Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(stat_APP_cluster, paste(path, "stat_GT_within_strain.txt", sep = ""), delim = "\t")

# find strain difference of WT (comparing to B6 in for WT in each strain)
stat_WT_strain_cluster <- stat_all %>% 
  filter(str_count(comparison, "WT")==2) %>% 
  mutate(Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(stat_WT_strain_cluster, paste(path, "stat_WT_between_strain.txt", sep = ""), delim = "\t") 

# find strain difference of APP/PS1 (comparing to B6 in for APP/PS1 in each strain)
stat_APP_strain_cluster <- stat_all %>% 
  filter(str_count(comparison, "APP/PS1")==2) %>% 
  mutate(Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(stat_APP_strain_cluster, paste(path, "stat_APP_between_strain.txt", sep = ""), delim = "\t") 
```


### check the statistics on nFeature, percent of mitochodria, and percent of ribosomal genes for each cluster

```r
## nFeature
aov_stat = aov(Med_nFeature ~ new_clusters, data=cd11b.meta.stat)
aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>% 
  mutate(comparison=paste(" ", comparison, sep = ""), Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(aov_table, path = paste(path, "Med_nFeature_comp_cluster.txt", sep = ""), delim = "\t")
filter(aov_table, Significance=="S")
```

```
##    comparison       diff        lwr         upr        p.adj Significance
## 1         6-H  1041.1724   796.0242  1286.32060 0.000000e+00            S
## 2         8-H -1153.9828 -1399.1309  -908.83457 0.000000e+00            S
## 3         9-H -1276.5862 -1521.7344 -1031.43802 0.000000e+00            S
## 4        11-H   747.2241   502.0759   992.37233 1.665335e-15            S
## 5        12-H -1121.1034 -1366.2516  -875.95526 0.000000e+00            S
## 6         7-6  -930.4138 -1175.5620  -685.26560 0.000000e+00            S
## 7         8-6 -2195.1552 -2440.3034 -1950.00698 0.000000e+00            S
## 8         9-6 -2317.7586 -2562.9068 -2072.61043 0.000000e+00            S
## 9        10-6  -831.6724 -1076.8206  -586.52423 0.000000e+00            S
## 10       11-6  -293.9483  -539.0965   -48.80009 7.249001e-03            S
## 11       12-6 -2162.2759 -2407.4241 -1917.12767 0.000000e+00            S
## 12        8-7 -1264.7414 -1509.8896 -1019.59319 0.000000e+00            S
## 13        9-7 -1387.3448 -1632.4930 -1142.19664 0.000000e+00            S
## 14       11-7   636.4655   391.3173   881.61371 2.650880e-12            S
## 15       12-7 -1231.8621 -1477.0103  -986.71388 0.000000e+00            S
## 16       10-8  1363.4828  1118.3346  1608.63095 0.000000e+00            S
## 17       11-8  1901.2069  1656.0587  2146.35509 0.000000e+00            S
## 18       10-9  1486.0862  1240.9380  1731.23440 0.000000e+00            S
## 19       11-9  2023.8103  1778.6622  2268.95853 0.000000e+00            S
## 20      11-10   537.7241   292.5759   782.87233 4.353195e-09            S
## 21      12-10 -1330.6034 -1575.7516 -1085.45526 0.000000e+00            S
## 22      12-11 -1868.3276 -2113.4758 -1623.17940 0.000000e+00            S
```

```r
## percent of mitochodria
aov_stat = aov(Med_percent_mt ~ new_clusters, data=cd11b.meta.stat)
aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>% 
  mutate(comparison=paste(" ", comparison, sep = ""),Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(aov_table, path = paste(path, "Med_percent_mt_comp_cluster.txt", sep = ""), delim = "\t")
filter(aov_table, Significance=="S")
```

```
##    comparison       diff         lwr        upr        p.adj Significance
## 1         8-H  1.5357783  1.12890765  1.9426489 0.000000e+00            S
## 2         9-H  1.4652269  1.05835624  1.8720975 0.000000e+00            S
## 3        12-H  0.5280024  0.12113175  0.9348730 2.418138e-03            S
## 4         8-6  1.1877000  0.78082939  1.5945707 9.436896e-15            S
## 5         9-6  1.1171486  0.71027798  1.5240193 1.495470e-13            S
## 6         8-7  1.4443520  1.03748139  1.8512227 0.000000e+00            S
## 7         9-7  1.3738006  0.96692998  1.7806713 0.000000e+00            S
## 8        12-7  0.4365761  0.02970549  0.8434468 2.580750e-02            S
## 9        10-8 -1.4571625 -1.86403315 -1.0502919 0.000000e+00            S
## 10       11-8 -1.2688316 -1.67570229 -0.8619610 6.883383e-15            S
## 11       12-8 -1.0077759 -1.41464653 -0.6009053 2.557010e-11            S
## 12       10-9 -1.3866111 -1.79348174 -0.9797405 0.000000e+00            S
## 13       11-9 -1.1982802 -1.60515088 -0.7914096 6.106227e-15            S
## 14       12-9 -0.9372245 -1.34409513 -0.5303538 6.195263e-10            S
## 15      12-10  0.4493866  0.04251597  0.8562573 1.907217e-02            S
```

```r
## percent of ribosomal
aov_stat = aov(Med_percent_ribo ~ new_clusters, data=cd11b.meta.stat)
aov_table <- TukeyHSD(aov_stat) %>% .$new_clusters %>% data.frame() %>% rownames_to_column(var = "comparison") %>% 
  mutate(comparison=paste(" ", comparison, sep = ""),Significance=ifelse(p.adj<0.05, "S", "NS"))
write_delim(aov_table, path = paste(path, "Med_percent_ribo_comp_cluster.txt", sep = ""), delim = "\t")
filter(aov_table, Significance=="S")
```

```
##    comparison       diff        lwr       upr        p.adj Significance
## 1         6-H   5.017224   3.472120  6.562327 6.661338e-16            S
## 2         8-H  -5.343313  -6.888416 -3.798209 0.000000e+00            S
## 3         9-H   6.706032   5.160929  8.251136 0.000000e+00            S
## 4        12-H  12.460227  10.915123 14.005330 0.000000e+00            S
## 5         7-6  -4.286663  -5.831766 -2.741559 9.137135e-14            S
## 6         8-6 -10.360536 -11.905640 -8.815433 0.000000e+00            S
## 7         9-6   1.688808   0.143705  3.233912 2.132229e-02            S
## 8        10-6  -5.116565  -6.661668 -3.571462 2.997602e-15            S
## 9        11-6  -4.585704  -6.130807 -3.040600 6.661338e-15            S
## 10       12-6   7.443003   5.897899  8.988106 0.000000e+00            S
## 11        8-7  -6.073874  -7.618977 -4.528770 0.000000e+00            S
## 12        9-7   5.975471   4.430368  7.520574 0.000000e+00            S
## 13       12-7  11.729665  10.184562 13.274769 0.000000e+00            S
## 14        9-8  12.049345  10.504241 13.594448 0.000000e+00            S
## 15       10-8   5.243971   3.698868  6.789075 0.000000e+00            S
## 16       11-8   5.774833   4.229729  7.319936 0.000000e+00            S
## 17       12-8  17.803539  16.258436 19.348643 0.000000e+00            S
## 18       10-9  -6.805373  -8.350477 -5.260270 0.000000e+00            S
## 19       11-9  -6.274512  -7.819616 -4.729409 0.000000e+00            S
## 20       12-9   5.754194   4.209091  7.299298 0.000000e+00            S
## 21      12-10  12.559568  11.014465 14.104671 0.000000e+00            S
## 22      12-11  12.028707  10.483603 13.573810 0.000000e+00            S
```

## Pseudotime analysis (diffusion map)
#### too many cells for diffusion map, need sampling


```r
library(destiny)
library(SingleCellExperiment)
library(scater)

cd11b.integrated$final_clusters <-  ifelse(cd11b.integrated$seurat_clusters %in% 0:5,"H",
                                           cd11b.integrated$seurat_clusters %>% as.character())

sampling <- cd11b.integrated@meta.data %>% 
  rownames_to_column(var = "cell_ID") %>% 
  group_by(Group) %>% 
  sample_n(1000) # take 1000 random cells from each group

mg.small <- subset(cd11b.integrated, cells=sampling$cell_ID)

mg.small <- as.SingleCellExperiment(mg.small)

# use diffusion map to calculate pseudotime
pca <- reducedDim(mg.small)
cellLabels <- mg.small$seurat_clusters

pca_tidy <- as.data.frame(pca) %>% rownames_to_column()

rownames(pca) <- cellLabels

dm <- DiffusionMap(pca)

dpt <- DPT(dm) #

mg.small$pseudotime_dpt <- rank(dpt$dpt) 

df <- colData(mg.small) %>% as.data.frame()

df$final_clusters<- ifelse(df$seurat_clusters %in% 0:5,"H", df$seurat_clusters %>% as.character())

ggplot(df, aes(pseudotime_dpt, fill=final_clusters)) +
  geom_histogram(binwidth = 100, color="grey",size=0.1)+ 
  facet_grid(Group~., switch="y")+ 
  scale_y_continuous("count", position="right") +
  labs(x="DAM <-  pseudotime  -> Homeostatic")+
  theme_bw()+
  theme(text = element_text(family = "Arial", size = 10),
        strip.text.y = element_text(size=5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y.right = element_blank(),
        legend.position = "null")
```

![](04b_microglia_statistics_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
ggsave(paste(path, "pseudotime.png", sep = ""), width = 3.5 , height = 5.3, units = "in", dpi = 600)
```




