# supplementary table showing how GWAS gene overlap with mouse scRNA-seq DE genes

library(tidyverse)
library(ggrepel)

output <- "11_GWAS_gene_table/"

# make overlaping gene with GWAS genes
# the pool of GWAS gene is from 
# paper1 Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk https://www.nature.com/articles/s41588-018-0311-9#Tab1  
# Table 1 Summary statistics of significantly associated regions identified in the genome-wide association analysis of AD case-control status, AD-by-proxy phenotype, and meta-analysis
# downloaded to "GWAS/Jansen_GWAS_table1.xlsx"

# paper 2  "table 2 from this paper: https://www.medrxiv.org/content/10.1101/19012021v1.full.pdf"

library(readxl)
GWAS1 <- read_excel(path = "GWAS/Jansen_GWAS_table1.xlsx", skip = 1) 
GWAS1 <- GWAS1 %>%
  mutate(Locus=as.integer(Locus)) %>% 
  #filter(Bolded=="N") %>% 
  #drop_na(Locus) %>% 
  rename(SYMBOL="Gene") %>% 
  select("SYMBOL") %>% 
  mutate(Source="Jansen et al.")

GWAS2 <- read_excel(path = "GWAS/19012021v1.full 27.xlsx", sheet = 2)
GWAS2 <- data.frame(SYMBOL=GWAS2$Gene %>% unique())
GWAS2$Source <- "de Rojas et al."

GWAS <- GWAS1 %>% full_join(GWAS2, by="SYMBOL")
GWAS <- GWAS %>% 
  mutate(Source=ifelse(!is.na(Source.x)&!is.na(Source.y), paste(Source.x, Source.y, sep = " & "),
                            ifelse(is.na(Source.x), Source.y, Source.x)))

GWAS <- GWAS %>% 
  select(-Source.x:-Source.y) %>% 
  mutate(symbol_mouse=str_to_title(SYMBOL))

# Find whether GWAS gene expressed in mouse
load("01_DE_gene/y_H_cluster.rda")
mouse_gene_ID<- data.frame(Ensemble_ID=rownames(y))
EnsembleID <- read_tsv("../documents/Ensembl2Symbol.tsv")
names(EnsembleID) <- c("Ensemble_ID", "Symbol")
mouse_gene_ID <- mouse_gene_ID %>% left_join(EnsembleID)

GWAS <- GWAS %>% left_join(mouse_gene_ID, by = c("symbol_mouse"="Symbol"))
GWAS <- GWAS %>% mutate(Overlap_with_mouse_scRNAseq=ifelse(is.na(Ensemble_ID), "No", "Yes"))

GWAS_x <- GWAS %>% left_join(EnsembleID, by = c("symbol_mouse"="Symbol")) # to check what genes are non-homologous in mice

GWAS <- GWAS %>% select(-Ensemble_ID) %>% arrange(SYMBOL)
GWAS <- GWAS %>% arrange(SYMBOL)

write_delim(GWAS, paste(output, "GWAS_overlap_mouse_table.txt"), delim = "\t")


## load GWAS gene DE gene intersection term for strain effect
load("10_GWAS_heatmap/gene_table_intersection_cluster_combined.rda")
#df <- select(df, Orig_Symbol, Intersections, Cluster)
#GWAS <- GWAS %>% left_join(df, by = c("symbol_mouse"="Orig_Symbol"))
