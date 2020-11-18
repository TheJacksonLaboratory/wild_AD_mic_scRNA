
what_dims <- function(cd11b.strain, path, strain, round){
  cd11b.strain <- AddModuleScore( object = cd11b.strain, features = ribo.genes, ctrl = 100, name = 'ribo_Features')
  
  cd11b.strain <- cd11b.strain %>% 
    NormalizeData() %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
    ScaleData(vars.to.regress = c("batch", "ribo.genes", "percent.mt", "nFeature_RNA")) %>% 
    RunPCA()
  
  cd11b.strain <- JackStraw(cd11b.strain, num.replicate = 30, dim=30)
  cd11b.strain <- ScoreJackStraw(cd11b.strain, dims = 1:30)
  ElbowPlot(cd11b.strain, ndims = 30) + ggtitle(label = paste(strain, round, sep=" "))
  ggsave(paste(path, strain, "_", round, "_", "ElbowPlot",  ".png", sep=""), units = "in", width = 7, height = 4,  dpi=150)
  
  print(cd11b.strain[["pca"]], dims = 1:30, nfeatures = 30)
  return(cd11b.strain)
}


# QC: violin plot of nfeatures_RNA, percent.mt, percent.ribo for each cluster
QC_plot <-function(data, y){
  p <- ggplot(data, aes_string(x="seurat_clusters", y=y, color="seurat_clusters")) +
    geom_violin() +
    geom_boxplot(width=0.07, outlier.shape = NA, color = "black", alpha=0.7) +
    theme_bw()+
    theme(legend.position = "none", axis.title.x = element_blank())
  return(p)
}


markers <- function(cd11b.strain, path, strain, round, res){
  # check dimension reduction
  DimPlot(cd11b.strain, reduction = "umap", label = TRUE, pt.size = 0.001) + 
    ggtitle(label = strain) + coord_fixed()
  ggsave(paste(path, strain, "_", round, "_", res, "_", "DimPlot1",  ".png", sep=""), units = "in", width = 7.3, height = 7,  dpi=150)
  
  DimPlot(cd11b.strain, reduction = "umap", label = FALSE, group.by="batch", pt.size = 0.001) + 
    ggtitle(label = strain) + coord_fixed()
  ggsave(paste(path, strain, "_", round, "_", res, "_", "DimPlot2",  ".png", sep=""), units = "in", width = 7.3, height = 7,  dpi=150)
  
  DimPlot(cd11b.strain, reduction = "umap", label = TRUE, pt.size = 0.001, split.by = "Genotype")+ ggtitle(label = strain)+ coord_fixed()
  ggsave(paste(path, strain, "_", round, "_", res, "_", "DimPlot3",  ".png", sep=""), units = "in", width = 8.6, height = 4.7,  dpi=150)
  
  DimPlot(cd11b.strain, reduction = "pca", label = TRUE, pt.size = 0.001) + 
    ggtitle(label = strain)
  ggsave(paste(path, strain, "_", round, "_", res, "_", "DimPlot4",  ".png", sep=""), units = "in", width = 10, height = 5,  dpi=150)
  
  # QC: violin plot of nfeatures_RNA, percent.mt, percent.ribo for each cluster
  # funcion: QC_plot(data, y) from source
  p_QC <- c("nFeature_RNA", "percent.mt", "percent.ribo") %>% map(~QC_plot(cd11b.strain@meta.data, .))
  p <- plot_grid(plotlist=p_QC, ncol=1, align="hv")
  title <- ggdraw() + draw_label(paste(strain, round, res, "QC", sep=" "), fontface='bold')
  plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
  ggsave(paste(path, strain, "_", round, "_", res, "_", "QC",  ".png", sep=""), units = "in", width = 10, height = 5,  dpi=150)
  
  # Find cluster markers
    cd11b.markers <- FindAllMarkers(cd11b.strain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident=300) #max.cells.per.ident
    cd11b.markers <- cd11b.markers %>% rownames_to_column(var="symbol")
  
  # save cell metadata and marker info into rda
    meta <- cd11b.strain@meta.data %>% select(-starts_with("ribo_"))
    save(meta, cd11b.markers, file=paste(path, strain, "_", round, "_", res, "_", "Meta_Marker.rda", sep=""))
    return(cd11b.markers)
  }

