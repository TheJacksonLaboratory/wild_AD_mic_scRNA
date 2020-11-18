##

## 

# plot individual genes in violin plot
my_vlnPlot <- function(gene, object, idents=NULL, title = FALSE, pt.size=0, slot="data"){
  if (title){
    VlnPlot(object = object, features = gene, pt.size=pt.size, idents = idents , slot = slot, ncol=1, log = TRUE) + 
      NoLegend() + 
      theme(axis.title = element_blank())
  }else{
    VlnPlot(object = object, features = gene, pt.size=pt.size, idents = idents , slot = slot, ncol=1, log = TRUE) + 
      NoLegend() + 
      theme(title = element_blank(), 
            axis.text.x = element_text(face="bold"),
            axis.text = element_text(size=15),
            axis.line.x.top = element_line(colour = "black"),
            axis.ticks.x = element_blank())
  }
}


## Violin plot
my_vlnPlot2 <- function(gene, object, idents=NULL, pt.size=0, slot="data", legend=FALSE){
  if(legend){
    VlnPlot(object = object, features = gene, pt.size=pt.size, idents = idents , slot = slot, ncol=1, log = TRUE) + 
      #NoLegend() + 
      theme(title = element_blank(), 
            axis.text.x = element_blank(),
            axis.text = element_text(size=15),
            axis.line.x.top = element_line(colour = "black"),
            axis.ticks.x = element_blank(),
            legend.position = "bottom")
  }else{
    VlnPlot(object = object, features = gene, pt.size=pt.size, idents = idents , slot = slot, ncol=1, log = TRUE) + 
      NoLegend() + 
      theme(title = element_blank(), 
            axis.text.x = element_blank(),
            axis.text = element_text(size=15),
            axis.line.x.top = element_line(colour = "black"),
            axis.ticks.x = element_blank())
  }
}



