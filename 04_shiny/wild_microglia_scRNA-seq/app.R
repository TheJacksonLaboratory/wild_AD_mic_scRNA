# setwd("~/code/2019_04_scRNA_CD11b/shiny/wild_microglia_scRNA-seq/")
library(shiny)
library(shinyWidgets)
library(Seurat)
library(tidyverse)
library(config)
library(patchwork)


#cd11b.integrated <- readRDS("~/code/2019_04_scRNA_CD11b/emase_2/output/Integration_mg/mg_int.rds") 
cd11b.integrated <- readRDS(config::get("myfile"))
DE_vsB6 <- read_delim("data/DE_vsB6_coef_cluster.txt", delim = "\t", col_types = "cccdcdddddcc") %>% 
  select(-Symbol, -logFC, -logCPM, -F) %>% 
  rename(Symbol="Orig_Symbol") %>% 
  mutate(FDR = formatC(FDR, format = "e", digits = 2),
         PValue = formatC(PValue, format = "e", digits = 2))

# the_file <- config::get(“myfile”)
# 
# readr::read_tsv(the_file)


DefaultAssay(cd11b.integrated) <- "RNA"
cd11b.integrated$Strain <- str_replace_all(cd11b.integrated$Strain, pattern = "B6J", replacement = "B6")
cd11b.integrated$Group <- str_replace_all(cd11b.integrated$Group, pattern = "B6J", replacement = "B6")
cd11b.integrated$Group <- factor(cd11b.integrated$Group, levels = c("B6_WT","B6_APP/PS1","CAST_WT", "CAST_APP/PS1", 
                                                                    "PWK_WT",  "PWK_APP/PS1", "WSB_WT", "WSB_APP/PS1"))
cd11b.integrated$final_clusters <- ifelse(cd11b.integrated$seurat_clusters %in% 0:5, "H", as.character(cd11b.integrated$seurat_clusters))
cd11b.integrated$final_clusters <- factor(cd11b.integrated$final_clusters, levels = c("H", "6", "7", "8", "9", "10", "11", "12"))

Idents(cd11b.integrated) <- "Group"
group <- levels(cd11b.integrated$Group)
split <- vector(mode = "list", length = length(group))
for (i in seq_along(group)){
  split[[i]] <- cd11b.integrated %>% subset(idents=group[i])
  Idents(split[[i]]) <- "final_clusters"
  # subset(subset=Group==group[i]) doesn't work
}



ui <- fluidPage(
  # Application title
  titlePanel("Microglia in wild-derived APP/PS1 strains"),
  
  # Sidebar with a select input for gene name 
  sidebarLayout(
    sidebarPanel(
      selectInput("gene",
                  "Select a gene",
                  unique(cd11b.integrated %>% GetAssay() %>% rownames()),
                  selected = "NULL"
      ),
      actionBttn("button", "Search")
    ),
    
    # Show a plot of the UMAP plot of the gene
    mainPanel(
      tabsetPanel(
        tabPanel("UMAP Plot", plotOutput("UMAP", height = "1600px")),
        tabPanel("Ridge Plot", plotOutput("Ridge", height = "1500px", width = "500px")),
        tabPanel("Violin Plot", plotOutput("Violin", height = "700px", width = "600px")),
        tabPanel("DE genes (wild vs B6)", DT::DTOutput('DE_table'))
      )
    )
  )
)


server <- function(input, output) {
  
  my_FeaturePlot <- function(object, group, feature){
    FeaturePlot(object = object, features = feature, label = TRUE, repel = TRUE) +
      coord_fixed()+
      theme(axis.title = element_blank(),
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5)) +
      ggtitle(group)
  }
  
  my_RidgePlot <- function(object, group, feature){
    RidgePlot(object = object, features = feature) +
      ggtitle(group) +
      theme(axis.title = element_blank(), 
            legend.position = "none", 
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  }
  
  my_Vln <- function(object, group, feature){
    VlnPlot(object = object, features = feature, pt.size=0, ncol=1, log = TRUE) + 
      NoLegend() + 
      ggtitle(group) +
      theme(axis.title = element_blank(), 
            legend.position = "none", 
            plot.title = element_text(face = "bold", size = 12, hjust = 0.5))
  }
  
  plot_fe <- eventReactive(input$button, {map2(split, group, my_FeaturePlot, feature=input$gene)})
  
  output$UMAP <- renderPlot({
    (plot_fe()[[1]]+plot_fe()[[2]])/(plot_fe()[[3]]+plot_fe()[[4]])/(plot_fe()[[5]]+plot_fe()[[6]])/(plot_fe()[[7]]+plot_fe()[[8]])+
    plot_annotation(title = input$gene, 
                    theme = theme(plot.title = element_text(size=16, face = "bold.italic", hjust = 0.5)))
  })
  
  plot_ri <- eventReactive(input$button, {map2(split, group, my_RidgePlot, feature=input$gene)})
  
  output$Ridge <- renderPlot({
    (plot_ri()[[1]]+plot_ri()[[2]])/(plot_ri()[[3]]+plot_ri()[[4]])/(plot_ri()[[5]]+plot_ri()[[6]])/(plot_ri()[[7]]+plot_ri()[[8]])+
      plot_annotation(title = input$gene, 
                      theme = theme(plot.title = element_text(size=16, face = "bold.italic", hjust = 0.5)))
  })
  
  plot_v <- eventReactive(input$button, {map2(split, group, my_Vln, feature=input$gene)})
  
  output$Violin <- renderPlot({
    (plot_v()[[1]]+plot_v()[[2]])/(plot_v()[[3]]+plot_v()[[4]])/(plot_v()[[5]]+plot_v()[[6]])/(plot_v()[[7]]+plot_v()[[8]]) +
      plot_annotation(title =input$gene, 
                      theme = theme(plot.title = element_text(size=16, face = "bold.italic", hjust = 0.5)))
  })
  
  
  tb <- eventReactive(input$button, {DE_vsB6 %>% filter(Symbol==input$gene)})
  output$DE_table <- DT::renderDT({
     tb()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

