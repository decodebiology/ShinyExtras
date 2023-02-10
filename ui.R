library(shiny, lib.loc = "/Users/subhash/R/")
library(Seurat, lib.loc = "/Users/subhash/R/") #BiocManager::install("Seurat", lib="/Users/subhash/R", force=T)
library(dittoSeq, lib.loc = "/Users/subhash/R/") #BiocManager::install("dittoSeq", lib="/Users/subhash/R", force=T)
library(ggplot2, lib.loc = "/Users/subhash/R/") #BiocManager::install("ggplot2", lib="/Users/subhash/R", force=T)
library(htmltools, lib.loc = "/Users/subhash/R/") 
library(patchwork)
#library(introdataviz, lib.loc = "/Users/subhash/R/") #geom_split_violin pre-packaged in this
SeuObj <- readRDS("/Users/subhash/iCloud Drive (Archive)/Documents/Projects/shiny_app/dataset/rds/mouse_AA_CREB_final_CellTypes.Rds")

#.libPaths( c("/Users/subhash/R", .libPaths()) )

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Split Violin/ Rain cloud Plots!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      selectizeInput(
        inputId = "gene", 
        selected = "Apoe",
        label = "Input gene names",
        multiple = FALSE,
        choices = c("gene" = "",rownames(SeuObj)),
        options = list(
          create = TRUE,
          placeholder = "Apoe",
          maxItems = '1'
          #onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
          #onType = I("function (str) {if (str === \"\") {this.close();}}")
        )
      ),
      selectizeInput(
        inputId = "cell_group", 
        selected = "Cell_type",
        label = "Input cell group",
        multiple = FALSE,
        choices = c("cell_group" = "",c("Cell_type","Cell_compartments")),
        options = list(
          create = TRUE,
          placeholder = "Cell_type",
          maxItems = '1'
          #onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
          #onType = I("function (str) {if (str === \"\") {this.close();}}")
        )
      ),
      checkboxGroupInput("samples",  label = "Select samples (max. 2)",
                         unique(SeuObj@meta.data$Sample), selected = c(unique(SeuObj@meta.data$Sample)[1:2])
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Histogram ----
      #plotOutput(outputId = "distPlot"),
      # Output: UMAP ----
      #plotOutput(outputId = "umapPlot"),
      # Output: Slitlot ----
      #plotOutput(outputId = "splitPlot"),
      # Output: Cell_type ----
      plotOutput(outputId = "ctPlot")
    )
  ),
 
)

