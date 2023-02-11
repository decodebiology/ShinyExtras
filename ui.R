library(shiny, lib.loc = "/Users/subhash/R/")
library(Seurat, lib.loc = "/Users/subhash/R/") #BiocManager::install("Seurat", lib="/Users/subhash/R", force=T)
library(dittoSeq, lib.loc = "/Users/subhash/R/") #BiocManager::install("dittoSeq", lib="/Users/subhash/R", force=T)
library(ggplot2, lib.loc = "/Users/subhash/R/") #BiocManager::install("ggplot2", lib="/Users/subhash/R", force=T)
library(htmltools, lib.loc = "/Users/subhash/R/") 
library(patchwork)
#library("Signac", lib.loc = "/Users/subhash/R")
library(introdataviz, lib.loc = "/Users/subhash/R/") #geom_split_violin pre-packaged in this
SeuObj <- readRDS("/Users/subhash/iCloud Drive (Archive)/Documents/Projects/shiny_app/dataset/rds/mouse_AA_CREB_final_CellTypes.Rds")
selected_markers <- "Ace2,Apob,Gda,Reg3g,Alpi,Car4,Ccl25,Pycard,Krt19,Gpx1,Lgr5,Olfm4,Slc12a2,Ascl2,Ly6a,S100a6,Fcgbp,Agr2,Clca1,Muc2,Spink4,Pcna,Mki67,Ube2c,Cenpm,Tk1,Ube2t,Tacc3,Neurod1,Chga,Chgb,Cck,Gip,Hck,Lrmp,Sh2d6,Cd24a,Klre1,Klrd1,Lag3,Trdc,Gzmk,Cd8a,Tcf7,Cd8a,Gzmk,Trdc,Cd79a,Ms4a1,Cd79b,Ighkc,Alpi,AY036118,Lars2"

#.libPaths( c("/Users/subhash/R", .libPaths()) )

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Split Violin/ Group  DotPlots!"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      selectizeInput(
        inputId = "gene", 
        selected = "Apoe",
        label = "1. Input gene name for split violins (max. 1):",
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
        label = "Input cell group:",
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
      checkboxGroupInput("samples",  label = "2. Select samples (max. 2 for violins & any numbers for dotplot):",
                         unique(SeuObj@meta.data$Sample), selected = c(unique(SeuObj@meta.data$Sample)[1:2])
      ),
      textAreaInput(
        inputId = "genelist",
        label = "3. Gene list for dotplot (comma separated list) & Select samples from section 2:",
        value = selected_markers,
        placeholder = "List of genes",
        resize = "vertical"
      ),
      selectizeInput(
        inputId = "ctype", 
        selected = "Enterocyte",
        label = "Input cell type name",
        multiple = FALSE,
        choices = c("ctype" = "",levels(SeuObj@meta.data$Cell_type)),
        options = list(
          create = TRUE,
          placeholder = "Enterocyte",
          maxItems = '1'
          #onDropdownOpen = I("function($dropdown) {if (!this.lastQuery.length) {this.close(); this.settings.openOnFocus = false;}}"),
          #onType = I("function (str) {if (str === \"\") {this.close();}}")
        )
      ),
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      plotOutput(outputId = "ctPlot"),
      br(),
      br(),
      br(),
      br(),
      plotOutput(outputId = "dotplot"),

   
    )
  ),
 
)

