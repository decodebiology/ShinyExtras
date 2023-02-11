# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  
  GeomSplitViolin <- ggplot2::ggproto(
    "GeomSplitViolin",
    ggplot2::GeomViolin,
    draw_group = function(self,
                          data,
                          ...,
                          # add the nudge here
                          nudge = 0,
                          draw_quantiles = NULL) {
      data <- transform(data,
                        xminv = x - violinwidth * (x - xmin),
                        xmaxv = x + violinwidth * (xmax - x))
      grp <- data[1, "group"]
      newdata <- plyr::arrange(transform(data,
                                         x = if (grp %% 2 == 1) xminv else xmaxv),
                               if (grp %% 2 == 1) y else -y)
      newdata <- rbind(newdata[1, ],
                       newdata,
                       newdata[nrow(newdata), ],
                       newdata[1, ])
      newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
      
      # now nudge them apart
      #newdata$x <- ifelse(newdata$group %% 2 == 1,
      #                   newdata$x - nudge,
      #                  newdata$x + nudge)
      
      if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
        
        stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
        
        quantiles <- ggplot2:::create_quantile_segment_frame(data,
                                                             draw_quantiles)
        aesthetics <- data[rep(1, nrow(quantiles)),
                           setdiff(names(data), c("x", "y")),
                           drop = FALSE]
        aesthetics$alpha <- rep(1, nrow(quantiles))
        both <- cbind(quantiles, aesthetics)
        quantile_grob <- ggplot2::GeomPath$draw_panel(both, ...)
        ggplot2:::ggname("geom_split_violin",
                         grid::grobTree(ggplot2::GeomPolygon$draw_panel(newdata, ...),
                                        quantile_grob))
      }
      else {
        ggplot2:::ggname("geom_split_violin",
                         ggplot2::GeomPolygon$draw_panel(newdata, ...))
      }
    }
  )
  
  geom_split_violin <- function(mapping = NULL,
                                data = NULL,
                                stat = "ydensity",
                                position = "identity",
                                # nudge param here
                                #nudge = 0,
                                ...,
                                draw_quantiles = NULL,
                                trim = TRUE,
                                scale = "area",
                                na.rm = FALSE,
                                show.legend = NA,
                                inherit.aes = TRUE) {
    
    ggplot2::layer(data = data,
                   mapping = mapping,
                   stat = stat,
                   geom = GeomSplitViolin,
                   position = position,
                   show.legend = show.legend,
                   inherit.aes = inherit.aes,
                   params = list(trim = trim,
                                 scale = scale,
                                 # don't forget the nudge
                                 #nudge = nudge,
                                 draw_quantiles = draw_quantiles,
                                 na.rm = na.rm,
                                 ...))
  }
  
  SeuObj <- readRDS("/Users/subhash/iCloud Drive (Archive)/Documents/Projects/shiny_app/dataset/rds/mouse_AA_CREB_final_CellTypes.Rds")
 
  output$ctPlot <- renderPlot({
    cell_group <- input$cell_group
    gene <- input$gene
    samples <- input$samples
    gEx <- FetchData(SeuObj, vars =  c("Ptprc",gene), slot="scale.data", assay="RNA")
    gE <- merge(gEx,SeuObj@meta.data,by="row.names")
    gE <- data.frame(gE)
    gE <- subset(gE,gE$Sample %in% c(samples))
    gE$Sample <- droplevels(factor(gE$Sample))
    gene <- gsub("-",".",gene)
    p1 <- ggplot(gE, aes(x="", gE[,gene], fill = Sample))+ scale_color_manual(values=c("black", "black"))+facet_wrap(~gE[,cell_group])+scale_fill_manual(values=c(alpha("#01A9DB",0.4), alpha("#B40431",0.4))) + ylab(paste0(gene," expression"))+  introdataviz::geom_split_violin( alpha = .4, trim = FALSE, scale="width") + theme(strip.text.x = element_text(size = 11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.text.y = element_text(size=10,color="black"),axis.title.y = element_text(size=12,color="black"),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "top",legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),legend.title = element_text(size=12),legend.text = element_text(size=12)) + stat_summary(fun.y = median, geom='point', shape = 21, size = 3) + scale_x_discrete(name ="")
    #
    #+scale_fill_brewer(palette = "Dark2")+ geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE)

    print(p1)
    
    
    
  },alt="Please select appropriate parameters from leftpanel",res=72, width=500, height=500)
 # output$dotplot <- renderPlot({
#    genes <- input$genelist
 #   ctype1 <- input$ctype
  #  sample_list <- input$samples
  #  glist <- c(strsplit(genes, ","))
  #  SeuObj_subset <- SeuObj
  #  Idents(SeuObj_subset) <- SeuObj_subset@meta.data$Sample
  #  SeuObj_subset1 <- subset(SeuObj_subset, subset = Cell_type == ctype1)
  #  p2 <- DotPlot(SeuObj_subset1, group.by="Sample", idents=c(sample_list),features=unique(unlist(glist)) )+coord_flip()+ggtitle(paste0(ctype1))+ scale_color_gradient2(low = "#4B088A", mid = "white",high = "#8A0829",midpoint = 0,name="Expression") + RotatedAxis()+ theme(axis.text.x=element_text(size=rel(1), angle=90, hjust = 1, vjust = 0.5), axis.text.y=element_text(size=rel(1), face="bold.italic") )

   # p2+plot_layout()

    
  #})
  plotInput <- reactive({
    genes <- input$genelist
    ctype1 <- input$ctype
    sample_list <- input$samples
    glist <- c(strsplit(genes, ","))
    SeuObj_subset <- SeuObj
    Idents(SeuObj_subset) <- SeuObj_subset@meta.data$Sample
    SeuObj_subset1 <- subset(SeuObj_subset, subset = Cell_type == ctype1)
    p2 <- DotPlot(SeuObj_subset1, group.by="Sample", idents=c(sample_list),features=unique(unlist(glist)) )+ggtitle(paste0(ctype1))+ scale_color_gradient2(low = "#4B088A", mid = "white",high = "#8A0829",midpoint = 0,name="Expression") + RotatedAxis()+ theme(axis.text.y=element_text(size=rel(1)), axis.text.x=element_text(size=rel(1), face="bold.italic", angle=90, hjust = 1, vjust = 0.5) )
    #+coord_flip()
    #p2
    
    
  })
  output$dotplot <- renderPlot({
    print(plotInput())
  },res=72, width=850, height=250)
  
 
}
