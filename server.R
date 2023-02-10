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
    p4 <- ggplot(gE, aes(factor(Sample), gE[,gene], fill = Sample))+facet_wrap(~gE[,cell_group]) + ylab(paste0(gene," expression"))+  geom_split_violin() +scale_fill_manual(values=c(alpha("#424242",0.75), alpha("#B40431",0.75)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),panel.background = element_blank(),legend.position = "top",legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),legend.title = element_text(size=10),legend.text = element_text(size=10)) + stat_summary(fun.y = median, geom='point', shape = 21, size = 3) + scale_x_discrete(name ="")+ scale_color_manual(values=c("black", "black"))
    #+ geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE)
    design <- "
  1##
"
    p4+plot_layout(design=design)
    
    
    
  })
    
  
}