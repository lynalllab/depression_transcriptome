#leading edge plot 
#adapted from mel
library(ggfittext)

# Function to make and save escarpment plot
ledge_plot <- function(gseaout, pathway, sample_extension, ymin, ymax, xmax, width=7, height=5, plot = T){
  if (pathway %in% gseaout@result$ID){
    leading <- gsub("/", " ", gseaout@result[gseaout@result$Description==pathway,"core_enrichment"])
    nes <- gseaout@result[gseaout@result$Description==pathway,"NES"]
    p <- gseaplot(gseaout, pathway, title=paste(pathway, "  NES =", round(nes,3), sep = " "), by="runningScore") + expand_limits(y=c(ymin, 0.1)) + 
      geom_fit_text(data=data.frame(x=1,y=1), aes(ymin = ymin, ymax = ymax, xmin = 0, xmax = xmax, label = paste0("Leading edge:\n",leading)), reflow=TRUE) +
      labs(subtitle = paste(ymin, ymax, sep = " ")) +
      theme(plot.subtitle = element_text(colour = "white"))
    if(plot == TRUE){
      ggsave(p, filename = paste0(sample_extension,"_leading_edge_",pathway,".pdf"), width = width, height = height)
      return(p)
    } else {
      return(p)
    }
  } else {
    print(paste0("Pathway ",pathway," not in gseaout"))
  }
}
