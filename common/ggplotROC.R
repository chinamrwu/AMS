ggptROC <- function(roc, showAUC = TRUE, interval = 0.2, breaks = seq(0, 1, interval)){
  require(pROC)
  if(class(roc) != "roc")
    simpleError("Please provide roc object from pROC package")
  plotx <- rev(roc$specificities)
  ploty <- rev(roc$sensitivities)
  
  ggplot(NULL, aes(x = plotx, y = ploty)) +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) + 
    geom_step() +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) + 
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    theme_bw() + 
    theme(axis.ticks = element_line(color = "grey80")) +
    coord_equal() + 
    annotate("text", x = interval/2, y = interval/2, vjust = 0, label = paste("AUC =",sprintf("%.3f",roc$auc)))
}