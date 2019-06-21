drawViolin <- function(value,sample1,strTitle="Violin plot"){ 
  library(ggplot2)
    df1 <- data.frame('probeId'=sample1 , value = value )
    ggplot(df, aes(x = probeId, y=value,fill= probeId)) + 
    geom_violin(trim=FALSE)+
    labs(ylab="Value",xlab="Violin")+
    geom_boxplot(width=0.1)+
    theme(legend.direction = 'horizontal',legend.position = 'top',
	      panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"))
}

