library(ggplot2)
p1=ggplot(tmp, aes(x=factor(x), y=y)) + geom_bar(stat='identity') 