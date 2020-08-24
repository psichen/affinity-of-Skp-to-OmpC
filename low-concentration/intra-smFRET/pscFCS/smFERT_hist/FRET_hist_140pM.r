sample <- ''

x <- read.csv(paste0('./',sample,'x_hist.txt'), header=FALSE)$V1
y <- read.csv(paste0('./',sample,'y_hist.txt'), header=FALSE)$V1
y <- y/sum(y)/.02

p <- ggplot(NULL)+
    geom_col(aes(x, y), fill='#DC0000FF')+
    geom_col(aes(x[x>.15&x<.4], y[x>.15&x<.4]), fill='#3C5488B2')+
    geom_col(aes(x[x>.65&x<.85], y[x>.65&x<.85]), fill='#F39B7FB2')+
    scale_x_continuous(name='FRET efficiency', breaks=seq(-.2,1.2,.2), limits=c(-.2, 1.2))+
    scale_y_continuous(name='Frequency', limits=c(0,6.5))+
    theme_bw()
