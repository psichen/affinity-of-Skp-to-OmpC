rm(list=ls())
x <- read.csv('./x_hist.txt', header=FALSE)$V1
a <- data.frame(y = read.csv('./sel_fret_a.txt', header=FALSE)$V1, molecule='a')
b <- data.frame(y = read.csv('./sel_fret_b.txt', header=FALSE)$V1, molecule='b')

d <- cbind(x, rbind(a, b))
d$y.n <- with(d, y/4394/.02)

p <- ggplot(d)+
    geom_col(aes(x,y.n,fill=molecule), col='black', position='stack')+
    ylim(c(0, 8.2))+
    labs(x='FRET efficiency', y='# events')+
    guides(fill=FALSE)+
    theme_bw()+
    scale_color_npg()