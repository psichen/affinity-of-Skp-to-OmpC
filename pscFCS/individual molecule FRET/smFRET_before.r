rm(list=ls())
x <- read.csv('./x_hist.txt', header=FALSE)$V1
a <- data.frame(y = read.csv('./fret_a.txt', header=FALSE)$V1, molecule='a')
b <- data.frame(y = read.csv('./fret_b.txt', header=FALSE)$V1, molecule='b')
c <- data.frame(y = read.csv('./fret_c.txt', header=FALSE)$V1, molecule='c')
d <- data.frame(y = read.csv('./fret_d.txt', header=FALSE)$V1, molecule='d')
e <- data.frame(y = read.csv('./fret_e.txt', header=FALSE)$V1, molecule='e')
f <- data.frame(y = read.csv('./fret_f.txt', header=FALSE)$V1, molecule='f')

d <- cbind(x, rbind(a, b, c, d, e, f))
d$y.n <- with(d, y/4394/.02)

p <- ggplot(d)+
    geom_col(aes(x,y.n,fill=molecule), col='black', position='stack')+
    ylim(c(0, 8.2))+
    labs(x='FRET efficiency', y='# events')+
    guides(fill=FALSE)+
    theme_bw()+
    scale_color_npg()