rm(list=ls())
d <- data.frame(
    BinDr = read.csv('./BinDr.txt', header=FALSE)$V1,
    BinAc = -1*read.csv('./BinAc.txt', header=FALSE)$V1,
    SelDr = read.csv('./BinDr_sel.txt', header=FALSE)$V1,
    SelAc = -1*read.csv('./BinAc_sel.txt', header=FALSE)$V1
)

d <- cbind(melt(d), time=seq(length(d$BinDr)))

p <- ggplot(d)+
    geom_line(aes(time, value, col=variable))+
    labs(x='time / ms', y='fluorescence counts')+
    guides(col=guide_legend(title='trace'))+
    scale_color_npg()+
    theme_bw()