y <- read.csv('./FRET_hist_1102.txt', header=FALSE)$V1
y <- y/sum(y)
x <- seq(-.19,1.19,.02)

c <- c(.34718, .55506)
w <- c(.20738, .23786)
A <- c(.00876, .01018)
o <- 7.54503e-4

f <- function(x, c, A, w) return(A/w/sqrt(pi/2)*exp(-2*((x-c)/w)^2))
xx <- seq(-.2, 1.18, .001)

peak.1 <- f(xx, c[1], A[1], w[1])
peak.2 <- f(xx, c[2], A[2], w[2])
peak <- peak.1 + peak.2

p <- ggplot(NULL)+
    geom_col(aes(x, y), fill='#DC0000FF')+
    geom_area(aes(xx, peak.1+o), size=1, fill='#B09C85FF', alpha=.8)+
    geom_area(aes(xx, peak.2+o), size=1, fill='#4DBBD5FF', alpha=.8)+
    geom_line(aes(xx, peak+o), size=1, col='#000000FF')+
    scale_x_continuous(breaks=seq(-.2,1.2,.2), limits=c(-.2, 1.2))+
    scale_y_continuous(limits=c(0,.05))+
    theme_bw()

ggsave('./E_x_hist.pdf', p, dpi=300, width=12.9, height=3)