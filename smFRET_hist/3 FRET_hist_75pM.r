sample <- '75pM/'

c <- c(0.01, 0.14, 0.73) #center
h <- c(5.67, 0.71, 0.52) #height
w <- c(0.04, 0.14, 0.22) #sigma 
o <- -0.08 #offset

x <- read.csv(paste0('./',sample,'x_hist.txt'), header=FALSE)$V1
y <- read.csv(paste0('./',sample,'y_hist.txt'), header=FALSE)$V1
y <- y - o
y[y<0] <- 0

f <- function(x, c, h, w) return(h*exp(-(x-c)^2/(2*w^2)))
xx <- seq(-.2, 1.18, .001)

peak.1 <- f(xx,c[1],h[1],w[1])
peak.2 <- f(xx,c[2],h[2],w[2])
peak.3 <- f(xx,c[3],h[3],w[3])
peak <- peak.1 + peak.2 + peak.3


p <- ggplot(NULL)+
    geom_col(aes(x, y), fill='#DC0000FF')+
    geom_area(aes(xx, peak.1), size=1, fill='#B09C85FF', alpha=.8)+
    geom_area(aes(xx, peak.2), size=1, fill='#4DBBD5FF', alpha=.8)+
    geom_area(aes(xx, peak.3), size=1, fill='#F39B7FFF', alpha=.8)+
    geom_line(aes(xx, peak), size=1, col='#000000FF')+
    scale_x_continuous(name='FRET efficiency', breaks=seq(-.2,1.2,.2), limits=c(-.2, 1.2))+
    scale_y_continuous(name='Frequency', limits=c(0,9))+
    theme_bw()
