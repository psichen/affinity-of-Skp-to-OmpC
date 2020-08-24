sample <- '0.35uM/'

c <- c(0.01, 0.14) #center
h <- c(6.36, 1.02) #height
w <- c(0.03, 0.12) #sigma 
o <- 0.1 #offset

x <- read.csv(paste0('./',sample,'x_hist.txt'), header=FALSE)$V1
y <- read.csv(paste0('./',sample,'y_hist.txt'), header=FALSE)$V1
y <- y - o
y[y<0] <- 0

f <- function(x, c, h, w) return(h*exp(-(x-c)^2/(2*w^2)))
xx <- seq(-.2, 1.18, .001)

peak.1 <- f(xx,c[1],h[1],w[1])
peak.2 <- f(xx,c[2],h[2],w[2])
peak <- peak.1 + peak.2


p <- ggplot(NULL)+
    geom_col(aes(x, y), fill='#DC0000FF')+
    geom_area(aes(xx, peak.1), size=1, fill='#B09C85FF', alpha=.8)+
    geom_area(aes(xx, peak.2), size=1, fill='#4DBBD5FF', alpha=.8)+
    geom_line(aes(xx, peak), size=1, col='#000000FF')+
    scale_x_continuous(breaks=seq(-.2,1.2,.2))+
    xlab('FRET efficiency')+
    ylab('Frequency')+
    theme_bw()
