c <- c(0.01, 0.32) #center
h <- c(5.99, 0.69) #height
w <- c(0.04, 0.26) #sigma 
o <- 0.0 #offset

x <- read.csv(paste0('./x_hist.txt'), header=FALSE)$V1
y <- read.csv(paste0('./OS6.txt'), header=FALSE)$V1

y <- y/sum(y)/0.02

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
