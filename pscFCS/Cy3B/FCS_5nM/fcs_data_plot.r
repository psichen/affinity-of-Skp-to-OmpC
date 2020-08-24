# parameters initialization
rm(list=ls())
path <- './Cy3B_5nM_'
fcs.model <- '2D1R'
fcs.time.start <- 1e-6
fcs.time.end <- 1e0

#fcs fitting equation
#2D1R function
fcs.2D1R <- function(x, n, tD, A, tR, offset){
    #n: molecule number
    #tD: diffusion time, us
    #T: triplet state ratio
    #tR: triplet life time, us
    y = offset+1+1/n/(1+x/tD*1e6)*(1+A*exp(-x/tR*1e6))
    return(y)
}

# file reading
read.sin <- function(path){
    fcs.skip <- grep('CorrelationFunction', readLines(path))[1]
    fcs.end <- grep('RawCorrelationFunction', readLines(path))

    d.fcs <- read.table(path, skip=fcs.skip, nrow=fcs.end-fcs.skip-2, col.names=c('time', 'AxB', 'BxA'))
    d.fcs <- filter(d.fcs, time>fcs.time.start&time<fcs.time.end)
    d.fcs <- mutate(d.fcs, G=.5*(AxB+BxA))
    return(d.fcs)
}

d.fcs <- data.frame()
for(i in seq(1)-1){
    sample <- paste0(path,i,'.sin')
    d.fcs <- rbind(d.fcs, read.sin(sample))
}

d.mean <- dcast(melt(d.fcs, id.vars='time', measure.vars='G'), time~variable, mean)

#fcs fitting
if(fcs.model == '2D1R'){
    fcs.value <- list(n=1, tD=500, A=.1, tR=10, offset=0)
    fcs.fit <- nlsLM(G~fcs.2D1R(time, n, tD, A, tR, offset), data=d.mean, start=fcs.value)
    fcs.fit.par <- list(
                        n=round(coef(fcs.fit)['n'], 1),
                        tD=round(coef(fcs.fit)['tD'], 0),
                        A=round(coef(fcs.fit)['A'], 2),
                        tR=round(coef(fcs.fit)['tR'], 1),
                        offset=round(coef(fcs.fit)['offset'], 3)
                         )
}

p <- ggplot()+
    geom_point(data=d.fcs, aes(time, G), size=1.5, alpha=.5, col='grey')+
    #geom_line(data=d.mean, aes(time, G))+
    geom_line(data=d.mean, aes(time, predict(fcs.fit)), col='red')+
    scale_x_log10(name='lag time / s', breaks=10^(-6:0), minor_breaks=as.numeric(1:10%o%10^(-6:0)), labels = trans_format("log10", math_format(10^.x)))+
    labs(y='correlation')+
    theme_bw()
