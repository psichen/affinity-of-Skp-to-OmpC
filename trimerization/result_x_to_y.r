rm(list=ls())
d<- read.csv('./summary2_plot.csv',sep=',')
d <- mutate(d, n=n/.85, n.err=n.err/.85, conc=conc*3169/1490*0.66) # correct for conc
d.fit <- subset(d, conc<200&conc>1)

f.fit <- function(n,k,k1,n1,n3){
    # k1: C1 <--> D + 2S
    # k : C0 <--> 3S
    # n : n3*f+n1*(1-f)
    # n1 : theoretic = 1
    # n3 : theoretic = 3
    # return: x, total concentration of unlabeled Skp monomer
    #S = sqrt(k1*(n-n1)/(n3-n))
    #x = S+3/k*S^3
    S2 = k1*(n-n1)/(n3-n)
    x = 3/k*S2^(3/2)+S2^(1/2)
    return(x)
}

#fit <- nlsLM(conc^2~f.fit(n,k,k1,n1,n3), data=d.fit, start=c(k=20000,k1=1000,n1=1,n3=3))
#poorly fitted when open n1 and n3

fit <- nlsLM(conc~f.fit(n,k,k1,n1=1,n3=3), data=d.fit, start=c(k=70000,k1=1000))
yy <- seq(1.0003,2.96,.01)
xx <- f.fit(yy, k=coef(fit)['k'],k1=coef(fit)['k1'],n1=1,n3=3)

p.n <- ggplot()+
    geom_point(data=d, aes(n, conc))+
    geom_errorbarh(data=d, aes(y=conc, xmin=n-n.err,xmax=n+n.err))+
    geom_line(data=NULL, aes(yy, xx), col='red')+
    scale_y_log10(breaks=c(.1,1,10,100,1000,10000,100000), minor_breaks=as.numeric(1:10 %o% 10 ^ (-1:4)), labels = trans_format("log10", math_format(10^.x)))+
    labs(x='relative oligomer state', y='skp monomer / nM')+
    theme(legend.position='none')+
    xlim(c(0.7,3.3))+
    theme_bw()+
    coord_flip()
