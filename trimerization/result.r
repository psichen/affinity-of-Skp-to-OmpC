d <- read.csv('./summary2.csv',sep=',')
d <- mutate(d, n=n/.85, n.err=n.err/.85, conc=conc*3169/1490*0.66, y0=y0*3169/1490*0.66) # correct for conc
d.fit <- subset(d, conc<200&conc>1)

#f.fit <- function(x,Kd,n1,n2,hill){
    #delta <- x^2*Kd^2/36+Kd^3/729
    #s <- (x*Kd/6+sqrt(delta))^(1/3)+abs(x*Kd/6-sqrt(delta))^(1/3)*sign(x*Kd/6-sqrt(delta))
    #f <- s^2/(s^2+Kd)
#    f <- x^hill/(x^hill+Kd)
#    return(n1*f+n2*(1-f))
#}

f.fit <- function(x,k,k1){
    # k1: C <--> D + 2S
    # k : C <--> 3S
    # x : n
    # xx: (n-1)/(3-n)
    # return: y0, conc of unlaebl monomer
    xx = (x-1)/(3-x)
    return(3*k1^(3/2)/k*xx^(3/2)+k1^(1/2)*xx^(1/2))
}

#n.fit<-nlsLM(n~f.fit(conc,Kd,n1,n2,hill),data=d,start=c(Kd=400,n1=3,n2=1,hill=2))
#n.fit<-nlsLM(n~f.fit(conc,Kd,n1,1,2),data=d,start=c(Kd=400,n1=3))

fit <- nlsLM(y0~f.fit(n,k,k1), data=d.fit, start=c(k=20000,k1=1000))

yy <- seq(1.0002, 2.96, .001)
xx <- f.fit(yy, coef(fit)['k'], coef(fit)['k1'])
s <- sqrt((yy-1)*coef(fit)['k']/(3-yy))
xx_k <- 3*s^3/coef(fit)['k']+s

#tD.fit<-nlsLM(tD~f.fit(conc,Kd,n1,n2,hill),data=subset(d,conc<1000),start=c(Kd=400,n1=710,n2=520,hill=2))
# 

dd <- read.csv('./summary2_plot.csv',sep=',')
dd <- mutate(dd, n=n/.85, n.err=n.err/.85, conc=conc*3169/1490*0.66) # correct for conc

p.n <- ggplot()+
    geom_line(data=NULL, aes(xx,yy), col='red')+
    geom_point(data=dd, aes(conc, n))+
    geom_errorbar(data=dd, aes(conc, ymin=n-n.err,ymax=n+n.err))+
    scale_x_log10(name=expression(Skp[1]*' / nM'), breaks=c(.1,1,10,100,1000,10000,100000), minor_breaks=as.numeric(1:10 %o% 10 ^ (-1:4)), labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(name='relative oligomer state', breaks=seq(1,3), limits=c(0.7, 3.3))+
    theme(legend.position='none')+
    theme_bw()

# 
# p.tD <- ggplot(d, aes(conc, tD))+
#     geom_point(aes(col=factor(lot)))+
#     geom_errorbar(aes(ymin=tD-tD.err,ymax=tD+tD.err,col=factor(lot)))+
#     scale_x_log10(breaks=c(.1,1,10,100,1000), minor_breaks=as.numeric(1:10 %o% 10 ^ (-1:4)))+
#     labs(x='skp subunits / nM',y=expression(diffusion~time~'/'~mu*s))+
#     theme(legend.position='none')
# 
# p.n.fit <- p.n+stat_function(fun=f.fit, args=coef(n.fit), xlim=c(-.5,3.1), col='red')
#p.n.fit <- p.n+stat_function(fun=f.fit, args=c(coef(n.fit),n2=1,hill=2), xlim=c(-.5,3.1), col='red')
#p.tD.fit <- p.tD+stat_function(fun=f.fit, args=coef(tD.fit), xlim=c(-.5,3.1), col='red')