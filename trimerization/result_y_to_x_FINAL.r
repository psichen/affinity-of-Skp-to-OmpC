rm(list=ls())
d <- read.csv('./summary2_plot.csv',sep=',')
d <- mutate(d, n=n/.85, n.err=n.err/.85, conc=conc*3169/1490*0.66) # correct for conc
#d.fit <- subset(d, conc<200&conc>1)

f.fit <- function(x,k,k1,H,n1,n3){
    # k1: C1 <--> D + 2S
    # k : C0 <--> 3S
    # x : S+3C0
    # H : hill coefficient
    # n1 : theoretic = 1
    # n3 : theoretic = 3
    # return: n, effective oligomer number
    p = k/3
    q = x*k/(-3)
    delta = (q/2)^2+(p/3)^3
    S = sign(q/(-2)+sqrt(delta))*abs(q/(-2)+sqrt(delta))^(1/3)+sign(q/(-2)-sqrt(delta))*abs(q/(-2)-sqrt(delta))^(1/3) 
    n = n3*(S^H/(k1+S^H))+n1*(1-S^H/(k1+S^H))
    return(n)
}

#fit <- nlsLM(n~f.fit(conc,k,k1,H,n1,n3), data=d, start=c(k=20000,k1=1000,H=2,n1=1,n3=3), weights=wfct(1/n.err^2))
#================manual run!!, reason unknown
fit <- nlsLM(n~f.fit(conc,k,k1,H=2,n1=1,n3=3), data=d, start=c(k=40000,k1=1000), weights=wfct(1/n.err^2))

p.n <- ggplot(data=d, aes(conc, n))+
    geom_point()+
    geom_errorbar(aes(conc, ymin=n-n.err,ymax=n+n.err))+
    scale_x_log10(name=expression(Skp[1]*' / nM'), breaks=10^(-1:3), minor_breaks=as.numeric(1:10 %o% 10 ^ (-1:3)), labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(name='relative oligomer state', breaks=seq(3), limits=c(.7, 3.3))+
    theme(legend.position='none')+
    theme_bw()+
    #stat_function(fun=f.fit, args=list(k=coef(fit)['k'],k1=coef(fit)['k1'],H=coef(fit)['H'],n1=coef(fit)['n1'],n3=coef(fit)['n3']), col='red')
    stat_function(fun=f.fit, args=list(k=coef(fit)['k'],k1=coef(fit)['k1'],H=2,n1=1,n3=3), color='red')