rm(list=ls())
#d <- read.csv('./area_data_fix_C2_position_width_6_unnormalized.csv')
#d <- read.csv('./area_data_fix_C1_position_6_unnormalized.csv')
d <- read.csv('./area_data_fix_position_6.csv')

#manual!!
#fit <- nlsLM(data=subset(d, conc>0.35), c.norm~conc/(kd+conc), start=c(kd=1), weights=wfct(1/c.norm.err^2))
fit <- nlsLM(data=subset(d, conc>0.35), c.norm~a*conc/(kd+conc), start=c(a=.1, kd=1))#, weights=wfct(1/c.norm.err^2))

p <- ggplot()+
    geom_point(data=d, aes(conc*3, c.norm))+
    geom_errorbar(data=d, aes(x=conc*3, ymax=c.norm+c.norm.err, ymin=c.norm-c.norm.err), width=.6)+
    #geom_line(data=NULL, aes(x=seq(0,13,.1), y=seq(0,13,.1)/(coef(fit)['kd']+seq(0,13,.1))), col='red')+
    geom_line(data=NULL, aes(x=seq(0,40,.1), y=coef(fit)['a']*seq(0,40,.1)/3/(coef(fit)['kd']+seq(0,40,.1)/3)), col='red')+
    #geom_line(data=NULL, aes(x=seq(0,13,.1), y=seq(0,13,.1)^coef(fit)['H']/(coef(fit)['kd']+seq(0,13,.1)^coef(fit)['H'])), col='red')+
    theme_bw()+
    xlab(expression('Skp'*' / '*mu*'M'))+
    ylab(expression('Normalized portion of '*OmpC(Skp[3])[2]))
