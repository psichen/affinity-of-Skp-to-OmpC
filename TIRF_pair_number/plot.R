data <- read.csv('./pair_num.csv', sep='\t')

f <- function(k, n, x) return(x^n/(x^n+k^n))

p <- ggplot(data)+
    geom_point(aes(conc, mean))+
    geom_errorbar(aes(x=conc, ymin=mean-err, ymax=mean+err), width=.1)+
    stat_function(fun=f, args=list(k=.54548,n=1.63226), xlim=c(0, 8), color='red')+
    scale_y_continuous(name='normalized proportion of bound-OmpC', breaks=seq(0,1,.1))+
    scale_x_continuous(name=expression('[Skp]'[0]~' / nM'), breaks=seq(0,8))+
    theme_bw()