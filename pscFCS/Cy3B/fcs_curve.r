d <- read.csv('./data.csv')

p <- ggplot(d)+
    geom_line(aes(time, norm, col=factor(pt)), size=1)+
    #geom_line(aes(time, norm.fit, col=factor(pt)), size=1)+
    scale_x_log10(name='lag time / s', breaks=10^(-6:0), minor_breaks=as.numeric(1:10%o%10^(-6:0)), labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(name='normalized correlation', limits=c(0, 1.15), breaks=seq(0,1,.25))+
    scale_color_npg()+
    theme_bw()