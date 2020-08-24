k1=6.7 #nM^-1s^-1
k2=6.6 #nM^-1s^-1
kD=4.55e4 #nM^2

skp = 10^seq(-3,3,.01)
x = 3/kD*skp^3+skp
skp3 = (x-skp)/3

p <- ggplot(NULL)+
    geom_line(aes(x, k1*skp, col='skp'))+
    geom_line(aes(x, k2*skp3, col='skp3'))+
    scale_x_log10(breaks=10^seq(-3,5), minor_breaks=as.numeric(1:10 %o% 10 ^ (-3:5)), labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks=10^seq(-13,6,2), labels = trans_format("log10", math_format(10^.x)))+
    labs(x='skp / nM', y='reaction rate / s^-1')+
    theme_bw()