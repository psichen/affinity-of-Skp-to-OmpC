k1 <- 1170
k <- 45470

y.u <- seq(.001, .995, .001)
x <- sqrt(k*y.u/3/(1-y.u)^3)
y.l <- x^2*(1-y.u)^2/(k1+x^2*(1-y.u)^2)

x.u.half <- sqrt(k*.5/3/(1-.5)^3)
x.l.half <- (3*k1^(3/2)+k*k1^(1/2))/k

p <- ggplot()+
geom_line(data=NULL, aes(x, y.u))+
#geom_line(data=NULL, aes(x, y.u, col='skp unlabel'))+
#geom_line(data=NULL, aes(x, y.l, col='skp label'))+
scale_x_log10(breaks=c(.1,1,1e1,1e2,1e3,1e4,1e5,1e6), minor_breaks=as.numeric(1:10 %o% 10 ^ (-1:6)), limits=c(1, 1e6), labels = trans_format("log10", math_format(10^.x)))+
labs(x=expression(Total~Skp~'/'~nM), y=expression(Fraction~of~Skp[3]))+
theme_bw()
