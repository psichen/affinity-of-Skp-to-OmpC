rm(list=ls())

d.box <- data.frame()
d.cur <- data.frame()

source('./fit_OmpC_Skp_1.r')
d.box <- rbind(
           d.box,
           data.frame(
                      tD = tD.csv$a*1e6,
                      lab = 'OmpC_Skp_1'
                      )
           )
d.cur <- rbind(
               d.cur,
               fcs.coef
               )

source('./fit_OmpC_Skp_2.r')
d.box <- rbind(
           d.box,
           data.frame(
                      tD = tD.csv$a*1e6,
                      lab = 'OmpC_Skp_2'
                      )
           )
d.cur <- rbind(
               d.cur,
               fcs.coef
               )

p <- ggplot(NULL)+
    geom_point(data=subset(d.cur, i=='15_40'), aes(pt, tD*1e6, col=factor(rep)))+
    geom_point(data=subset(d.cur, i=='65_85'), aes(pt, tD*1e6, col=factor(rep)))+
    geom_line(data=subset(d.cur, i=='15_40'), aes(pt, fit*1e6, col=factor(rep)))+
    geom_line(data=subset(d.cur, i=='65_85'), aes(pt, fit*1e6, col=factor(rep)))+
    scale_x_continuous(name='peak threshold', breaks=seq(40,70,10))+
    scale_y_continuous(name=expression('apparent diffusion time /'~mu*s))+
    scale_color_npg()+
    theme_bw()
