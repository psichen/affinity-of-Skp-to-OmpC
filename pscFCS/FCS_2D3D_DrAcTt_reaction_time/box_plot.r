rm(list=ls())
d.td <- rbind(
              cbind(
                    subset(read.csv('./data_2d.csv'), select = -g0),
                    model = '2D'
                    ),
              cbind(
                    subset(read.csv('./data_3d.csv'), select = -c(g0, w2)),
                    model = '3D'
                    )
              )

d.ratio <- rbind(
    data.frame(
        ratio = subset(d.td, model=='3D'&fret=='0.5_0.7'&chl=='Tt')$td/subset(d.td, model=='3D'&fret=='0.2_0.4'&chl=='Tt')$td,
        time.react = rep(c('0us','100us','200us','500us'), each = 20),
        model = '3D'
    ),
    data.frame(
        ratio = subset(d.td, model=='2D'&fret=='0.5_0.7'&chl=='Tt')$td/subset(d.td, model=='2D'&fret=='0.2_0.4'&chl=='Tt')$td,
        time.react = rep(c('0us','100us','200us','500us'), each = 20),
        model = '2D'
    )
)

p.fret0.3 <- ggplot(subset(d.td, fret==c('0.2_0.4')))+
    geom_boxplot(aes(chl, td*1e3, col=model, fill=time.react))+
    scale_color_npg()+
    theme_bw()+
    xlab('autocorrelations of channels')+
    ylab('diffusion time / ms')

p.fret0.6 <- ggplot(subset(d.td, fret==c('0.5_0.7')))+
    geom_boxplot(aes(chl, td*1e3, col=model, fill=time.react))+
    scale_color_npg()+
    theme_bw()+
    xlab('autocorrelations of channels')+
    ylab('diffusion time / ms')

p.ratio <- ggplot(d.ratio)+
    geom_boxplot(aes(time.react, ratio, col=model))+
    xlab('reaction time')+
    ylab('diffusion time ratio')+
    ylim(c(.375,.525))+
    scale_color_npg()+
    theme_bw()