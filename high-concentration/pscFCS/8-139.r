#20180615 - modify p.fcs

rm(list=ls())
FCS.2D <- function(offset, G0, tD, time) offset+G0/(1+time/tD)
FCS.2D1R <- function(offset, G0, tD, A, tR, time) offset+G0/(1+time/tD)*(1+A*exp(-time/tR))
tDcurve <- function(a, b, x) a+b*x^2

sam<-'./8-139/'
reps<-c(2,3,4,6)

interval <- c('40_80')
pt.seq <- seq(40,70)
accc <- 'tt'
model = '2D'
fcs.start.time = 1e-5
fcs.save.each = 0
p.tD.save.each = 0
p.tD.all.save = 0
tD.csv.write = 0

tD.csv <- data.frame()
fcs.coef <- data.frame()
p.fcs <- list()

for(rep in reps){
    data <- paste0(sam,rep)
    tDout <- data.frame()
    ptCurve <- data.frame()
    fcsCurve <- data.frame()
    fcs.predict <- data.frame()

    fcs <- sapply(paste0('./', data, '/', interval, '/', accc, '/p', pt.seq, '.txt'), read.table)
    names(fcs) <- pt.seq
    fcs <- melt(fcs)
    names(fcs) <- c('g', 'pt')
    fcs$time <- read.table(paste0('./', data, '/lag_time.txt'))$V1
    fcs <- transform(fcs, pt=factor(pt, levels=pt.seq))
    fcs <- filter(fcs, time>fcs.start.time)
    #----------fit FCS curves----------
    fcs.temp <- data.frame() #store G0, tD series

    for(pt.val in pt.seq){
        if(model == '2D'){
            fcs.fit <- with(filter(fcs, pt==pt.val), nlsLM(g~FCS.2D(offset, G0, tD, time), start=c(offset=1, G0=1, tD=1e-3)))
        }
        if(model == '2D1R'){
            fcs.fit <- with(filter(fcs, pt==pt.val), nlsLM(g~FCS.2D1R(offset, G0, tD, A, tR, time), start=c(offset=1, G0=1, tD=1e-3, A=.1, tR=1e-5)))
        }
        fcs.temp <- rbind(fcs.temp, data.frame(as.list(coef(fcs.fit)), pt=pt.val))

        if(!exists('fcs.predict')){
            fcs.predict <- data.frame(predict=predict(fcs.fit))
        }else{
            fcs.predict <- rbind(fcs.predict, data.frame(predict=predict(fcs.fit)))
        }
    }
    fcs <- cbind(fcs, fcs.predict)

    tD.fit <- nls(tD~tDcurve(a, b, pt), data=fcs.temp, start=list(a=1e-3, b=1e-6), algorithm='port', lower=list(b=0))
    fcs.temp$fit <- predict(tD.fit)
    fcs.temp$sam <- sam
    fcs.temp$rep <- rep
    fcs.temp$i <- interval
    fcs.coef <- rbind(fcs.coef, fcs.temp)
    if(p.tD.save.each){
        p.tD <- ggplot(fcs.temp)+geom_point(aes(pt, tD))+geom_line(aes(pt, fit))+labs(x='peak threshold', y='diffusion time / s')
        ggsave(paste0('./', sam, '_', rep,  '_', interval, '.png'), p.tD)
    }

    tD.temp <- as.data.frame(t(coef(tD.fit)))
    tD.temp$sam <- sam
    tD.temp$rep <- rep
    tD.csv <- rbind(tD.csv, tD.temp)
    
    if(!exists('fcs.all')){
        fcs$rep <- rep
        fcs.all <- fcs
    }
    else{
        fcs$rep <- rep
        fcs.all <- rbind(fcs.all, fcs)
    }

    p.fcs[[rep]] <- ggplot(fcs)+geom_line(aes(time, g))+geom_line(aes(time, predict), color='red')+facet_wrap(~pt)+scale_x_log10()+labs(x='lag time / s', y='correlation')
}

p <- ggplot(fcs.coef)+
    geom_point(aes(pt, tD*1e6, col=factor(rep), shape=factor(rep)), size=3)+
    geom_line(aes(pt, fit*1e6, col=factor(rep)))+
    scale_x_continuous(name='peak threshold', breaks=seq(40,90,10))+
    scale_y_continuous(name=expression('apparent diffusion time /'~mu*s), limits=c(450, 960))+
    scale_color_npg()+
    theme_bw()
