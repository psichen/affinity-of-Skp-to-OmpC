sample <- './1.4uM/'
for(rep in 1:3){
x <- seq(-.19, 1.17, .02)
y <- read.table(paste0(sample,'y_hist.txt'))$V1
m <- sample(x, size=10000, replace=T, prob=y)
n <- hist(m, breaks=seq(-.2, 1.2, .02))

write.table(n$density, paste0(sample, 'bootstrap_', rep, '.txt'), quote=F, row.names=F, col.names=F)
}
