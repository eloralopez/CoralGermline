#Exercise Set 3-1
xvals<-anscombe$x1
meanx<-mean(xvals)

yvals<-anscombe$y1
meany<-mean(yvals)

bslope<-sum(((xvals-meanx)*(yvals-meany)))/(sum((xvals-meanx)^2))    
bslope
aintercept<- meany - (bslope*meanx) 
aintercept
mod.fit <- lm(anscombe$y1 ~ anscombe$x1)
summary(mod.fit)

