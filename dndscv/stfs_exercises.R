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

#Exercise Set 5-1
#2
histfunction<- function(samplesize) { 
  samp.size<-samplesize
  n.samps<- 1000
  samps<-rnorm(samp.size * n.samps, mean=0, sd = 1)
  samp.mat<-matrix(samps, ncol=n.samps)
  samp.means<-colMeans(samp.mat)
  return(hist(samp.means, xlim=c(-5,5),ylim=c(0,250)))
}
par(mfrow=c(2,3))

s1<-histfunction(1)
s5<-histfunction(5)
s20<-histfunction(20)
s100<-histfunction(100)
s1000<-histfunction(1000)

histfunctionEXP<- function(samplesize) { 
  samp.size<-samplesize
  n.samps<- 1000
  samps<-rexp(samp.size * n.samps, rate=1)
  samp.mat<-matrix(samps, ncol=n.samps)
  samp.means<-colMeans(samp.mat)
  return(hist(samp.means, xlim=c(0,5),ylim=c(0,500)))
}
par(mfrow=c(2,3))

s1<-histfunctionEXP(1)
s5<-histfunctionEXP(5)
s20<-histfunctionEXP(20)
s100<-histfunctionEXP(100)
s1000<-histfunctionEXP(1000)
