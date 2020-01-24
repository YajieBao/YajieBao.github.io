setwd('F:\\STAT\\Time Series')
load("F:/STAT/Time Series/.RData")
library(forecast)
library(tseries)
candy=read.csv('candy_production.csv')
colnames(candy)=c('date','production')
production=candy$production
acf(production)
ts.plot(production)
adf.test(production)
desstat=function(data){
  library(moments)
  library(tseries)
  adft=adf.test(data)
  bt=Box.test(data)
  a=c(min(data),median(data),max(data),mean(data),
      sd(data),skewness(data),kurtosis(data),
      jarque.bera.test(data)$statistic,jarque.bera.test(data)$p.value,
      as.numeric(adft$statistic),adft$p.value,
      as.numeric(bt$statistic),bt$p.value)
  names(a)=c('Min','Median','Max','Mean','Sd.','Skewness',
             'Kurtosis','JB Stat.','p-value(JB)','DF Stat.','p-value(DF)',
             'Box Stat.','p-value(Box)')
  return(a)
}
x=seq(40,140,by=0.01)
xn=dnorm(x,mean(production),sd(production))
plot(density(production),main='Kernel Density Plot',lwd=2,ylim=c(0,0.025))
lines(x,xn,col='blue3')


prodiff=diff(production,lag=12)
diff.arima=auto.arima(prodiff)
summary(diff.arima)
acf(diff.arima$residuals,main='ACF Plot of Residuals')
pacf(diff.arima$residuals,main='PACF Plot of Residuals')
diff.res=diff.arima$residuals
qqnorm(diff.res)
qqline(diff.res)
Box.test(diff.res)
Box.test(diff.res^2,lag = 1)
library(rugarch)
garch.spec=ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=F),
                      variance.model = list(model=('sGARCH'),garchOrder=c(1,1)),
                      distribution.model ='std')
garch.spec1=ugarchspec(mean.model=list(armaOrder=c(2,0),include.mean=F),
                       variance.model = list(model=('sGARCH'),garchOrder=c(1,1)),
                       distribution.model ='norm')

diff.garch=ugarchfit(garch.spec,data = prodiff)
diff.garch
diff.garch1=ugarchfit(garch.spec1,data=prodiff)
plot(diff.garch1,which=9)
garchfor=ugarchforecast(diff.garch,n.ahead = 12)
plot(garchfor)
garchfore=garchfor@forecast[["seriesFor"]]+production[537:548]
garchupper=garchfore-garchfor@forecast[["sigmaFor"]]
garchlower=garchfore+garchfor@forecast[["sigmaFor"]]
op=par(cex.main=0.9,cex.lab=0.9,cex.axis=0.8,mar=c(5,4,4,2)-0.3)
ts.plot(production[500:548],lwd=2,xlab='Time',ylab='Production',
        xlim=c(1,62),main='Forecast Series \nwith unconditional 1-sigma bonds')
lines(c(49:60),garchfore,col='red3',lwd=2)
lines(c(49:60),garchupper,col='blue3',lwd=1)
lines(c(49:60),garchlower,col='blue3',lwd=1)
grid(nx=NULL,ny=NULL)
par(op)
plot(diff.garch)


n=length(production)
trainprodiff=prodiff[1:524]
testpro=production[537:548]
garchtrain=ugarchfit(garch.spec,data = trainprodiff)
testfore=ugarchforecast(garchtrain,n.ahead = 12)
garchtf=as.numeric(testfore@forecast[["seriesFor"]])
x12=production[525:536]
garchtf=garchtf+x12  ###还原x
sqrt(sum((garchtf-testpro)^2)/12)
## rolling sample
garchroll=NULL
for(i in 1:36){
  sam=prodiff[i:(i+499)]
  gar=ugarchfit(garch.spec,data = sam)
  fore=ugarchforecast(gar,n.ahead = 1)
  garchroll[i]=fore@forecast[["seriesFor"]]+production[(500+i)]
}
grolle1=sqrt(sum((garchroll-production[537:548])^2)/36)

garchroll12=matrix(0,ncol = 3,nrow=12)
for(i in 1:3){
  sam=prodiff[((i-1)*12+1):((i-1)*12+500)]
  gar=ugarchfit(garch.spec,data = sam)
  fore1=ugarchforecast(gar,n.ahead = 12)
  garchroll12[,i]=fore1@forecast[["seriesFor"]]+production[(500+12*(i-1)+1):(512+12*(i-1))]
}
grolle2=sqrt(sum((garchroll12-production[537:548])^2)/36)

prots=ts(production,frequency = 12,start = c(1972,1))


## season index
seasonindex=NULL
totalmean=mean(production)
for(i in 1:12){
  index=seq(i,548,12)
  seasonindex[i]=mean(production[index])/totalmean
}
op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8,mar=c(5,4,4,2)-0.3)
plot(seasonindex,type='o',main='Seasonal Index',xaxt='n',
     ylab = 'Production',xlab = '',col='red3',pch=15)
axis(1,at=c(1:12),labels = month.abb)


proholt=HoltWinters(prots,seasonal = "multiplicative")
plot(proholt)
grid(nx=NULL,ny=NULL)
holtfit=proholt$fitted[,1]
sqrt(sum((production[-c(1:12)]-holtfit)^2)/536)#98.37636
proholt$coefficients

holtres=production[-c(1:12)]-holtfit
ts.plot(holtres)
Box.test(holtres,lag = 1)
holtfore=forecast(proholt,h=12)
plot(holtfore)
grid(nx=NULL,ny=NULL)


library(Rssa)
candy.ssa=ssa(production,L=274)
plot(candy.ssa)
w <- wcor(candy.ssa, groups = 1:20)
plot(w, grid = c(2,4, 5,7))
r <- reconstruct(candy.ssa, groups = list(Trend = c(1, 4:5),
                                          Seasonality = c(2:3, 6:11)))

ssares=production-r$Trend-r$Seasonality
sqrt(sum(ssares^2)/548)
ts.plot(ssares)
plot(r,type='cumsum')

# error convergence curve
in.error=NULL
out.error=NULL
for(i in 8:230){
  r1=reconstruct(candy.ssa, groups = list(Trend = c(1, 4),
                                          Seasonality = c(2:3, 5:7,8:i)))
  in.error[i-7]=sqrt(sum((production-r1$Trend-r1$Seasonality)^2)/548)
}
op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8)
plot(c(8:230),in.error,type='l',col='blue3',xlab='Number of Components',
     ylab='RMSE',main='RMSE Curve of SSA',lwd=2)
grid(nx=NULL,ny=NULL)
par(op)

n=length(production)
trainpro=production[1:536]
testpro=production[537:548]
trainssa=ssa(trainpro)
set.seed(1012)
f1 <- bforecast(trainssa, groups = list(1:11), len = 12, R = 500)
sqrt(sum((f1[,1]-testpro)^2)/12)
traints=ts(trainpro,frequency = 12)
h1=HoltWinters(traints,seasonal = "multiplicative")
h11=as.numeric(forecast(h1,h=12)$mean)
sqrt(sum((h11-testpro)^2)/12)
# in-sample forcast comparison
op=par(cex.main=0.9,cex.lab=0.9,cex.axis=0.8)
ts.plot(production[537:548],xlim=c(1,13),ylab='Production',ylim=c(95,135),
        main='12-Steps Ahead Forecast of 3 Models')
lines(c(1:12),f1[,1],col='red3')
lines(c(1:12),h11,col='blue3')
lines(c(1:12),garchtf,col='magenta3')
grid(nx=NULL,ny=NULL,col='orange')
legend('topright',c('Actual','SSA','Holt-Winters','AR-GARCH'),cex=0.7,
       col=c('black','red3','blue3','magenta3'),lty=rep(1,4),inset = 0.01)

par(op)

for(i in 8:230){
  f <- bforecast(trainssa, groups = list(1:i), len = 12, R = 50)
  out.error[i-7]=sqrt(sum((f[,1]-testpro)^2)/12)
}

write.csv(out.error,'F:\\STAT\\Time Series\\out.error.csv')
out.error=read.csv('F:\\STAT\\Time Series\\out.error.csv')$x
oe=NULL
for(i in 4:7){
  f <- bforecast(trainssa, groups = list(1:i), len = 12, R = 50)
  oe[i-3]=sqrt(sum((f[,1]-testpro)^2))
}
out.error=c(oe,out.error)
out.error=sqrt(out.error^2/12)

op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8,mar=c(5,4,4,2))

plot(c(4:77),out.error[1:74],type='l',col='blue3',lwd=2,
     xlab = 'Number of Components',ylab='RMSE',main='Forecast RMSE Curve')
#points(c(4:77),out.error[1:74],col='blue3')
text(12,out.error[8],'11',cex=0.8,col='red3')
grid(nx=NULL,ny=NULL)

f <- bforecast(candy.ssa, groups = list(1:11), len = 1, R = 50)


htmean=holtfore$mean
htlower=holtfore$lower[,2]
htupper=holtfore$upper[,2]
op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8,mar=c(5,4,4,2)-0.3)
ts.plot(production[(n-35):n],xlim=c(1,52),ylim=c(80,150),lwd=2,
        ylab='Production',main='Forecasts and CI')
lines(c(37:48),f[,1],col='blue3',lwd=1)
lines(c(37:48),f[,2],col='blue3',lwd=1)
lines(c(37:48),f[,3],col='blue3',lwd=1)
lines(c(37:48),htmean,col='red3',lwd=1)
lines(c(37:48),htlower,col='red3',lwd=1)
lines(c(37:48),htupper,col='red3',lwd=1)
legend('topleft',c('Holt-Winters','SSA'),col=c('red3','blue3'),
       lty=c(1,1),inset = 0.01,cex=0.8,bg='grey90')
grid(nx=NULL,ny=NULL)


## 1-step forecast
ssaf=matrix(0,ncol=12,nrow = 48)
for(j in 1:12)
{
  for(i in 1:48){
  sam=production[i:(500+i-1)]
  ssam=ssa(sam)
  ssaf[i,j]=bforecast(ssam, groups = list(1:(2+j)), len = 1, R = 100)[1]
}
}
write.csv(ssaf,'F:\\STAT\\Time Series\\ssaf.csv')
ssafe1=NULL

for(i in 1:12){
  ssafe1[i]=sqrt(sum((production[513:548]-ssaf[,i][13:48])^2)/36)
}


htf=NULL
for(i in 1:48){
  sam=production[i:(500+i-1)]
  sam.ts=ts(sam,frequency = 12)
  htm=HoltWinters(sam.ts,seasonal = "multiplicative")
  htf[i]=forecast(htm,h=1)$mean
}
e3=sqrt(sum((production[513:548]-htf[13:48])^2)/36)
e3

op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8,mar=c(5,4,4,2)-0.3)
plot(c(3:14),ssafe1,ylim = c(1,31),col='red3',
     type='o',pch=15,xlab='Number of Components',ylab='RMSE',
     main='Rolling-sample 1-Step Forecast RMSE ')
abline(h=e3,col='blue3',lty=2)
abline(h=grolle1,col='magenta3',lty=2)
grid(nx=NULL,ny=NULL)
legend('topleft',c('SSA','Holt-Winters','AR-GARCH'),col=c('red3','blue3','magenta3'),
       lty=c(1,2,2),inset = 0.01,cex=0.7)
text(4,4,'5.5531',col='blue3',cex=0.8)
text(9,4,'4.9402',col='red3',cex=0.8)
text(4,8,'6.4905',col='magenta3',cex=0.8)
par(op)


## 12-step forecast
set.seed(1015)
ssaf12=matrix(0,ncol=12,nrow = 36)
for(j in 1:12)
  {
    for(i in 1:3){
      sam=production[(12*(i-1)+13):(512+(i-1)*12)]
      ssam=ssa(sam)
       ssaf12[((i-1)*12+1):(i*12),j]=bforecast(ssam,groups = list(1:(2+j)),
                                                     len = 12, R = 100)[,1]
       }
}

ssafe=NULL
for(i in 1:12){
  ssafe[i]=sqrt(sum((production[513:548]-ssaf12[,i])^2)/36)
}

htf12=NULL
for(i in 1:5){
  sam=production[(12*(i-1)+1):(500+(i-1)*12)]
  sam.ts=ts(sam,frequency = 12)
  htm=HoltWinters(sam.ts,seasonal = "multiplicative")
  htf12[((i-1)*12+1):(i*12)]=as.numeric(forecast(htm,h=12)$mean)
}

e2= sqrt(sum((production[513:548]-htf12[13:48])^2)/36)
e2

op=par(cex.main=0.9,cex.lab=0.8,cex.axis=0.8,mar=c(5,4,4,2)-0.3)
plot(c(3:14),ssafe,ylim = c(1,35),col='red3',
     type='o',pch=15,xlab='Number of Components',ylab='RMSE',
     main='Rolling-sample 12-Step Forecast RMSE ')
abline(h=e2,col='blue3',lty=2)
abline(h=grolle2,col='magenta3',lty=2)
grid(nx=NULL,ny=NULL)
legend('topleft',c('SSA','Holt-Winters','AR-GARCH'),col=c('red3','blue3','magenta3'),
       lty=c(1,2,2),inset = 0.01,cex=0.7)
text(4,11.5,'10.4605',col='blue3',cex=0.8)
text(9,5.5,'5.0580',col='red3',cex=0.8)
text(4,7,'8.9363',col='magenta3',cex=0.8)
par(op)

sediff= (production[537:548]-htf12[13:48])^2-(production[537:548]-ssaf12[,6])^2
garchf12=c(garchroll12[,1],garchroll12[,2],garchroll12[,3])
sediff2=(production[537:548]-garchf12)^2-(production[537:548]-ssaf12[,6])^2
ts.plot(sediff)
hist(sediff,breaks = 15)
plot(density(sediff))
elm=lm(sediff~1)
summary(elm)

###  把一个整数x转化成k位 2进制的数  ##########
binary<-function(x, k)
{ tmp=NULL
y=x
if(x<2^k) 
{ for(i in k-1:k)
{ a = floor(y/2^i)
#     print(c(i, a))
tmp=c(tmp, a)
y = y-a*2^i
}
}
2*(tmp-0.5)
}

######  计算一组置换样本对应的统计量  ##############
ppmr<-function(x)    ###  paired permutation
{ ind=rbinom(length(x), 1, 0.5)
ind=(ind-0.5)*2
mean(x*ind) 
}

x1=sediff
x1=rank(abs(x1))*sign(x1)
Obs1=mean(x1)

x2=sediff2
x2=rank(abs(x2))*sign(x2)
Obs2=mean(x2)
####  随机置换检验     ############################
#12-step Wilcoxon test
rep=10000
set.seed(1205)
results1 <- replicate(rep, ppmr(x1))
p.value1 <- length(results1[results1>=Obs1])/rep
print(p.value1)

results2 <- replicate(rep, ppmr(x2))
p.value2 <- length(results2[results2>=Obs2])/rep
print(p.value2)

#1step Wilcoxon test
y1=(production[513:548]-htf[13:48])^2-(production[537:548]-ssaf[,6][13:48])^2
ts.plot(y1)
yr1=rank(abs(y1))*sign(y1)
Obsy1=mean(y1)
set.seed(1026)
resultsy1 <- replicate(rep, ppmr(yr1))
p.valuey1 <- length(resultsy1[resultsy1>=Obsy1])/rep
print(p.valuey1)

y2=(production[513:548]-garchroll)^2-(production[537:548]-ssaf[,6][13:48])^2
ts.plot(y2)
yr2=rank(abs(y2))*sign(y2)
Obsy2=mean(y2)
set.seed(1026)
resultsy2 <- replicate(rep, ppmr(yr2))
p.valuey2 <- length(resultsy2[resultsy2>=Obsy2])/rep
print(p.valuey2)


pro500=production[1:536]
pro500ts=ts(pro500,frequency = 12)
htm=HoltWinters(pro500ts,seasonal = "multiplicative")
htf=forecast(htm,h=12)$mean
sqrt(sum((production[537:548]-htf)^2)/36)

