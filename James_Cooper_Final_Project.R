library(tseries)
library(zoo)
library(fBasics)
library(lmtest)
library(tibble)
library(ggplot2)
library(xts)
library(astsa)
library(forecast)
library(devtools)
library(ggfortify)
library(repr)
install_github('sinhrks/ggfortify')

  

setwd('C:\\Users\\James Cooper\\Desktop\\DePaul\\Time Series Analysis\\Final Project')
emp_data = read.table("Employment.csv", header=T, sep=",")
head(emp_data)
anyDuplicated(emp_data$Month)

#################################################
#        Using zoo for time series              #
#################################################

#Turning the Construction data into a time series object using zoo
constr_time_series = zoo(emp_data$Construction, as.Date(as.character(emp_data$Month),
                                                     format = "%m/%d/%Y"))

head(constr_time_series)
plot(constr_time_series, main = "US Construction Employment by Year")
autolayer(constr_time_series, colour='Blue') + ggtitle('US Construction Employment by Year')
abline(reg=lm(constr_time_series~time(constr_time_series)))
help(autoplot.ts)


cdiff = diff(constr_time_series)
acf2(emp_data$Construction)

plot(diff(emp_data$Construction), type='l')
acf2(diff(emp_data$Construction))

sconstr = diff(cdiff)
plot(sconstr)

#################################################
#       Method 2 - using ts function            #
#################################################

constr_ts = ts(emp_data$Construction, start=c(2001,4), freq=12)
constr_ts
autoplot(constr_ts, xlab='Year', ylab='Number Employed (in thousands)',colour = 'orange', size=1.5, linetype = ) +
  ggtitle('US Construction Employment by Year') + 
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") 
emp_data[,1:2]

############################################# 
#           Exploratory Analysis            #
#############################################
head(constr_ts)
summary(emp_data$Construction)
start(constr_ts)
end(constr_ts)
frequency(constr_ts)
basicStats(constr_ts)
hist(constr_ts)
boxplot(constr_ts~cycle(constr_ts), xlab='Month', ylab='Total Employment (in thousands)', main='Seasonal Distribution')
normalTest(constr_ts,method=c('jb'))
adf.test(constr_ts)

#Histograms
qplot(emp_data$Construction, geom='histogram', xlab='Employment Totals', bins=12, fill=I("slategrey"), 
      col=I("orange"))
# add approximating normal density curve
xfit<-seq(min(emp_data$Construction),max(emp_data$Construction),length=500)
yfit<-dnorm(xfit,mean=mean(emp_data$Construction),sd=sd(emp_data$Construction))
lines(xfit, yfit, col="blue", lwd=2) 
# creates 2 by 2 display for 4 plots
hist(emp_data$Construction, xlab="Construction Employment", prob=TRUE, main="Histogram", col='slategrey')
# add approximating normal density curve
xfit<-seq(min(emp_data$Construction),max(emp_data$Construction),length=40)
yfit<-dnorm(xfit,mean=mean(emp_data$Construction),sd=sd(emp_data$Construction))
lines(xfit, yfit, col="blue", lwd=2) 

#QQ Plots
qqnorm(emp_data$Construction)
qqline(emp_data$Construction, col=2)

#Box Tests
Box.test(constr_ts,lag=2,type='Ljung')
adf.test(constr_ts, alternative = 'stationary')

#Plot with mean
plot(constr_ts, main = 'US Construction Monthly Employment')
abline(reg=lm(constr_ts~time(constr_ts)))

#General Trend
plot(aggregate(constr_ts, FUN=mean))
#plot(constr_ts, type="l", lwd=2, col="red", ylab= "# employed",xlim=c(2001,2018),axes=F)
#axis(1,at=2001:2018,labels=2001:2018);axis(2);box()

#Returns
rets = log(constr_ts/lag(constr_ts, -1))
ret = coredata(rets)
ret
basicStats(rets)
par(mfrow=c(3,1))
plot(rets, type='l')
plot(rets^2, type='l')
acf2(rets)
Box.test((rets^2),lag=2,type='Ljung')


#####################################
#     Differenced Data 1x and 2x    #
#####################################

#Plotting 1x differenced data, difference of log data, and 2x differenced data
#We can see that 2x differencing makes the data stationary
par(mfrow = c(1,1))
plot(constr_ts)
plot(diff(constr_ts))
plot(log(constr_ts))
plot(diff(log(constr_ts)))
plot(diff(diff(constr_ts)))
adf.test(diff(log(constr_ts), alternative='stationary'))
adf.test(diff(diff(log(constr_ts)), alternative='stationary'))

#Exploratory of 2x diff
diff = diff(constr_ts)
d_diff = diff(diff, 2)
plot(d_diff)
abline(reg=lm(d_diff~time(d_diff)))
summary(d_diff)
basicStats(d_diff)
acf2(d_diff)
qqnorm(d_diff)
qqline(d_diff, col=2)

#Histogram of the 2x differenced data "ddiff"
hist(ddiff, xlab="Construction Employment", main="Histogram")
xfit<-seq(min(ddiff),max(ddiff),length=40) 
yfit<-dnorm(xfit,mean=mean(ddiff),sd=sd(ddiff)) 
lines(xfit, yfit, col="blue", lwd=4) 

#Using log 1x difference and turning and using arima model with difference order = 1
log_ts = diff(log(constr_ts))
plot(log_ts)
hist(log_ts, xlab="Construction Employment", main="Histogram")
hist(diff(log_ts), xlab="Construction Employment", main="Histogram")
acf(diff(log_ts))
pacf(diff(log_ts))
normalTest(log_ts,method=c('jb'))

adf.test((log_ts), alternative='stationary')


#####################################
#         Testing Models            #
#####################################

auto.arima(log_ts, trace=T)

log_model = Arima(log_ts, c(0,1,2), method="ML")
log_model2 = Arima(log_ts, c(1,1,0), method="ML")
log_model3 = Arima(log_ts, c(0,1,1), method="ML")
log_seasonal = Arima(log_ts, order = c(0,1,1), seasonal = list(order=c(0,0,1), period = 12))
log_seasonal
coeftest(log_seasonal)

sarima(log_ts, 0,1,1)
sarima(log_ts, 2,1,0)

pm1 = backtest(log_model, log_ts, orig=199, h=1)
pm2 = backtest(log_model2, log_ts, orig=199, h=1)
pm3 = backtest(log_model3, log_ts, orig=199, h=1)
pm4 = backtest(log_seasonal, log_ts, orig=199, h=1)

#####################################
#          Forecasting              #
#####################################

f = forecast(log_model3, h=5)
f
plot(f, include=200)
lines(ts(c(f$fitted, f$mean), frequency=12,start=c(2001,4)), col="blue")
log_model3

sarima.for(constr_ts, n.ahead=12, 0,1,1,0,0,3,12)

d.arima <- auto.arima(log_ts, trace = T)
d.arima
d.forecast <- forecast(d.arima, level = c(95), h = 50)
autoplot(d.forecast)

coeftest(log_model)
m = forecast(log_model, h=2)
plot(m)

#Testing for Normality
qqnorm(diff(log_ts)) 
qqline(diff(log_ts), col = 2)

normalTest(d_diff,method=c("jb"))  

Box.test(d_diff, lag=1, type='Ljung')

#We can see at lag 2 we can reject the null hypothesis
Box.test(d_diff, lag=2, type='Ljung')


#Fitting a model
ddiff = diff(diff(constr_ts), lag=12)
acf2(ddiff)

acf2(d_diff, max.lag=60)

sarima(ddiff, 2, 1, 0, 1, 0, 0, 3)

m = sarima(constr_ts, 2, 1, 0, 1, 0, 0, 3)
backtest(m, log_ts)

#Forecast
sarima.for(constr_ts, n.ahead=12, 2, 1, 0, 1, 0, 0, 3)
plot(constr_ts)



#Testing
plot(acf2(constr_ts))

sarima(constr_ts, 3, 2, 3)

sarima(d_diff, 0, 0, 3, 1, 0, 0, 12)
pred = predict(s, n.ahead=10*12)


kpss.test(constr_ts)
adf.test(log_ts, alternative = 'stationary')
adf.test(diff(log_ts), alternative='stationary')


######################################## 
#           KAGGLE Function            #
########################################

ts.decomp <- function(df, col = 'Series Name', span = 0.13, Mult = TRUE)
{
  if(Mult) temp = log(df)  else temp = df
  spans = span * length(temp)  
  fit <- stl(temp, s.window = "periodic", t.window = spans)
  plot(fit, main = paste('Decompositon of',col,'with loess span = ', as.character(span)))
  fit$time.series
}

Constr.decomp <- ts.decomp(log(constr_ts), col = 'Employment', Mult = FALSE)
str(Constr.decomp)

#Plot the ACF & PACF of the remainder
options(repr.pmales.extlot.width=8, repr.plot.height=6)
acf(Constr.decomp[, 3], is.df = FALSE)
pacf(Constr.decomp[, 3], is.df = FALSE)




ts.model = function(ts, col = 'remainder', order = c(0,0,1))
{
  mod = arima(ts, order = order, include.mean = FALSE)
  print(mod)
  mod
}

arima.estimate1 <- ts.model(Constr.decomp[, 3], order = c(1,0,1))#ARIMA(1,0,1) model
arima.estimate2 <- ts.model(Constr.decomp[, 3], order = c(0,1,1))#ARIMA(0,1,1) model
arima.estimate3 <- ts.model(Constr.decomp[, 3], order = c(1,1,3))#ARIMA(1,1,3) model 
arima.estimate4 <- ts.model(Constr.decomp[, 3], order = c(0,1,2))#ARIMA(0,1,5) model
arima.estimate5 <- ts.model(Constr.decomp[, 3], order = c(2,1,3))#ARIMA(2,1,3) model
arima.estimate6 <- ts.model(Constr.decomp[, 3], order = c(3,0,4))#ARIMA(3,0,4) model
arima.estimate7 <- ts.model(Constr.decomp[, 3], order = c(4,0,5))#ARIMA(4,0,5) model
arima.estimate8 <- ts.model(Constr.decomp[, 3], order = c(4,1,5))#ARIMA(4,1,5) model

cat(paste('Sigma^2 of the original series = ', as.character(var(log(constr_ts)))))


acf(arima.estimate7$resid[-1], is.df = FALSE)
pacf(arima.estimate7$resid[-1], is.df = FALSE)


fit_constr <- auto.arima((constr_ts), max.p = 5, max.q = 5, max.P = 5,
                           max.Q = 5, max.order = 5, max.d = 0, max.D = 5, start.p = 0,
                           start.q = 0, start.P = 0, start.Q = 0)
summary(fit_constr)


constr_forecast <- forecast(fit_constr, h=30)
summary(constr_forecast)
plot(constr_forecast, lwd=3)
lines(ts(c(constr_forecast$fitted, constr_forecast$mean), frequency=12,start=c(2001,4)), col="blue", lwd=1)

model_1 = arima(log(constr_ts), c(0,1,1))
model_1_forecast = forecast(model_1, h=5)
model_1_forecast
plot(model_1_forecast, include=10)
lines(ts(c(f$fitted, f$mean), frequency=12,start=c(2001,4)), col="blue", lwd=3)
log_model3
