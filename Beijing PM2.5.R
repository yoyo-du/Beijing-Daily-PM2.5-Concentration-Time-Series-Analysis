library(TSA)
library(corrplot)
library(dplyr)
library(ggplot2)
library(forecast)
library(tseries)
library(readxl)
library(lubridate)
library(urca)
library(forecast)
library(stats)

data <- read.csv("PRSA_Data_Wanliu_20130301-20170228.csv", header = TRUE)
data = as.data.frame(data)
str(data)
summary(data)

## check missing values
sum(is.na(data$PM2.5))/35064
sum(is.na(data$PRES))/35064
sum(is.na(data$TEMP))/35064
sum(is.na(data$WSPM))/35064
sum(is.na(data$DEWP))/35064

attach(data)

## check the distribution of the original data
par(mfrow=c(1,2))
boxplot(PM2.5)
boxplot(log(PM2.5))

par(mfrow=c(2,3))
boxplot(WSPM)
boxplot(log(WSPM))
boxplot(DEWP)
boxplot(TEMP)
boxplot(PRES)
boxplot(RAIN)


par(mfrow=c(2,3))
hist(PM2.5)
hist(WSPM)
hist(DEWP)
hist(TEMP)
hist(PRES)

### fill missing value

PM2.5.filled <- ifelse(is.na(PM2.5), median(na.omit(PM2.5)), PM2.5)
TEMP.filled <- ifelse(is.na(TEMP), median(na.omit(TEMP)), TEMP)
WSPM.filled <- ifelse(is.na(WSPM), median(na.omit(WSPM)), WSPM)
PRES.filled <- ifelse(is.na(PRES), median(na.omit(PRES)), PRES)
DEWP.filled <- ifelse(is.na(DEWP), median(na.omit(DEWP)), DEWP)

PM2.5 <- PM2.5.filled
TEMP <- TEMP.filled
WSPM <- WSPM.filled
PRES <- PRES.filled
DEWP <- PRES.filled


### daily aggregate
data$date <- paste(year, month, day, hour, sep = "/")
data = transform(data, date = as.POSIXct(as.character(date), format = "%Y/%m/%d/%H"))
data = as.data.frame(data)
attach(data)

logPM2.5 <- log(PM2.5)
ag1 <- as.data.frame(aggregate(logPM2.5 ~ floor_date(data$date, unit = "day"), data = data, mean))
logdaily_PM2.5 <- ts(ag1$logPM2.5, start = decimal_date(min(ag1$`floor_date(data$date, unit = "day")`)), frequency = 365)

logWSPM <- log(WSPM+0.000001)
ag1 <- as.data.frame(aggregate(logWSPM ~ floor_date(data$date, unit = "day"), data = data, mean))
logdaily_WSPM <- ts(ag1$logWSPM, start = decimal_date(min(ag1$`floor_date(data$date, unit = "day")`)), frequency = 365)

ag1 <- as.data.frame(aggregate(DEWP ~ floor_date(data$date, unit = "day"), data = data, mean))
daily_DEWP <- ts(ag1$DEWP, start = decimal_date(min(ag1$`floor_date(data$date, unit = "day")`)), frequency = 365)

ag1 <- as.data.frame(aggregate(TEMP ~ floor_date(data$date, unit = "day"), data = data, mean))
daily_TEMP <- ts(ag1$TEMP, start = decimal_date(min(ag1$`floor_date(data$date, unit = "day")`)), frequency = 365)

ag1 <- as.data.frame(aggregate(PRES ~ floor_date(data$date, unit = "day"), data = data, mean))
daily_PRES <- ts(ag1$PRES, start = decimal_date(min(ag1$`floor_date(data$date, unit = "day")`)), frequency = 365)

### correlation
par(mfrow=c(2,2))

ccf(logdaily_PM2.5, logdaily_WSPM)
ccf(logdaily_PM2.5, daily_TEMP)
ccf(logdaily_PM2.5, daily_DEWP)
ccf(logdaily_PM2.5, daily_PRES)

ccf(daily_TEMP, logdaily_PM2.5)
ccf(logdaily_WSPM,logdaily_PM2.5)
ccf(daily_DEWP, logdaily_PM2.5)
ccf(daily_PRES, logdaily_PM2.5)

prewhiten(logdaily_PM2.5, logdaily_WSPM, ylab = "CCF between logdaily_PM2.5 and logdaily_WSPM")
prewhiten(logdaily_PM2.5, daily_TEMP, ylab = "CCF bewtween logdaily_PM2.5 and daily_TEMP")
prewhiten(logdaily_PM2.5, daily_DEWP, ylab = "CCF between logdaily_PM2.5 and daily_DEWP")
prewhiten(logdaily_PM2.5, daily_PRES, ylab = "CCF bewtween logdaily_PM2.5 and daily_PRES")


## distribution of daily data
par(mfrow=c(2,3))
hist(logdaily_PM2.5)
hist(logdaily_WSPM)
hist(daily_DEWP)
hist(daily_TEMP)
hist(daily_PRES)

par(mfrow=c(2,2))
ts.plot(logdaily_PM2.5)
acf(logdaily_PM2.5, 365)
pacf(logdaily_PM2.5)

ts.plot(logdaily_WSPM)
acf(logdaily_WSPM, 365)
pacf(logdaily_WSPM)

ts.plot(daily_TEMP)
acf(daily_TEMP,365)
pacf(daily_TEMP)

ts.plot(daily_DEWP)
acf(daily_DEWP, 365)
pacf(daily_DEWP)

ts.plot(daily_PRES)
acf(daily_PRES, 365)
pacf(daily_PRES)

### fit best arima
logPM2.5_train <- data.frame(logdaily_PM2.5[1:1430])
fit.best.arima = function(logPM2.5_train, maxord=c(1,1,1)){
  best.aic = Inf
  n        = length(logPM2.5_train)
  for(p in 0:maxord[1]){
    for(d in 0:maxord[2]){
      for(q in 0:maxord[3]){
        mod = arima(logPM2.5_train, order=c(p,d,q),method="ML")
        #fit.aic = -2*mod$loglik + (log(n)+1)*length(mod$coef)  # consistent AIC
        fit.aic = mod$aic
        if(fit.aic<best.aic){best.aic=fit.aic
        best.fit = mod
        best.model= c(p,d,q)}
      }
    }
  }
  list(best.aic,best.fit,best.model)
}

fit1 = fit.best.arima(logPM2.5_train)

### fit best arima (ARMA(1,1))
logPM2.5_model1 = arima(logPM2.5_train, order = c(1,0,1), method = "ML")
logPM2.5_model1 

plot(logdaily_PM2.5[1:1430], type="l")
lines(logdaily_PM2.5[1:1430]-residuals(logPM2.5_model1), ylab="logdaily_PM2.5", col="red")
title(main="ARMA(1,1)")


par(mfrow=c(2,2))
plot(residuals(logPM2.5_model1))
acf(residuals(logPM2.5_model1), 100)
pacf(residuals(logPM2.5_model1), 100)
shapiro.test(residuals(logPM2.5_model1))
qqnorm(residuals(logPM2.5_model1))
qqline(residuals(logPM2.5_model1))

detectIO(logPM2.5_model1)
detectAO(logPM2.5_model1)

LB.test(logPM2.5_model1)

### fit linear regression model (daily)
data_filled <- data[,c(6, 12:14, 17)]
attach(data_filled)
logPM2.5_model <- lm(logdaily_PM2.5 ~ logdaily_WSPM+daily_TEMP+daily_DEWP+daily_PRES)

summary(logPM2.5_model)
par(mfrow=c(2,2))
plot(residuals(logPM2.5_model), type="l")
acf(residuals(logPM2.5_model), 100)
pacf(residuals(logPM2.5_model), 100)
eacf(residuals(logPM2.5_model))
adf.test(residuals(logPM2.5_model1))

### ARIMAX
logPM2.5_model2=arima(logdaily_PM2.5[1:1430],order=c(2,1,1),
                      xreg=data.frame(daily_DEWP[1:1430], logdaily_WSPM[1:1430], 
                                      daily_TEMP[1:1430], daily_PRES[1:1430]))

plot(logdaily_PM2.5[1:1430], type="l")
lines(logdaily_PM2.5[1:1430]-residuals(logPM2.5_model2), ylab="logdaily_PM2.5", col="red")
title(main="ARIMAX")

logPM2.5_model2
par(mfrow=c(2,2))
plot(residuals(logPM2.5_model2))
acf(residuals(logPM2.5_model2),100)
pacf(residuals(logPM2.5_model2),100)
qqnorm(residuals(logPM2.5_model2))
qqline(residuals(logPM2.5_model2))

detectAO(logPM2.5_model2)
detectIO(logPM2.5_model2)

### Forecast
## model1
forecast1 = predict(logPM2.5_model2, newxreg = data.frame(daily_DEWP[1431:1461], logdaily_WSPM[1431:1461], 
                                                          daily_TEMP[1431:1461], daily_PRES[1431:1461]), n.ahead = 31)

par(mfrow=c(2,2))
plot(logdaily_PM2.5[1281:1430], type="l", xlim=c(0,180), ylim=c(1,7.5))
points(151:181, logdaily_PM2.5[1431:1461], pch=5)
lines(151:181, forecast1$pred, col="blue")
lines(151:181, c(2.973741, 1.272025, 1.839234, 1.625463, 2.588281, 2.991227,
                 3.419626, 2.678687, 2.090556, 3.415541, 2.364693, 1.471679,
                 1.311773, 2.229217, 2.494843, 2.567485, 3.64573, 3.805754,
                 2.93695, 1.542632, 2.647353, 2.873888, 1.269223, 2.840613,
                 3.394138, 2.442489, 2.661151, 2.713197, 2.994619, 3.235433,
                 2.210632), col="red", lty=3)

lines(151:181, c(4.981096, 3.740352, 4.449873, 4.31206, 5.343256, 5.817141,
                 6.317679, 5.648002, 5.12948, 6.522403, 5.537981, 4.710016,
                 4.613879, 5.593886, 5.920934, 6.053917, 7.191476, 7.409837,
                 6.598441, 5.260645, 6.421042, 6.702443, 5.151869, 6.776606,
                 7.382765, 6.483065, 6.753015, 6.855716, 7.187179, 7.477446,
                 6.501527), col="red", lty=3)

title(main="Prediction of logdaily_PM2.5 (ARIMAX)")


## model2
forecast2 = predict(logPM2.5_model1, n.ahead = 31)

par(mfrow=c(2,2))
plot(logdaily_PM2.5[1281:1430], type="l", xlim=c(0,180))
points(151:181, logdaily_PM2.5[1431:1461], pch=5)
lines(151:181, forecast2$pred, col="blue")
lines(151:181, c(3.649564, 2.493387, 2.209602, 2.130245, 2.107183, 2.100401,
                 2.0984, 2.097809, 2.097634, 2.097582, 2.097567, 2.097563,
                 2.097561, 2.097561, 2.097561, 2.097561, 2.097561, 2.097561,
                 2.097561, 2.097561, 2.097561, 2.097561, 2.097561, 2.097561,
                 2.097561, 2.097561, 2.097561, 2.097561, 2.097561, 2.097561,
                 2.097561), col="red", lty=3)

lines(151:181, c(6.712726, 6.094753, 5.854204, 5.7786, 5.755865, 5.749112,
                 5.747113, 5.746522, 5.746347, 5.746296, 5.74628, 5.746276,
                 5.746274, 5.746274, 5.746274, 5.746274, 5.746274, 5.746274,
                 5.746274, 5.746274, 5.746274, 5.746274, 5.746274, 5.746274,
                 5.746274, 5.746274, 5.746274, 5.746274, 5.746274, 5.746274, 
                 5.746274), col="red", lty=3)

title(main="Prediction of logdaily_PM2.5 ARMA(1,1)")



