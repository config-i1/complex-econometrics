x <- complex(real=Mcomp::M3[[2830]]$x,imaginary=Mcomp::M3[[2831]]$x)

obs <- 1000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=1.5*x0+rnorm(length(x0),0,10))
x <- complex(real=rnorm(1000,10,10),imaginary=rnorm(1000,10,10))
y <- complex(real=rnorm(1000,10,10),imaginary=rnorm(1000,10,10))
save(x, y, file="x.Rdata")

cor(data.frame(y1=Re(y),y2=Im(y),x1=Re(x),x2=Im(x)))
ccor(x,y)

test <- cacf(x)
which(abs(1-Re(test$acf))<0.2)
plot(test,1)
plot(test,2)

layout(matrix(c(1,2,1,3),2,2))
test <- ccf(Re(x),Im(x))
test <- acf(Re(x))
test <- acf(Im(x))


x <- c(1+2i, 2+3i, 3+4i, 5+6i, 7+8i)
a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,10))

plot(data.frame(y1=Re(y),y2=Im(y),x1=Re(x),x2=Im(x)))
# plot(y)
# text(Re(y), Im(y), which(y %in% y), pos=2)

ccor(x,y)

ccov(x,y) / cvar(x)
ccov(x,y) / cvar(y)

mean(y) - ccov(x,y) / cvar(x) * mean(x)

xreg <- as.matrix(data.frame(Intercept=1, x=x))

# OLS
invert(t(xreg) %*% xreg) %*% t(xreg) %*% y

complexData <- cbind(y,x)

invert(t(xreg) %*% xreg) * sigma(test)

test <- clm(y~x, complexData, subset=c(1:100))
summary(test)

lossFunction <- function(actual,fitted,B,xreg){
    error <- actual-fitted
    return(Re(sum(error * Conj(error))));
}

lossFunction <- function(actual,fitted,B,xreg){
    error <- actual-fitted
    return(sum(abs(error)));
}

test <- clm(y~x, complexData, loss="CLS")
test1 <- clm(y~x, complexData, loss="OLS")
test2 <- clm(y~x, complexData, subset=c(1:100), loss=lossFunction, print_level=3)
test3 <- clm(y~x, complexData, subset=c(1:100), loss=lossFunction, print_level=3)

ccov(complexData)[1,2] / cvar(complexData)[2,2]

b <- mean((complexData[,1]-mean(complexData[,1]))*Conj((complexData[,2]-mean(complexData[,2])))) /
    mean(abs(complexData[,2]-mean(complexData[,2]))^2)
a <- mean(complexData[,1]) - b*mean(complexData[,2])

errors <- complexData[,1] - (a+ b*complexData[,2])

mean(resid(test))
sd(Re(resid(test)))
sd(Im(resid(test)))
hist(Re(resid(test)))
hist(Im(resid(test)))

mean(errors)
sd(Re(errors))
sd(Im(errors))
hist(Re(errors))
hist(Im(errors))

summary(Re(resid(test)))
summary(Re(errors))

qqnorm(Re(resid(test)))
qqline(Re(resid(test)))

qqnorm(Re(errors))
qqline(Re(errors))




#### Do a simulation experiment to check whether the estimates of parameters are efficient and unbiased for the two methods!
### Uncorrelated xr and xi, Same Variances
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=rnorm(length(x0),0,10))

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersUCSV <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersUCSV[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersUCSV[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersUCSV[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersUCSV[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCSV[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersUCSV[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCSV[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersUCSV[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCSV[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


### Uncorrelated xr and xi, Different Variances
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=rnorm(length(x0),0,20))

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersUCDV <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersUCDV[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersUCDV[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersUCDV[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersUCDV[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCDV[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersUCDV[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCDV[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersUCDV[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersUCDV[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")



### Perfectly correlated xr and xi
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=1.5*x0)

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersPC <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersPC[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersPC[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersPC[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersPC[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPC[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersPC[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPC[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersPC[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPC[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")



### Perfectly correlated xr and xi, but with an intercept shift
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=1.5*x0+2)

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 1.5*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersPCInt <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersPCInt[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersPCInt[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersPCInt[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersPCInt[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPCInt[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersPCInt[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPCInt[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersPCInt[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersPCInt[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


### Correlated, High variance
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=0.5*x0 + rnorm(length(x0),0,15))

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 10*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersHV <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersHV[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersHV[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersHV[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersHV[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHV[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersHV[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHV[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersHV[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHV[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


### Functional relation
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=0.5*x0 + rnorm(length(x0),0,15))

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x

complexData <- cbind(y,x)

nsim <- 9980
parametersFR <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersFR[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersFR[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersFR[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersFR[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersFR[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersFR[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(1.5,2.5))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersFR[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-2,-1))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersFR[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersFR[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


### High Variance of error
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=0.5*x0 + rnorm(length(x0),0,15))

a0 <- 10 + 15i
a1 <- 2-1.5i
y <- a0 + a1 * x + 100*complex(real=rnorm(length(x),0,1), imaginary=rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersHVError <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersHVError[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersHVError[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersHVError[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersHVError[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHVError[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-3,0))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersHVError[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHVError[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-3,0))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersHVError[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersHVError[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


### Correlated errors
obs <- 10000
x0 <- rnorm(obs,10,10)
x <- complex(real=x0,imaginary=0.5*x0 + rnorm(length(x0),0,15))

a0 <- 10 + 15i
a1 <- 2-1.5i
error <- rnorm(length(x),0,1)
y <- a0 + a1 * x + 10*complex(real=error, imaginary=0.8*error+rnorm(length(x),0,1))

complexData <- cbind(y,x)

nsim <- 9980
parametersCorError <- array(NA, c(nsim,3,2), dimnames=list(NULL, c("CLS","OLS","Likelihood"), c("a0","a1")))

for(i in 1:nsim){
    test <- clm(y~x, complexData, loss="CLS", subset=sample(c(1:obs),20+i))
    parametersCorError[i,1,] <- coef(test)
    test <- clm(y~x, complexData, loss="OLS", subset=sample(c(1:obs),20+i))
    parametersCorError[i,2,] <- coef(test)
    test <- clm(y~x, complexData, loss="likelihood", subset=sample(c(1:obs),20+i))
    parametersCorError[i,3,] <- coef(test)
}


par(mfcol=c(2,3))
plot(c(21:obs), Re(parametersCorError[,1,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="CLS", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersCorError[,1,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="CLS", ylim=c(-3,0))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersCorError[,2,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="OLS", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersCorError[,2,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="OLS", ylim=c(-3,0))
abline(h=Im(a1), col="red")

plot(c(21:obs), Re(parametersCorError[,3,2]), type="l", ylab="Re(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(0,4))
abline(h=Re(a1), col="red")
plot(c(21:obs), Im(parametersCorError[,3,2]), type="l", ylab="Im(a1)", xlab="Sample size",
     main="Likelihood", ylim=c(-3,0))
abline(h=Im(a1), col="red")


save(parametersUCSV, parametersUCDV, parametersPC, parametersPCInt,
     parametersHV, parametersFR, parametersHVError, parametersCorError,
     file="./data/clmExperimentsParametersNew.Rdata")


#### Correct vcov ####
test <- clm(y~x, complexData, loss="CLS", subset=c(1:30))
test1 <- clm(y~x, complexData, loss="OLS", subset=c(1:30))
test3 <- clm(y~x, complexData, loss="likelihood", subset=c(1:30))

vcov(test)
sigma(test)
summary(test)
summary(test1)
summary(test3)

test2 <- predict(test, tail(complexData,10), interval="pred")
test3 <- predict(test1, tail(complexData,10), interval="pred")
plot(test2)
plot(test3)

### Loss experiments

test <- clm(y~x, complexData, loss="CLS", subset=c(1:30))
Re(test$lossValue) + Im(test$lossValue)
B <- coef(test)
Bnew <- B
Bnew[1] <- Bnew[1] * (1-0.02)
Bnew[2] <- Bnew[2] * (1+0.01)
test1 <- clm(y~x, complexData, loss="CLS", subset=c(1:30), parameters=Bnew)
Re(test1$lossValue)*Im(test1$lossValue)
plot(resid(test1))

lossFunction <- function(actual,fitted,B,xreg){
    CFValue <- sum((actual - fitted)^2);
    CFValue <- Re(CFValue) + Im(CFValue);
}
test2 <- clm(y~x, complexData, loss=lossFunction, subset=c(1:30))
test2$lossValue
plot(resid(test2))
