library(carData)
library(fitdistrplus)
library(logspline)
library(bootstrap)
library("tidyverse")
set.seed(29)
dep <- Ginzberg$adjdep

# Cullen & Frey Graph
descdist(dep, discrete = FALSE,boot = 10000)


gamma <- fitdist(dep, "gamma")
lognormal <- fitdist(dep, "lnorm")
normal <- fitdist(dep, "norm")
weibull <- fitdist(dep, "weibull")

# Fit Analysis
plot(gamma)
plot(lognormal)
plot(normal)
plot(weibull)

# AIC Test
gamma$aic
lognormal$aic
normal$aic
weibull$aic

gamma

# Kolmogorov-Smirnov Test 
dep.shape <- gamma$estimate[1]
dep.rate <- gamma$estimate[2]
n <- length(dep)
ks.sim <- rep(0,10000)
for(i in 1:10000)
{ 
  resample.dep <- rgamma(n, shape = dep.shape, rate = dep.rate)
  ks.sim[i] <- as.numeric(ks.test(resample.dep,"pgamma",
                                  shape = dep.shape,rate = dep.rate)$statistic)
} 


kstest.fit <- logspline(ks.sim)
1 - plogspline(ks.test(resample.dep,"pgamma",
                       shape = dep.shape,rate = dep.rate)$statistic
               ,kstest.fit)
par(mfrow = c(1,1))
plot(ecdf(ks.sim), 
     las = 1,main = "KS-test statistic simulation (CDF)", col = "blue")

# ECDF and Plotting of Confidence Interval
dep.ecdf <- ecdf(dep)
plot(dep.ecdf, 
     las = 1,
     main = "ECDF",
     col = "green",
     cex=0.2
)
Alpha=0.05
Eps=sqrt(log(2/Alpha)/(2*n))
x<-seq(min(dep),max(dep), length.out = 10000)
lines(x, pmin(dep.ecdf(x)+Eps,1),col='red')
lines(x, pmax(dep.ecdf(x)-Eps,0),col='blue')


sum(pgamma(x,dep.shape,dep.rate) >=pmax(dep.ecdf(x)-Eps,0) & 
      pgamma(x,dep.shape,dep.rate) <=pmin(dep.ecdf(x)+Eps,1))/100

#Parametric and Non Parametric Bootstrappin for parameter(mean) estimation
mu.cap <- dep.shape/dep.rate
sd.cap <- sd(dep)
B <- 10000
mu.cap.star.p <- rep(0,B)
mu.cap.star.np <- rep(0,B)
for(i in 1:B) { 
  x.p <- rgamma(n,shape = dep.shape, rate = dep.rate)
  pfit <- fitdist(x.p, "gamma")
  mu.cap.star.p[i] <- pfit$estimate[1] / pfit$estimate[2]
  x.np <- sample(dep,size=n,replace=TRUE)
  npfit <- fitdist(x.np, "gamma")
  mu.cap.star.np[i] <- npfit$estimate[1] / npfit$estimate[2]
}
se.cap.boot.p <- sd(mu.cap.star.p)
se.cap.boot.np <- sd(mu.cap.star.np)

theta.np <- as.data.frame(mu.cap.star.np)
theta.p <- as.data.frame(mu.cap.star.p)
names(theta.np) <- "val"
names(theta.p) <- "val"
theta.np$type <- 'Parametric Bootstrap'
theta.p$type <- 'Non-Parametric Bootstrap'
dis <- rbind(theta.p,theta.np)
dis <- as.data.frame(dis)
#Bootstrap distribution of estimate (mean)
ggplot(dis, aes(val, fill = type)) + geom_density(alpha = 1) +
  scale_fill_grey() +
  ggtitle("Parametric vs Non-Parametric Distribution of mean")

normal.p<-c(mu.cap-2*se.cap.boot.p, mu.cap+2*se.cap.boot.p)
#normal ci, at 95%
normal.np<-c(mu.cap-2*se.cap.boot.np, mu.cap+2*se.cap.boot.np) 
#pivotal ci at 95%
pivotal.np<-c(2*mu.cap-quantile(mu.cap.star.np,0.975),
              2*mu.cap-quantile(mu.cap.star.np,0.025))
#quantile ci at 95%
quantile.np<-quantile(mu.cap.star.np, c(0.025, 0.975))


x1 <- Ginzberg$adjdep[Ginzberg$adjsimp <= 1]
x2 <- Ginzberg$adjdep[Ginzberg$adjsimp > 1]

#permutation test
perm.matrix <- t(replicate(500, sample(72)))


perm.T<-apply(
  perm.matrix,1,
  function(x){
    abs(mean(Ginzberg$adjdep[x[1:46]])-mean(Ginzberg$adjdep[x[47:72]]))
  })
perm.T


p.value<-mean(perm.T>abs(mean(x1)-mean(x2)))
p.value

#Wilcox test
wilcox.test(x1, x2,conf.int = T,exact=F)

#Simple Linear Regression
lm.fit <- lm(depression~adjsimp, Ginzberg)
summary(lm.fit)






install.packages("DT")