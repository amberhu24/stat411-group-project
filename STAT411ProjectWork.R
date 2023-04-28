rm(list = ls())



# FIT 4 Linear Models




library(R2jags)
library(coda)

hist(y,xlab='Life Expectancy (Years)')

data=read.csv("C:\\Users\\17039\\Downloads\\Life Expectancy Data.csv",header=TRUE,stringsAsFactors=FALSE)
data=data[data$Year==2014,]
data=na.omit(data)
data$Hepatitis.B=100-data$Hepatitis.B
data$Polio=100-data$Polio
data$Diphtheria=100-data$Diphtheria




data$percentage.expenditure=log(data$percentage.expenditure)
data$GDP=log(data$GDP)
data$Population=log(data$Population)
data$thinness..1.19.years=log(data$thinness..1.19.years)
data$thinness.5.9.years=log(data$thinness.5.9.years)
data$Hepatitis.B=log(data$Hepatitis.B)
data$Polio=log(data$Polio)
data$Diphtheria=log(data$Diphtheria)

country=data$Country
data=subset(data, select = -c(Country,Year,Status) )





y=data$Life.expectancy
x=data.frame(subset(data, select = -c(Life.expectancy)))

N = length(y)
K = ncol(x)


model=lm(y~as.matrix(x))
summary(model)
confint(model,level=0.95)
informationMLR=c(AIC(model),BIC(model))
names(informationMLR)=c('AIC','BIC')

informationMLR

library(glmnet)
cv_model = cv.glmnet(as.matrix(x), y, alpha = 1)
lassomodel=glmnet(as.matrix(x),y,alpha=1,lambda=cv_model$lambda.min)
summary(lassomodel)
coef(lassomodel)

lassomodelMLR=lm(y~x$Adult.Mortality+x$Alcohol+x$percentage.expenditure+x$Measles+x$BMI+x$Total.expenditure+x$Diphtheria+x$HIV.AIDS+x$GDP+x$Population+x$thinness.5.9.years+x$Income.composition.of.resources)
summary(lassomodelMLR)
informationlassomodelMLR=c(AIC(lassomodelMLR),BIC(lassomodelMLR))
names(informationlassomodelMLR)=c('AIC','BIC')

informationlassomodelMLR


y_predicted <- predict(lassomodel, s = cv_model$lambda.min, newx = as.matrix(x))
rsq = 1 - (sum((y_predicted - y)^2)/sum((y - mean(y))^2))
rsq

SSX=sum(as.matrix(x)^2)
g=5
Xo=sqrt(SSX/g)

x=subset(x,select=c(Adult.Mortality,Alcohol,percentage.expenditure,Measles,BMI,Total.expenditure,Diphtheria,HIV.AIDS,GDP,Population,thinness.5.9.years,Income.composition.of.resources)
)
N = length(y)
K = ncol(x)
jagsData = with(data, list(y = y, x = as.matrix(x), nX =K , n = N,SSX=sum(as.matrix(x)^2)))

jags.MLRinits <- function(){ 
  return(list(beta0 = runif(1,0,1), beta = runif(K,0,1), tau = runif(1,0.01,0.05)))
}
jags.MLRparams = c("beta0", "beta", "sigma2")
MLRmodel <- function(){
  for (i in 1:n) {
    y[i]~dnorm(mu[i],tau)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  beta0 ~ dnorm(0,1.0E-6)
  for (j in 1:nX) {
    beta[j] ~ dnorm(0,1.0E-6)
  }
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
}

MLRmodelSensitive <- function(){
  for (i in 1:n) {
    y[i]~dnorm(mu[i],tau)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  beta0 ~ dnorm(0,1.0E-3)
  for (j in 1:nX) {
    beta[j] ~ dnorm(0,1.0E-6)
  }
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
}




jags.MLRinitsGPrior <- function(){ 
  return(list(beta0 = runif(1,0,1), beta = runif(K,0,1), tau = runif(1,0.01,0.05),lambda=rgamma(1,0.5,0.5)))
}
jags.MLRparamsGPrior = c("beta0", "beta", "sigma2",'lambda')
MLRmodelGPrior <- function(){
  for (i in 1:n) {
    y[i]~dnorm(mu[i],tau)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  
  beta0 ~ dnorm(0,tau*lambda*SSX/n)
  
  for (j in 1:nX) {
    beta[j] ~ dnorm(0,tau*lambda*SSX/n)
  }
  lambda~dgamma(0.5,0.5)
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
}

MLRmodelGPriorSensitive <- function(){
  for (i in 1:n) {
    y[i]~dnorm(mu[i],tau)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  
  beta0 ~ dnorm(0,2*tau*lambda*SSX/n)
  
  for (j in 1:nX) {
    beta[j] ~ dnorm(0,tau*lambda*SSX/n)
  }
  lambda~dgamma(0.5,0.5)
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
}



jags.MLRinitsTDistribution <- function(){ 
  return(list(beta0 = runif(1,0,1), beta = runif(K,0,1), tau = runif(1,0.01,0.05),nuMinusOne=20))
}
jags.MLRparamsTDistribution=c('beta0','beta','sigma2','nu')
MLRmodelTDistribution <- function(){
  for (i in 1:n) {
    y[i]~dt(mu[i],tau,nu)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  beta0 ~ dnorm(0,1.0E-6)
  for (j in 1:nX) {
    beta[j] ~ dnorm(0,1.0E-6)
  }
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
  nuMinusOne~dexp(1/29)
  nu=nuMinusOne+1
}



MLRmodelLASSO <- function(){
  for (i in 1:n) {
    y[i]~dnorm(mu[i],tau)
    mu[i] = beta0 + inprod(beta[],x[i,])
  }
  beta0 ~ dnorm(0,1.0E-6)
  for (j in 1:nX) {
    beta[j] ~ ddexp(0,lambda)
  }
  lambda~dunif(0.001,10)
  tau~dgamma(0.001, 0.001)
  sigma2 = 1/tau
}









jagsfitMLR = jags(
  data = jagsData, 
  inits = jags.MLRinits, 
  parameters.to.save = jags.MLRparams, 
  model.file = MLRmodel,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)

jagsfitMLRLASSO = jags(
  data = jagsData, 
  inits = jags.MLRinits, 
  parameters.to.save = jags.MLRparams, 
  model.file = MLRmodelLASSO,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)

jagsfitMLRGPrior = jags(
  data = jagsData, 
  inits = jags.MLRinitsGPrior, 
  parameters.to.save = jags.MLRparamsGPrior, 
  model.file = MLRmodelGPrior,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)

jagsfitMLRTDistribution = jags(
  data = jagsData, 
  inits = jags.MLRinitsTDistribution, 
  parameters.to.save = jags.MLRparamsTDistribution, 
  model.file = MLRmodelTDistribution,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)


# ----------------------------

jagsfitMLRSensitive = jags(
  data = jagsData, 
  inits = jags.MLRinits, 
  parameters.to.save = jags.MLRparams, 
  model.file = MLRmodelSensitive,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)

jagsfitMLRGPriorSensitive = jags(
  data = jagsData, 
  inits = jags.MLRinitsGPrior, 
  parameters.to.save = jags.MLRparamsGPrior, 
  model.file = MLRmodelGPriorSensitive,
  n.chains = 1, 
  n.iter = 11000, 
  n.burnin = 1000, 
  n.thin = 1)


#VARIABLES FOR LASSO
names(data)[2:length(names(data))]


jagsfitMLR
jagsfitMLRLASSO
jagsfitMLRTDistribution
jagsfitMLRGPrior

jagsfitMLR$BUGSoutput$summary[,1]/jagsfitMLRSensitive$BUGSoutput$summary[,1]

jagsfitMLRGPrior$BUGSoutput$summary[,1]/jagsfitMLRGPriorSensitive$BUGSoutput$summary[,1]

notMLR=jagsfitMLR$BUGSoutput$summary[jagsfitMLR$BUGSoutput$summary[,7]>0 & jagsfitMLR$BUGSoutput$summary[,3]<0,c(1,2,3,7)]
notMLR
notMLRLASSO=jagsfitMLRLASSO$BUGSoutput$summary[jagsfitMLRLASSO$BUGSoutput$summary[,7]>0 & jagsfitMLRLASSO$BUGSoutput$summary[,3]<0,c(1,2,3,7)]
notMLRLASSO

jagsfitMLRLASSO$BUGSoutput$summary[,1]/jagsfitMLR$BUGSoutput$summary[,1]

library("bayesplot")
library("ggplot2")
library("rstanarm") 

library(mcmcplots)

jagsfit=jagsfitMLR
jagsfit.mcmc = as.matrix(as.mcmc(jagsfit))
del_col = which(colnames(jagsfit.mcmc)=="deviance")
output = jagsfit.mcmc[,-del_col]
outputMCMC=as.mcmc(output)

mcmc_intervals(outputMCMC, pars = c('beta[12]','beta0','sigma2'))
mcmc_intervals(outputMCMC, pars = c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]'))


denplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))
traplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))

#------------------------

jagsfit=jagsfitMLRLASSO
jagsfit.mcmc = as.matrix(as.mcmc(jagsfit))
del_col = which(colnames(jagsfit.mcmc)=="deviance")
output = jagsfit.mcmc[,-del_col]
outputMCMC=as.mcmc(output)

mcmc_intervals(outputMCMC, pars = c('beta[12]','beta0','sigma2'))
mcmc_intervals(outputMCMC, pars = c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]'))


denplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))
traplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))

#-----------------------

jagsfit=jagsfitMLRTDistribution
jagsfit.mcmc = as.matrix(as.mcmc(jagsfit))
del_col = which(colnames(jagsfit.mcmc)=="deviance")
output = jagsfit.mcmc[,-del_col]
outputMCMC=as.mcmc(output)

mcmc_intervals(outputMCMC, pars = c('beta[12]','beta0','sigma2'))
mcmc_intervals(outputMCMC, pars = c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]'))


denplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))
traplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))

#------------------------

jagsfit=jagsfitMLRGPrior
jagsfit.mcmc = as.matrix(as.mcmc(jagsfit))
del_col = which(colnames(jagsfit.mcmc)=="deviance")
output = jagsfit.mcmc[,-del_col]
outputMCMC=as.mcmc(output)

mcmc_intervals(outputMCMC, pars = c('beta[12]','beta0','sigma2'))
mcmc_intervals(outputMCMC, pars = c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]'))



denplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))
traplot(outputMCMC,parms=c('beta[1]','beta[2]','beta[3]','beta[4]','beta[5]','beta[6]','beta[7]','beta[8]','beta[9]','beta[10]','beta[11]','beta[12]','beta0','sigma2'))

#----------------------
samples=jags.samples(jagsfitMLRGPrior$model,c('WAIC','deviance'),type='mean',n.iter=5000,n.burnin=1000,n.thin=1)
samples$p_waic=samples$WAIC
samples$WAIC=samples$deviance+samples$p_waic
temp=sapply(samples,sum)
effective=temp[3]
aic=temp[2]+2*K
bic=temp[2]+log(N)*K
waic=temp[1]
dic=temp[2]
vars=c(aic,bic,waic,dic)
names(vars)=c('AIC','BIC','WAIC','DIC')
print(vars)

mcmc1=jagsfitMLR$BUGSoutput$sims.matrix

info1=infoJAGS(jagsfitMLR)
info2=infoJAGS(jagsfitMLRLASSO)
info3=infoJAGS(jagsfitMLRTDistribution)
info4=infoJAGS(jagsfitMLRGPrior)


library(BAS)

xx=x
LMaic=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='AIC')
LMaic
summary(LMaic)
coef(LMaic)
confint(coef(LMaic))
plot(LMaic,ask=F)

LMbic=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='BIC')
LMbic
summary(LMbic)
coef(LMbic)
confint(coef(LMbic))
plot(LMbic,ask=F)

LMbic=bas.lm(as.vector(unlist(y))~x$Adult.Mortality+x$percentage.expenditure+x$HIV.AIDS+x$GDP+x$Income.composition.of.resources,method='MCMC',MCMC.iterations=20000,prior='BIC')
LMbic=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='BIC')

LMbic
summary(LMbic)
coef(LMbic)
confint(coef(LMbic))
plot(LMbic,ask=F)

gpriorLM=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='g-prior',alpha=N)
gpriorLM
summary(gpriorLM)
coef(gpriorLM)
confint(coef(gpriorLM))
plot(gpriorLM,ask=F)

gpriorLMbigN=bas.lm(as.vector(unlist(y))~x$Adult.Mortality+x$percentage.expenditure+x$Total.expenditure+x$HIV.AIDS+x$GDP+x$Income.composition.of.resources,method='MCMC',MCMC.iterations=20000,prior='g-prior',alpha=N*10)
gpriorLMbigN=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='g-prior',alpha=N*10)
gpriorLMbigN
summary(gpriorLMbigN)
coef(gpriorLMbigN)
confint(coef(gpriorLMbigN))
plot(gpriorLMbigN,ask=F)

gpriorLMsmallN=bas.lm(as.vector(unlist(y))~as.matrix(x),method='MCMC',MCMC.iterations=20000,prior='g-prior',alpha=N/10)
gpriorLMsmallN
summary(gpriorLMsmallN)
coef(gpriorLMsmallN)
confint(coef(gpriorLMsmallN))
plot(gpriorLMsmallN,ask=F)





image(LMaic,rotate=F)
image(LMbic,rotate=F)
image(gpriorLM,rotate=F)
image(gpriorLMbigN,rotate=F)
image(gpriorLMsmallN,rotate=F)
