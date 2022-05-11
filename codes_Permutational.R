####################  USEFUL CODE FOR APPLIED STATISTICS #################### 
########################## PERMUTATIONAL STATISTICS #########################

rm(list=ls())

#### setting directories ####

setwd('/Users/gildamatteucci/OneDrive - Politecnico di Milano/UNIVERSITA/APPLIED_STATISTICS/LABORATORY')
setwd('C:/Users/dario/OneDrive - Politecnico di Milano/LABORATORY')


#### Libraries ####

#### Test for the Mean of a Univariate Population ####

# We want to test:
# \[H_0: E[Y]=\mu_0 \quad \vs \quad H_1: [Y] \neq \mu_0\]
# we add the assumption that data distribution f is symmetric (with respect to its mean \mu), namely:
# $f(y) = f(2\mu − y)$


x1 <- .... #************************************************************************************************

n <- dim(x1)[1]


mu0      <- c(...,...,...)  # dim=p #************************************************************************************************

x.mean   <- mean(x1)

# choose the test statistics
T20 <- as.numeric( abs(x.mean) - mu0 )


# Estimating the permutational distribution under H0
B <- 100000 #************************************************************************************
T2 <- numeric(B) 

for(perm in 1:B){
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1
  x_perm <- matrix(mu0,nrow=n,ncol=p,byrow=T) + (x1 - mu0) * matrix(signs.perm,nrow=n,ncol=p,byrow=FALSE)
  x.mean_perm <- colMeans(x_perm)

  T2[perm] <- as.numeric( abs(x.mean_perm-mu0) )
}

# plotting the permutational distribution under H0
hist(T2,xlim=range(c(T2,T20)),breaks=100)
abline(v=T20,col=3,lwd=4)
# if the line falls on the right of the curve, we can assume that the two distributions 
# are different. 

plot(ecdf(T2))
abline(v=T20,col=3,lwd=4)

# p-value
p_val <- sum(T2>=T20)/B
p_val
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that $E[Y]$ is different from mu0.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that $E[Y]$ is equal to mu0.










#### Test for the Mean of a Multivariate Population ####

# We want to test:
# \[H_0: E[Y]=\mu_0 \quad \vs \quad [Y] \neq \mu_0\]
# We have to suppose that the p-variate data distribution $f(y1,...,yp)$ 
# is symmetric with respect to the mean $\mu = (\mu_1 \dots \mu_p), 
# i.e. $f(y_1 \dots y_p) = f(2\mu_1 − y_1 \dots 2\mu_p − y_p)$


x1 <- .... #************************************************************************************************

p <- dim(x1)[2]
n <- dim(x1)[1]


mu0      <- c(...,...,...)  # dim=p #************************************************************************************************

x.mean   <- colMeans(x1)
x1.cov <- cov(x1)

# choose the test statistics
# T20 <- as.numeric((x.mean-mu0) %*% (x.mean-mu0) ) # if you have few data
T20 <- as.numeric((x.mean-mu0) %*% solve(x1.cov) %*% (x.mean-mu0) )


# Estimating the permutational distribution under H0
B <- 100000 #************************************************************************************
T2 <- numeric(B) 

for(perm in 1:B){
  signs.perm <- rbinom(n, 1, 0.5)*2 - 1
  x_perm <- matrix(mu0,nrow=n,ncol=p,byrow=T) + (x1 - mu0) * matrix(signs.perm,nrow=n,ncol=p,byrow=FALSE)
  x.mean_perm <- colMeans(x_perm)
  x.cov_perm <- cov(x1)
  
  # T2[perm]  <- (x.mean_perm-mu0)  %*% (x.mean_perm-mu0) # if you have few data
  T2[perm] <- as.numeric((x.mean_perm-mu0) %*% solve(x1.cov_perm) %*% (x.mean_perm-mu0) )
}

# plotting the permutational distribution under H0
hist(T2,xlim=range(c(T2,T20)),breaks=100)
abline(v=T20,col=3,lwd=4)
# if the line falls on the right of the curve, we can assume that the two distributions 
# are different. 

plot(ecdf(T2))
abline(v=T20,col=3,lwd=4)

# p-value
p_val <- sum(T2>=T20)/B
p_val
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that $E[Y]$ is different from mu0.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that $E[Y]$ is equal to mu0.















#### Two Independent Univariate Populations ####

# We want to test: 
# \[ H_0: X_1 \stackrel{d}{=} X_2 \quad vs \quad X_1 \stackrel{d}{\neq} X_2\]

x1 <- .... #************************************************************************************************
x2 <- .... #************************************************************************************************

n1 <- length(x1)
n2 <- length(x2)
n <- n1 + n2


## First visualization of data:
par(mfrow=c(1,2))
boxplot(x1,x2,main='Original data')


# # Number of distinct values of T*:
# factorial(n)/(2*factorial(n1)*factorial(n2))
# 
# # Minimun achieveable p-value:
# 1/(factorial(n)/(2*factorial(n1)*factorial(n2)))


# Test statistic: absolute difference between the two means
T0 <- abs(mean(x1) - mean(x2))

## setting parameters:
B <- 100000          # Number of permutations #*********************************************
T_stat <- numeric(B) # Vector where we will store the values of T*

for(perm in 1:B){
  # permutation:
  permutation <- sample(1:n)
  x_perm <- x_pooled[permutation]
  x1_perm <- x_perm[1:n1]
  x2_perm <- x_perm[(n1+1):n]
  # test statistic:
  T_stat[perm] <- abs(mean(x1_perm) - mean(x2_perm))
}

## Permutational distribution of T
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)
# if the line falls on the right of the curve, we can assume that the two ditributions 
# are different. 

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=2)

# p-value
p_val <- sum(T_stat>=T0)/B
p_val
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that the two distributions are different.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the two distributions are equal.








#### Two Independent Multivariate Populations ####

# We want to test: 
# \[ H_0: X_1 \stackrel{d}{=} X_2 \quad vs \quad X_1 \stackrel{d}{\neq} X_2\]

x1 <- .... #************************************************************************************************
x2 <- .... #************************************************************************************************

n1 <- dim(x1)[1]
n2 <- dim(x2)[2]
n <- n1 + n2


x1.mean <- colMeans(x1)
x2.mean <- colMeans(x2)
x1.cov <- cov(x1)
x2.cov <- cov(x2)
Sp <- ((n1-1)*x1.cov + (n2-1)*x2.cov)/(n1+n2-2) # Spooled matrix 

## Test Statistics:
# T20 <- as.numeric((x1.mean-x2.mean) %*% (x1.mean-x2.mean)) # if you have a small number of data
T20 <-((n1*n2)/(n))* as.numeric((x1.mean-x2.mean) %*% solve(Sp) %*% (x1.mean-x2.mean))


# # number of possible data point permutations 
# factorial(n)
# # number of different values of the test statistic
# choose(n,n1)


## Estimating the permutational distribution under H0
B <- 100000
T2 <- numeric(B)

for(perm in 1:B){
  x_pooled <- rbind(x1,x2)
  permutation <- sample(n)
  x_perm <- x_pooled[permutation,]
  x1_perm <- x_perm[1:n1,]
  x2_perm <- x_perm[(n1+1):n,]
  
  x1.mean_perm <- colMeans(x1_perm)
  x2.mean_perm <- colMeans(x2_perm)
  x1.cov <- cov(x1_perm)
  x2.cov <- cov(x2_perm)
  Sp <- ((n1-1)*x1.cov + (n2-1)*x2.cov)/(n1+n2-2) # Spooled matrix 
  T2[perm]  <- ((n1*n2)/(n))*(x1.mean_perm-x2.mean_perm) %*% solve(Sp)%*% (x1.mean_perm-x2.mean_perm)
  # T2[perm]  <- (x1.mean_perm-x2.mean_perm) %*% (x1.mean_perm-x2.mean_perm) # if you have a small number of data
}


# plotting the permutational distribution under H0
hist(T2,xlim=range(c(T2,T20)),breaks=1000)
abline(v=T20,col=3,lwd=4)
# if the line falls on the right of the curve, we can assume that the two ditributions 
# are different. 

plot(ecdf(T2))
abline(v=T20,col=3,lwd=4)

# p-value
p_val <- sum(T2>=T20)/B
p_val
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that the two distributions are different.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the two distributions are equal.


#### Two Paired Multivariate Populations ####

# We want to test: 
# \[ H_0: X_1 \stackrel{d}{=} X_2 \quad vs \quad X_1 \stackrel{d}{\neq} X_2\]

x1 <- .... #************************************************************************************************
x2 <- .... #************************************************************************************************

p <- dim(x1)[2]
n1 <- dim(x1)[1]
n2 <- dim(x2)[1]
n <- n1 + n2

delta.0 <- rep(0, times=p) #**********************************************************************

diff <- x1-x2
diff.mean <- colMeans(diff)
diff.cov <- cov(diff)
diff.invcov <- solve(diff.cov)

# choose the proper test statistics:
# T20 <- as.numeric(n1 * (diff.mean-delta.0)  %*% (diff.mean-delta.0))
# T20 <- as.numeric(n1 * (diff.mean-delta.0) %*% solve(diag(diag(diff.cov))) %*% (diff.mean-delta.0))
T20 <- as.numeric(n1 * (diff.mean-delta.0) %*% diff.invcov %*% (diff.mean-delta.0))


# # number of possible data point reflections 
# 2^n1


# Estimating the permutational distribution under H0
B <- 10000
T2 <- numeric(B)

for(perm in 1:B){
  # obs: exchanging data within couples means changing the sign of the difference
  signs.perm <- rbinom(n1, 1, 0.5)*2 - 1
  
  diff_perm <- diff * matrix(signs.perm,nrow=n1,ncol=p,byrow=FALSE)
  diff.mean_perm <- colMeans(diff_perm)
  diff.cov_perm <- cov(diff_perm)
  diff.invcov_perm <- solve(diff.cov_perm)
  
  # choose the proper test statistics:
  #T2[perm] <- as.numeric(n1 * (diff.mean_perm-delta.0) %*% (diff.mean_perm-delta.0))
  #T2[perm] <- as.numeric(n1 * (diff.mean_perm-delta.0) %*% solve(diag(diag(diff.cov_perm))) %*% (diff.mean_perm-delta.0))
  T2[perm] <- as.numeric(n1 * (diff.mean_perm-delta.0) %*% diff.invcov_perm %*% (diff.mean_perm-delta.0))
}

# plotting the permutational distribution under H0
hist(T2,xlim=range(c(T2,T20)),breaks=100)
abline(v=T20,col=3,lwd=4)
# if the line falls on the right of the curve, we can assume that the two ditributions 
# are different. 

plot(ecdf(T2))
abline(v=T20,col=3,lwd=4)


# p-value
p_val <- sum(T2>=T20)/B
p_val
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that the two distributions are different.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the two distributions are equal.



#### ANOVA ####

#I perform the Test:
#\[H_0 :  Y_1 \stackrel{d}{=} Y_2 \stackrel{d}{=} \dots \stackrel{d}{=} Y_g \quad vs \quad H_1 :  \exists  j_1, j_2 \in \{1, \dots, g\} \, st. \, Y_{j_1} \stackrel{d}{\neq} Y_{j_2}.\]

data <- as.data.frame(data)
n <- dim(data)[1]



factor1 <- factor(data$...)   # Categorical Variable ******************************************************************************************
x <- data$...   # Numerical Variable **********************************************************************************************************
g <- nlevels(factor1)

# first visualization of data:

plot(factor1, x, xlab='treat',col=rainbow(g),main='Data')



fit <- aov(x ~ factor1)
summary(fit)

# Permutation test:
# Test statistic: F stat
T0 <- summary(fit)[[1]][1,4]


# Permutations
B <- 10000
T_stat <- numeric(B)

for(perm in 1:B){
  # Permutation:
  permutation <- sample(1:n)
  x_perm <- x[permutation]
  fit_perm <- aov(x_perm ~ factor1)
  
  # Test statistic:
  T_stat[perm] <- summary(fit_perm)[[1]][1,4]
}


layout(1)
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)
# if the line falls on the right of the curve, we can assume that the two ditributions 
# are different. 

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)
  
pvalAN <- sum(T_stat>=T0)/B  
pvalAN
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that the two distributions are different.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the two distributions are equal.


#### MANOVA ####

data <- as.data.frame(data)

#I perform the Test:
#\[H_0 :  Y_1 \stackrel{d}{=} Y_2 \stackrel{d}{=} \dots \stackrel{d}{=} Y_g \quad vs \quad H_1 :  \exists  j_1, j_2 \in \{1, \dots, g\} \, st. \, Y_{j_1} \stackrel{d}{\neq} Y_{j_2}.\]

# labels:
factor1 <- factor(data$...) # Categorical Var. ****************************************************************************************************
# var numerica:
x <- data[, c(1, , , ...)]  # Numerical Var. (TOGLIERE CATEGORICA) ********************************************************************************

n <- dim(x)[1]
p <- dim(x)[2]
g <- nlevels(factor1)

#BoxPlot
par(mfrow=c(1,dim(x)[2]))
for(i in 1:dim(x)[2]){
  boxplot(x[,i]~factor1, main=colnames(x)[i], ylim=c(min(x[,i]),max(x[,i])), col = rainbow(g))
}


fit <- manova(x ~ factor1)
summary.manova(fit,test="Wilks")

# Permutation test:
# Test statistic: Wilks Lambda
T0 <- -summary.manova(fit,test="Wilks")$stats[1,2]


# Permutations
B <- 10000                              #***********************************************************************************************
T_stat <- numeric(B)

for(perm in 1:B){
  # choose random permutation
  permutation <- sample(1:n)
  factor1.perm <- factor1[permutation]
  fit.perm <- manova(x ~ factor1.perm)
  T_stat[perm] <- -summary.manova(fit.perm,test="Wilks")$stats[1,2]
}


layout(1)
hist(T_stat,xlim=range(c(T_stat,T0)),breaks=30)
abline(v=T0,col=3,lwd=2)
# if the line falls on the right of the curve, we can assume that the two ditributions 
# are different. 

plot(ecdf(T_stat))
abline(v=T0,col=3,lwd=4)

pvalMAN <- sum(T_stat>=T0)/B  
pvalMAN
# if pvalue is SMALL: we have evidence to reject H0, namely we can assume that the two distributions are different.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the two distributions are equal.


#### Two-Ways ANOVA ####

data <- as.data.frame(data)

# The model we have to test is the following:
# \[ x = \mu + \alpha factor_1 + \beta factor_2 + \gamma factor_1*factor_2 \]
# We have 3 different tests:
#
# 1. $factor_1  \quad (H_0: \alpha=0)$
# 2. $factor_2  \quad  (H_0: \beta=0)$
# 3. $factor_1 : factor_2 (Interaction)  \quad  (H_0: \gamma=0)$

p <- 1

# setting dataset:
# labels:
factor1 <- factor(data$......) # ******************************************************************************
factor2 <- factor(data$......) # ******************************************************************************
factor12 <- factor(paste(factor1, factor2, sep='-'))

# var numerica:
x <- data$.....        # ******************************************************************************

g <- length(levels(factor1))   # number of factor.1 
b <- length(levels(factor2))   # number of factor.2 
gb <- length(levels(factor12)) # number of factor.12 
n <- length(x)/(g*b)           # number of data

M         <- mean(x)                     # overall mean
Mfactor1  <- tapply(x, factor1, mean)    # mean for the types of factor1
Mfactor2  <- tapply(x, factor2, mean)    # mean for the types of factor2
Mfactor12 <- tapply(x, factor12, mean)   # mean for the all possible combinations factor1+factor2


#I check Interaction:
summary.aov(aov(x ~ factor1 + factor2 + factor1:factor2))
T0_interaction <- summary.aov(aov(x ~ factor1 + factor2 + factor1:factor2))[[1]][3,4]


# the idea is to permute the residuals under H0: $ x = \mu + \alpha factor_1 + \beta factor_2 $
# additive model

aov.H0interaction <- aov(x ~ factor1 + factor2)
aov.H0interaction
residuals.H0interaction <- aov.H0interaction$residuals

B <- 1000 # *********************************************************************************************
T_interaction <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  residuals.H0interaction <- residuals.H0interaction[permutation]
  x.perm.H0interaction <- aov.H0interaction$fitted + residuals.H0interaction
  T_interaction[perm] <- summary.aov(aov(x.perm.H0interaction ~ factor1 + factor2 + factor1:factor2))[[1]][3,4]
}

# p-value
pint <- sum(T_interaction >= T0_interaction)/B
pint
# if pvalue is SMALL: we have evidence to reject H0, namely we can not assume that the interaction is not significant.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the interaction is not significant.


# TEST OF FACTOR1   (H0: alpha=0)
T0_factor1 <- summary.aov(aov(x ~ factor1 + factor2))[[1]][1,4]
# residuals under H0:
# x = mu + beta*factor2
aov.H0factor1 <- aov(x ~ factor2)
residuals.H0factor1 <- aov.H0factor1$residuals

# TEST OF FACTOR2   (H0: beta=0)
T0_factor2 <- summary.aov(aov(x ~ factor1 + factor2))[[1]][2,4]
# residuals under H0:
# x = mu + alpha*factor1
aov.H0factor2 <- aov(x ~ factor1)
residuals.H0factor2 <- aov.H0factor2$residuals

# TEST OF FACTOR1 AND TEST OF FACTOR2
# p-values
B <- 10000      #***********************************************************************************************************
T_factor2 <- T_factor1 <- numeric(B)
for(perm in 1:B){
  permutation <- sample(n)
  
  x.perm.H0factor1 <- aov.H0factor1$fitted + residuals.H0factor1[permutation]
  T_factor1[perm] <- summary.aov(aov(x.perm.H0factor1 ~ factor1 + factor2))[[1]][1,4]
  
  x.perm.H0factor2 <- aov.H0factor2$fitted + residuals.H0factor2[permutation]
  T_factor2[perm] <- summary.aov(aov(x.perm.H0factor2 ~ factor1 + factor2))[[1]][2,4]
}

pfactor1 <- sum(T_factor1 >= T0_factor1)/B
pfactor2 <- sum(T_factor2 >= T0_factor2)/B
pfactor1
pfactor2
# if pvalue is SMALL: we have evidence to reject H0, namely we can not assume that the interaction is not significant.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the interaction is not significant.


#Then you have to perform a last ANOVA test on the remaining factor

#### Regression ####

data <- as.data.frame(data)



Y <- data$... #*******************************************************************************
n <- dim(data)[1]
reg1 <- data$... #*****************************************************************
reg2 <- data$... #*****************************************************************
reg3 <- data$... #*****************************************************************


## Model:
# \[ y = \beta_0 + \beta_1 * reg1 + \beta_2 * reg2 + \beta_3 * reg3 + \varepsilon \]


#Represent the data:
plot(reg1,Y,pch=16)
plot(reg2,Y,pch=16)
plot(reg3,Y,pch=16)


# permutation inference
# we want to perform different tests

# A) Overall model
# H0: beta1 = beta2 = beta3 = 0
# test statistic
T0_glob <- summary(result)$f[1]


# B) Test on variable reg1
# H0: beta1 = 0

# test statistic
T0_reg1 <- abs(summary(result)$coefficients[2,3])
# reduced model:
# Y = beta0 + beta2*reg2 + beta3*reg3
regr.H01 <- lm(Y ~ reg2 + reg3)
residui.H01 <- regr.H01$residuals

# C) Test on variable reg2
# H0: beta2 = 0

# test statistic
T0_reg2 <- abs(summary(result)$coefficients[3,3])
# reduced model:
# Y = beta0 + beta1*reg1 + beta3*reg3
regr.H02 <- lm(Y ~ reg1 + reg3)
residui.H02 <- regr.H02$residuals

# D) Test on variable reg3
# H0: beta3 = 0

# test statistic
T0_reg3 <- abs(summary(result)$coefficients[4,3])
# reduced model:
# Y = beta0 + beta1*reg1 + beta2*reg2
regr.H03 <- lm(Y ~ reg1 + reg2)
residui.H03 <- regr.H03$residuals


# p-values of the tests
B <- 1000
T_H0glob <- T_H01 <- T_H02 <- T_H03 <- numeric(B)

for(perm in 1:B){
  permutazione <- sample(n)
  
  Y.perm.glob <- Y[permutazione]
  T_H0glob[perm] <- summary(lm(Y.perm.glob ~ reg1 + reg2 + reg3))$f[1]
  
  residui.H01.perm <- residui.H01[permutazione]
  Y.perm.H01 <- regr.H01$fitted + residui.H01.perm
  T_H01[perm] <- abs(summary(lm(Y.perm.H01 ~ reg1 + reg2 + reg3))$coefficients[2,3])
  
  residui.H02.perm <- residui.H02[permutazione]
  Y.perm.H02 <- regr.H02$fitted + residui.H02.perm
  T_H02[perm] <- abs(summary(lm(Y.perm.H02 ~ reg1 + reg2 + reg3))$coefficients[3,3])
  
  residui.H03.perm <- residui.H03[permutazione]
  Y.perm.H03 <- regr.H03$fitted + residui.H03.perm
  T_H03[perm] <- abs(summary(lm(Y.perm.H03 ~ reg1 + reg2 + reg3))$coefficients[4,3])
  
}

pvalglob <- sum(T_H0glob>=T0_glob)/B

pvalreg1 <- sum(T_H01>=T0_reg1)/B
pvalreg2 <- sum(T_H02>=T0_reg2)/B
pvalreg3 <- sum(T_H03>=T0_reg3)/B

pvalglob
# if pvalue is SMALL: we have evidence to reject H0, namely we can not assume that all the regressors are not significant.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that all the regressors are not significant.

pvalreg1
pvalreg2
pvalreg3
# if pvalue is SMALL: we have evidence to reject H0, namely we can not assume that the regressor is not significant.
# if pvalue is LARGE: we do NOT have evidence to reject H0, namely we can assume that the regressor is not significant.

#### ciao ####