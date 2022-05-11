####################  USEFUL CODE FOR APPLIED STATISTICS #################### 
########################## FUNCTIONAL DATA ANALYSIS #########################

rm(list=ls())

#### setting directories ####

setwd('/Users/gildamatteucci/OneDrive - Politecnico di Milano/UNIVERSITA/APPLIED_STATISTICS/LABORATORY')
setwd('C:/Users/dario/OneDrive - Politecnico di Milano/LABORATORY')


#### Libraries ####

library(fda)             # regression splines, smoothing splines
library(KernSmooth)      # locpoly() function (local polynomial regression)
library(rgl)             # 3d plot
library(fdakma)          # fda kma

load('Functional Data Analysis/Lab_FDA/FDAfunctions.RData')
load('FDAfunctions.RData')

#### Smoothing: Central Finite Differences ####

x <- data$... #************************************************************************************************
y <- data$... #************************************************************************************************
n <- length(x)

# if you have available the true curve:
x_true <- data$... #************************************************************************************************
y_true <- data$... #************************************************************************************************
y1_true <- data$... #************************************************************************************************
y2_true <- data$... #************************************************************************************************


# graphical representation:
plot(x,y,xlab="x",ylab="observed data", cex=0.7, pch=16)

# Central Finite Differences:
y_fd1 <- (y[3:n]-y[1:(n-2)])/(x[3:n]-x[1:(n-2)])
y_fd2 <- ((y[3:n]-y[2:(n-1)])/(x[3:n]-x[2:(n-1)]) - 
         (y[2:(n-1)]-y[1:(n-2)])/(x[2:(n-1)] - x[1:(n-2)]))*2/(x[3:(n)]-x[1:(n-2)])

par(mfrow=c(1,3))
plot(x,y,xlab="t",ylab="observed data")
# points(x_true,y_true,type='l',col="orange",lwd=3) # true curve 
plot(x[2:(n-1)], y_fd1,xlab="t",ylab="first differences x",type="l")
# points(x_true, y1_true,type='l',col="orange",lwd=3) # true first derivative
plot(x[2:(n-1)], y_fd2,xlab="t",ylab="second differences x",type="l")
# points(x_true, y2_true,type='l',col="orange",lwd=3) # true second derivative










#### Smoothing: Regression Splines ####

x <- data$... #************************************************************************************************
y <- data$... #************************************************************************************************
n <- length(x)
## NOTA: se y sono piu funzioni, y=matrice in cui OGNI COLONNA e' UNA FUNZIONE

# Set parameters
m <- ...       # spline order     #************************************************************************************************
degree <- m-1  # spline degree    #************************************************************************************************
nbasis <- ...  # number of basis  #************************************************************************************************
nbreaks <- ... # number of breaks #************************************************************************************************
breaks <- ...  # OPTIONAL, locations of the knots (if not provided, they are equispaced) #*****************************************
## you can set nbasis OR alternatively breaks. Default:
# nbasis = nbreaks + norder - 2,
# nbreaks = nbasis - norder + 2,

# create bspline basis:
basis <- create.bspline.basis(rangeval=range(x), nbasis=nbasis, norder=m) #breaks = breaks
# plot(basis) # plot the basis

# # Evaluate the basis on the grid of abscissa
# basismat <- eval.basis(x, basis)
# basismat1 <- eval.basis(x, basis, Lfdobj=1) # evaluate the first derivate
# basismat2 <- eval.basis(x, basis, Lfdobj=2) # evaluate the second derivate
# 
# # compute the coefficients of the splines:
# coef_splines <- lsfit(basismat, y, intercept=FALSE)$coef
# coef_splines1 <- lsfit(basismat1, y, intercept=FALSE)$coef # first derivative
# coef_splines2 <- lsfit(basismat2, y, intercept=FALSE)$coef # second derivative
#
# # evaluation of the regression splines
# y_splines <- basismat %*% coef_splines
# y_splines1 <- basismat1 %*% coef_splines1 # first derivative
# y_splines2 <- basismat2 %*% coef_splines2 # second derivative

# Smooth the data (evaluate the basis system):
fd.object <- smooth.basis(argvals=x, y=y, fdParobj=basis)
y_splines <- eval.fd(x, fd.object$fd)         # the curve smoothing the data
y_splines1 <- eval.fd(x, fd.object$fd, Lfd=1) # first derivative
y_splines2 <- eval.fd(x, fd.object$fd, Lfd=2) # second derivative
df <- fd.object$df                            # degrees of freedom in the smoothing curve  

# graphical representation of regression splines:
plot(x,y,xlab="t",ylab="observed data", cex=0.7)
points(x,y_splines ,type="l",col="blue",lwd=2)
abline(v=basis$params)

# plot the derivatives:
par(mfrow=c(2,2))
plot(x,y,xlab="t",ylab="observed data")
points(x,y_splines ,type="l",col="blue",lwd=2)
plot(x,y_splines1 ,type="l",col="blue",lwd=2, xlab="t",ylab="first differences x")
# points(x[2:(n-1)],y_fd1,type="l") # add first finite difference
plot(x,y_splines2 ,type="l",col="blue",lwd=2, xlab="t",ylab="second differences x")
# points(x[2:(n-1)],y_fd2,type="l") # add second finite difference
plot(basis)



## Approximate pointwise confidence intervals
basismat <- eval.basis(x, basis)
S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) # projection operator 
sigmahat <- sqrt(sum((y_splines-y)^2)/(n-nbasis))           # estimate of sigma
lb <- y_splines - qnorm(0.975)*sigmahat*sqrt(diag(S))              # lower bound of CI
ub <- y_splines + qnorm(0.975)*sigmahat*sqrt(diag(S))              # upper bound of CI
# plot of the estimated smooth line
plot(x,y_splines,type="l",col="blue",lwd=2,ylab="")
points(x,lb,type="l",col="dodgerblue",lty="dashed")
points(x,ub,type="l",col="dodgerblue",lty="dashed")
# points(x, y_true, type="l")


## generalized cross-validation to choose the optimal nbasis:
m <- ... # spline order  (m = degree + 1) #***************************************************
nbasis <- 6:30 #******************************************************************************
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(range(x), nbasis[i], m)
  gcv[i] <- smooth.basis(x, y, basis)$gcv
}
nbasis_opt <- nbasis[which.min(gcv)] 
plot(nbasis,gcv)
points(nbasis_opt, min(gcv), pch=19)
nbasis_opt
# choose the nbasis which minimizes the gcv.


## Bias-Variance tradeoff (if you know the true curve):
basismat <- eval.basis(x, basis)
S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat) # projection operator 
sigmahat <- sqrt(sum((y_splines-y)^2)/(n-nbasis))           # estimate of sigma
sigma <- sigmahat                         # true sigma: estimated before as sigmahat
nbasis <- 9:15    #**************************************************************************
integrationinterval <- 10:(length(x)-10)  # exclude the extremes of the domain for boundary effects
bias <- rep(NA,len=length(nbasis))
var <- rep(NA,len=length(nbasis))
for (j in 1:length(nbasis)){
  basis <- create.bspline.basis(range(x), nbasis[j], m)
  basismat <- eval.basis(x, basis)
  S <- basismat%*%solve(t(basismat)%*%basismat)%*%t(basismat)
  bias[j] <- sum((y_true-S%*%y_true)[integrationinterval])
  var[j] <- (sigma^2)*sum(diag(S[integrationinterval,integrationinterval]))
}
mse <- var+bias^2                 # mean square error
# graphical representation
plot(nbasis,bias^2,ylim=c(0,max(mse)),type="l",ylab="",main="Bias-Variance tradeoff")
points(nbasis,var,col="red",type="l")
points(nbasis,mse,col="green",type="l",lwd=3)
legend('topright', c("Bias","Var","MSE"), col=c("black","red","green"), 
       lty=1, cex=.5)






## Smoothing Splines with Penalization term: 
der_penalized <- 3 #  order of the derivative to be penalized #**********************************
lambda <- ... # smoothing parameter #************************************************************
functionalPar <- fdPar(fdobj=basis, Lfdobj=der_penalized, lambda=lambda)  
smoothing <- smooth.basis(x, y, functionalPar) # smoothing splines estimation

# evaluation of the smoothing splines and its derivates
y_smoothing <- eval.fd(x, smoothing$fd, Lfd=0)
y_smoothing1 <- eval.fd(x, smoothing$fd, Lfd=1)
y_smoothing2 <- eval.fd(x, smoothing$fd, Lfd=2)

# graphical representation
par(mfrow=c(2,2))
plot(x,y,xlab="t",ylab="observed data")
points(x,y_smoothing ,type="l",col="blue",lwd=2)
plot(x,y_smoothing1 ,type="l",col="blue",lwd=2, xlab="t")
# points(x[2:(n-1)],y_fd1,type="l",ylab="first differences x") first finite difference
plot(x,y_smoothing2 ,type="l",col="blue",lwd=2, xlab="t",)
# points(x[2:(NT-1)],rappincX2,type="l", ylab="second differences x")
plot(basis)

# choose lambda via generalized cross-validation:
der_penalized <- 3 #  order of the derivative to be penalized #**********************************
lambda <- c(1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12) #*********************************************
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=der_penalized, lambda=lambda[i])  
  gcv[i] <- smooth.basis(x, y, functionalPar)$gcv
}
lambda_opt <- lambda[which.min(gcv)]
plot(log10(lambda),gcv)
# choose lambda which minimizes gcv



#### SMOOTHING: Local Polynomial Regression ####

x <- data$... #************************************************************************************************
y <- data$... #************************************************************************************************
n <- length(x)

# Set parameters
m <- ...       # order of the polynomial   #**************************************************
degree <- m-1  # degree of the polynomial  #**************************************************
bw <- 0.05     # bandwidth #******************************************************************
## NB: try to use different values of bandwidth

localpoly <- locpoly(x, y, degree=degree, bandwidth=bw, gridsize=length(x), range.x=range(x))
y_localpoly <- localpoly$y
localpoly1 <- locpoly(x, y, degree=degree, bandwidth=bw, gridsize=length(x), range.x=range(x), 
                      drv=1) # first derivative
y_localpoly1 <- localpoly1$y # first derivative
localpoly2 <- locpoly(x, y, degree=degree, bandwidth=bw, gridsize=length(x), range.x=range(x), 
                      drv=2) # second derivative
y_localpoly2 <- localpoly2$y # second derivative

# plot
par(mfrow=c(1,3))
plot(x,y,xlab="t",ylab="observed data")
points(x,y_localpoly ,type="l",col="blue",lwd=2)
plot(x,y_localpoly1 ,type="l",col="blue",lwd=2, xlab='t')
# points(x[2:(n-1)],y_fd1,ylab="first differences x",type="l") # first finite difference
plot(x,y_localpoly2 ,type="l",col="blue",lwd=2, xlab='t')
# points(x[2:(n-1)],y_fd2,ylab="second differences x",type="l") # second finite difference



#### SMOOTHING: Constrained Functions - Positive Curves ####

# \[ y_j = exp(w(t_j)) + e_j\]
# \[ f(t) = exp(w(t)) \]
# where the function w(t) is unconstrained and the function f(t) is monotone increasing.
# w(t) is modeled via a basis expansion: numerical methods are used to compute the coefficients of the basis expansion.

# You can use the function: smooth.pos().

x <- data$... #************************************************************************************************
y <- data$... #************************************************************************************************
n.x <- length(x)
n.y <- 1
## N.B.: if you have more than 1 curve:  #*********************************************************************
# y <- data[,c(...)]
# n.y <- dim(y)[2]

x_range <- range(x)
n_fine <- ... #************************************************************************************
x_fine <- seq(x_range[1], x_range[2], length=n_fine)

norder <- ... #************************************************************************************
nbasis <- n.x - 2 + norder
wbasis <- create.bspline.basis(rangeval = x_range, nbasis = nbasis,
                               norder = norder, breaks = x)

# construct the functional parameter with penalization of the third derivative
Lfdobj <- 3 # derivative to penalize  #************************************************************************************
lambda <- 10^(-0.5)  #************************************************************************************

cvecf <- matrix(0, nbasis, n.y) # this is used as initial value
# for the numerical techniques

Wfd0 <- fd(coef = cvecf, basisobj = wbasis)
yfdPar <- fdPar(fdobj = Wfd0, Lfdobj = Lfdobj, lambda = lambda)

# carry out a monotone smoothing
yPos <- smooth.pos(argvals = x, y = y, WfdParobj = yfdPar)

Wfd  <- yPos$Wfd
beta <- yPos$beta 
y_hatfd <- yPos$yhatfd

y1_pos <- deriv.fd(expr = y_hatfd, Lfdobj = 1) 
y1_pos.mean <- mean.fd(y1_pos)

y2_pos <- deriv.fd(expr = y_hatfd, Lfdobj = 2)
y2_pos.mean <- media.fda(y2_pos)


par(mfrow=c(2,2))
plot(y_hatfd, xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="y")
plot(y1_pos,  xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="First Derivative")
plot(y2_pos,  xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="Second Derivative")
plot(wbasis)






#### SMOOTHING: Constrained Functions - Monotone Curves ####

# I need to build a model for monotone curves:
#   \[ f(t) = \int_{t_0}^t exp(w(u)) du\]
#   \[y_j = b_0 + b_1 * f(t_j) + e_j\]
# where: 
#   
# * The function w(t) is unconstrained
# * The function f(t) is monotone increasing
# * $b_1$>0 for monotone increasing functions
# * $b_1$<0 for monotone decreasing functions
# * $b_0$ is the value of the function at $t_0$
#   
#   w(t) is modeled via a basis expansion numerical methods are used to compute the coefficients of the basis expansion, as well as $b_0$, $b_1$.

x <- data$... #************************************************************************************************
y <- data$... #************************************************************************************************
n.x <- length(x)
n.y <- 1
## N.B.: if you have more than 1 curve:  #*********************************************************************
# y <- data[,c(...)]
# n.y <- dim(y)[2]

x_range <- range(x)
n_fine <- ... #************************************************************************************
x_fine <- seq(x_range[1], x_range[2], length=n_fine)

norder <- ... #************************************************************************************
nbasis <- n.x - 2 + norder
wbasis <- create.bspline.basis(rangeval = x_range, nbasis = nbasis,
                               norder = norder, breaks = x)

# construct the functional parameter with penalization of the third derivative
Lfdobj <- 3 # derivative to penalize  #************************************************************************************
lambda <- 10^(-0.5)  #************************************************************************************

cvecf <- matrix(0, nbasis, n.y) # this is used as initial value
                                # for the numerical techniques

Wfd0 <- fd(coef = cvecf, basisobj = wbasis)
yfdPar <- fdPar(fdobj = Wfd0, Lfdobj = Lfdobj, lambda = lambda)

# carry out a monotone smoothing
yMon <- smooth.monotone(argvals = x, y = y, WfdParobj = yfdPar)

Wfd  <- yMon$Wfd
beta <- yMon$beta 
y_hatfd <- yMon$yhatfd

y1_mon <- deriv.fd(expr = y_hatfd, Lfdobj = 1) 
y1_mon.mean <- mean.fd(y1_mon)

y2_mon <- deriv.fd(expr = y_hatfd, Lfdobj = 2)
y2_mon.mean <- media.fda(y2_mon)


par(mfrow=c(2,2))
plot(y_hatfd, xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="y")
plot(y1_mon,  xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="First Derivative")
plot(y2_mon,  xlim=x_range, lty=1, lwd=2, cex=2, xlab="x", ylab="Second Derivative")
plot(wbasis)


#### SMOOTHING: Fourier Basis ####

# Choose Fourier Basis if data are periodic (i.e. years...)

x <- data$...  #********************************************************************************
y <- data$...  #************************************************************************************************
n.x <- length(x)
n.y <- 1
## N.B.: if you have more than 1 curve:  #*********************************************************************
# y <- data[,c(...)]
# n.y <- dim(y)[2]

nbasis <- ... #********************************************************************************
fourier.basis <- create.fourier.basis(rangeval=range(x), nbasis=nbasis) 
y_fd <- Data2fd(y = y, argvals = x,basisobj = fourier.basis)
plot.fd(y_fd)


## Estimate the Mean and the Covariance Kernel:
par(mfrow=c(1,2))
# mean
plot.fd(y_fd) 
lines(mean.fd(y_fd), lwd=3) 
# covariance
y.fourier <- eval.fd(x,y_fd) 
image.plot(x,x,(cov(t(y.fourier))[1:n.x,]))


#### Extension to Multidimensional Curves #####

# See end of lav Smoothing










#### Functional PCA ####



n.x <- ...  #********************************************************************************
n.y <- ...  #********************************************************************************

x <- as.numeric(c(1:n.x))
#x <- data$...  #********************************************************************************
y <- t(data[1:n.y,]) #Se ogni funzione ? una riga
y<- data[, 1:n.y] # Se ogni funzione ? una colonna
#y <- data$...  #************************************************************************************************

## N.B.: if you have more than 1 curve:  #*********************************************************************
# y <- data[,c(...)]
# n.y <- dim(y)[2]

nbasis <- ... #********************************************************************************
## Choose which basis you want (bspline or fourier):
# basis <- create.fourier.basis(rangeval=range(x), nbasis=nbasis) 
basis <- create.bspline.basis(rangeval=range(x), nbasis=nbasis, norder=m) #breaks = breaks
y_fd <- Data2fd(y = y, argvals = x,basisobj = basis)
plot.fd(y_fd)
y.fit <- eval.fd(x,y_fd) 

# performing PCA:
nharm <- ... # number of harmonics/components to compute (ne bastano 5) #*************************************************
pca <- pca.fd(y_fd, nharm=nharm, centerfns=TRUE) # set centerfns=T if you want to center functions

par(mfrow=c(1,2))
# eigenvalue
plot(pca$values[1:n.y],xlab='j',ylab='Eigenvalues')
# screeplot 
plot(cumsum(pca$values)[1:n.y]/sum(pca$values),xlab='j',ylab='CPV',ylim=c(0.8,1))
# choose the number of components to keep (elbow/threshold...)

#Proportion of Variance explained by each PC
pca$varprop
#Sum of the Proportion of variance explain by the first k PCs:
#(cumsum(pca$values)[1:n.y]/sum(pca$values))[k]

k <- ... # number of PC to keep  #*************************************************
# first k FPCs (loadings) (eigenfunctions)
par(mfrow=c(1,k))
for (i in 1:k){
  plot(pca$harmonics[i,], ylab=paste('FPC',i), col=rainbow(k)[i])
  lines(pca$harmonics[i,], ylab=paste('FPC',i), col=rainbow(k)[i])
  abline(h=0, lty = 2)
}

# to understand what PC are dealing with, I plot them with Sample Mean:
par(mfrow = c(1,k))
media <- mean.fd(y_fd)
for (i in 1:k){
  plot(media,lwd=2,main=paste('FPC',i)) 
  lines(media+pca$harmonics[i,]*sqrt(pca$values[i]), col=rainbow(3)[2]) 
  lines(media-pca$harmonics[i,]*sqrt(pca$values[i]), col=rainbow(3)[3])
}

# alternative plot:
par(mfrow=c(1,k))
plot.pca.fd(pca, nx=100, pointplot=TRUE, harm=1:k, expand=0, cycle=FALSE)


# scatter plot of the scores
#par(mfrow=c(1,2))
#plot(pca$scores[,1],pca$scores[,2],xlab="Scores FPC1",ylab="Scores FPC2",lwd=2)  #Without identifier
plot(pca$scores[,1],pca$scores[,2],type="n",xlab="Scores FPC1", ylab="Scores FPC2")
text(pca$scores[,1],pca$scores[,2], colnames(y), cex=1)
# from this plot you can spot outliers.

outliers <- c(,,,) # columns that seem outliers #*********************************************
matplot(y.fit,type='l')
lines(y.fit[,outliers],lwd=4, col=2) # profiles for outliers





#### K-means Clustering ####

x <- data$... # abscissas #*************************************************************************************
#x <- data.frame(1:80) Oppure se non ? nel data
y0 <- data[,c(,,,)]# evaluations of original functions #*************************************************************************************
y1 <- data[,c(,,,)]# evaluations of original functions’ first derivatives #*************************************************************************************


# graphical exploration of data
par(mfrow=c(1,2))
# Plot of original functions
matplot(x,y0, type='l', xlab='x', ylab='orig.func') #if it does not plot use t(x) and t(y0)
title('Original functions')
# Plot of original function first derivatives 
matplot(x,y1, type='l', xlab='x', ylab='orig.deriv') #if it does not plot use t(x) and t(y1)
title('Original function first derivatives')


n.clust <- ... #************************************************************************************
warping.method <- 'NOalignment' # 'affine', 'shift', 'dilation' #**************************
similarity.method <- 'd0.pearson' # 'd1.pearson', 'd0.L2', 'd1.L2', 'd0.L2.centered', 'd1.L2.centered'
#(d0.pearson is the correlation between the curves)
center.method <- 'k-means'   # 'k-medoids'  
set.seed(1) 

fdakma <- kma( x=t(x), y0=t(y0), #y1=y1,     # give y1 ONLY IF warping.method is d1
               n.clust=n.clust,             #Se da errore sulle dimensioni prova x=x, y0=y0
               warping.method = warping.method,
               similarity.method = similarity.method,
               center.method = center.method,
               #,seeds = c(1,11,21) # you can give a little help to the algorithm...)
)

my.kma.show.results(fdakma)

# access to labels assigned by kma:
fdakma$labels


# with my.kma.compare() you can set different warping methods amd different values of k
n.clust <- 1:3 #******************************************************************************
warping.method <- c('NOalignment') #,'affine', 'shift', 'dilation' #**************************
similarity.method <- 'd0.pearson' # 'd1.pearson', 'd0.L2', 'd1.L2', 'd0.L2.centered', 'd1.L2.centered'
center.method <- 'k-means'   # 'k-medoids'  
fdakma.compare <- my.kma.compare( x=x, y0=y0, #y1=y1,     # give y1 ONLY IF warping.method is d1
                                  n.clust=n.clust, 
                                  warping.method = warping.method,
                                  similarity.method = similarity.method,
                                  center.method = center.method,
                                  #,seeds = c(1,11,21) # you can give a little help to the algorithm...)
                                  plot.graph=1
                                  )



#### ciao ####