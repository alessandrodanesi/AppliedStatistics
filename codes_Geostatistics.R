####################  USEFUL CODE FOR APPLIED STATISTICS #################### 
################### GEOSTATISTICS AND FUNCTIONAL ANALYSIS ################### 

rm(list=ls())

#### setting directories ####

setwd('/Users/gildamatteucci/OneDrive - Politecnico di Milano/UNIVERSITA/APPLIED_STATISTICS/LABORATORY')
setwd('C:/Users/dario/OneDrive - Politecnico di Milano/LABORATORY')


#### Libraries ####

library(sp)       # meuse data
library(geoR)
library(gstat)

#### Functions for graphics ####
v.f <- function(x, ...){100-cov.spatial(x, ...)}
v.f.est <- function(x,C0, ...){C0-cov.spatial(x, ...)}


#### setting data ####
data(meuse)
data <- meuse

# select variable on which you want to make prediction:
data$target <- data$... #****************************************************************************
data$dist   <- data$... #****************************************************************************
# cancella le variabili che hai sovrascritto (?)




# define coordinates of the data:
coordinates(data) <- c('x','y')

# first visualization of data:
bubble(data,'target',do.log=TRUE,key.space='right', maxsize = 2.5, main = '')

# # if you have polygon:
# data.lst <- list(Polygons(list(Polygon(data.pol)), "data.pol")) 
# data.sr <- SpatialPolygons(data.lst)


## exploratory analysis:
# histogram: 
# par(mfrow=c(1,2))
hist(data$target, breaks=16, col="grey", main='',
     prob = TRUE, xlab = 'target variable')
# hist(log(data$target), breaks=16, col="grey", main='', 
#      prob = TRUE, xlab = 'log of target variable')
## if you have a concentration (or in general a non symmetric histogram),
## maybe it's better to take the log of it:
# data$target <- log(data$target)

# scatterplot of target variable with respect to 
xyplot(target ~ dist, as.data.frame(data), pch = 16)


#### variogram analysis ####

# Stationary Model:
# \[ F(s_i) = m + \delta_{s_i}, \quad s_i \in D \]
# where: $m = E[F_{s_i}]$ is the drift and $\delta_{s_i}$ is the stochastic residual. 

# NON-Stationary Model:
# \[ F(s_i) = m_{s_i} + \delta_{s_i}, \quad s_i \in D \]
# where: $m_{s_i} = E[F_{s_i}]$ is the drift and $\delta_{s_i}$ is the stochastic residual. 


# Note: if you have as regressor a dummy variable, you can predict the coefficients 
# of your model by making prediction (BLUE = T) on a random point with certain value of dummy.
# i.e.: beta0 = drift when dummy = 0
#       beta1 = (drift when dummy = 1) - (drift when dummy = 0)




# sampling the variogram:
v <- variogram(target ~ 1, data)
# try to sample the stationary model and non stationary model to see what is the best.
# in a good model I expect an asymptot.

# plot the estimated variogram:
plot(v, main = 'Sample Variogram', pch = 19)
## NB: you can set
#     boundaries = vector of points you want to estimate
#     width = distance between two estimated points
#     cutoff = maximum distance up to which variogram values are computed
#               (American - default: bbox diagonal / 3)
#               (French: smaller side of bbox / 2)
# 

# if you want to check anisotropy (alpha is the angle you want to check):
# plot(variogram(target ~ 1, data, alpha = c(0, 45, 90, 135)),pch=19) # alpha is 



# modeling the variogram:
v.mod <- vgm(psill=..., model=..., range=..., nugget=...)
# sill = value of horizontal asymptot 
# model = 'Exp': linear near 0, then it stabilyzes
#         'Sph': linear near 0, then it stabilyzes
#         'Gau': quadratic near 0
#         'Mat': big family (don't use it.)
# range = value from which variogram start to converge
# nugget = intercept of the variogram (jump of variogram in the origin)
# add.to = you can add the vgm to a specific variog (type the existing variog)
# anis = 	anisotropy parameters. You have to provide c(par1, par2) where:
#           par1 is angle for the principal direction of continuity 
#           par2 is anisotropy ratio

## NB:
# pure nugget: flat variogram (no spatial dependence)
# linear model: variogram is a line, it has NO horizontal asymptot (range is infinity)
## NB2: if you add a nugget, r will nest two models: one Nug and the other as you set


# fit the variogram:
v.fit <- fit.variogram(v, v.mod) # fit.method = ...
# other useful arguments:
#     fit method: default is weighted least squares. see details of fit.variogram
#     fit.sills = c(FALSE,TRUE): first logical referred to Nug (if present), second to your specific model.
#                               if you do NOT have set nugget, just put ONE logical
#                               if TRUE it will take your sill as initial value and make it converges.
#                               if FALSE il will take your sill as "fixed"
#     fit.range = c(FALSE,TRUE): first logical referred to Nug (if present), second to your specific model.
#                                if you do NOT have set nugget, just put ONE logical
#                                if TRUE it will take your range as initial value and make it converges.
#                                if FALSE il will take your range as "fixed"
v.fit
## NB: check if sill or nugget are negative. In this case, change initial values. 


# fit the variogram through maximum likelihood: 
v.fit.ml <- fit.variogram.reml(formula, data, model = v.mod)
# formula is the formula you have used to define v (i.e. target ~ 1)

# plot 
plot(v, v.fit, pch = 19)


# # fitting method: non linear regression with minimization of weighted 
# # sum of squares error. final value of the minimum
# attr(v.fit, 'SSErr')


#### prediction in a new location ####

# coordinates of the new location
s0.new <- data.frame(x=..., y=...) #, dist=...) # it can be a grid of points #************************************************************
coordinates(s0.new)=c('x','y')

# plot 
plot(data)
plot(s0.new, add = TRUE, col = 'red', lwd = 2)


# # if your model is NOT stationary and you want to compute distance for example:
# s0.vec <- as.vector(slot(s0.new,'coords'))
# # distance to the river: calculate the distance between s0 and s0.star 
# # (it can be a specific point or a grid of point - cfr i.e. meuse points)
# s0.dist <- min(rowSums(scale(s0.star, s0.vec)^2))
# s0.new <- as.data.frame(c(s0.new,s0.dist))
# names(s0.new) <- c('x','y','dist')
# coordinates(s0.new) <- c('x','y')
# s0.new <- as(s0.new, 'SpatialPointsDataFrame')
# s0.new


# Create a gstat object setting a spherical (residual) variogram.
# formula is the formula you have used to define v (i.e. target ~ 1)
g.tr <- gstat(formula = ..., data = data, model = v.mod) #*****************************************

# pointwise prediction:
s0.pred <- predict(g.tr, s0.new, BLUE = FALSE)
# it returns: 
# universal/ordinary kriging, new coordinates, prediction on new coordinates, kriging variance (var of the pred error)

## N.B. 
# if you are using a stationary model, this is (Simple/)Ordinary Kriging. 
# if you are using a NON stationary model, this is Universal Kriging. 
# if you set BLUE = FALSE, you are making prediction on the given point. 
# if you set BLUE = TRUE, you are making prediction on the mean at that location (using nearby locations)
# (basically you are estimating the drift. If the model is stationary it's constant on D).
# I can use BLUE = TRUE if I want to use other points to predict a target in a new location 
# in a situation that is independent from the one I used to build the model. 
# If you are computing UNIVERSAL kriging, you do NOT have to look at the kriging variance! 
# If s0.new is a grid of points, if you set BLUE=TRUE you are estimating the mean over the grid.


# if you are predicting a grid of points, you can plot them: 
spplot(s0.pred[,1], main = 'prediction')
spplot(s0.pred[,2], main = 'variance')



## comparing stationary variogram and non stationary variogram of residuals:
model <- c(...,...) # 'Gau', 'Sph', 'Exp', 'Nug'
sostituisci <- function(model){
   n <- length(model)
   model.new <- rep(0, times = n)
   if (n ==1){
      if (model == 'Gau'){
         model.new <- 'gaussian'
      }
      if (model == 'Exp'){
         model.new <- 'exponential'
      }
      if (model == 'Sph'){
         model.new <- 'spherical'
      }
      if (model == 'Nug'){
         model.new <- 'pure.nugget'
      }
   }
   if (n ==2){
      if (model[1]=='Gau' | model[2]=='Gau'){
         model.new[which(model=='Gau')] <- 'gaussian'
      }
      if (model[1]=='Exp' | model[2]=='Exp'){
         model.new[which(model=='Exp')] <- 'exponential'
      }
      if (model[1]=='Sph' | model[2]=='Sph'){
         model.new[which(model=='Sph')] <- 'spherical'
      }
      if (model[1]=='Nug' | model[2]=='Nug'){
         model.new[which(model=='Nug')] <- 'pure.nugget'
      }
   }
   return(model.new)
}

v.stat   <- variogram(target ~ 1, data)
v.nostat <- variogram(target ~ dist, data)
v.mod   <- vgm(psill=..., model=..., range=..., nugget=...) #************************************************************
v.stat.fit   <- fit.variogram(v.stat, v.mod)
v.nostat.fit <- fit.variogram(v.nostat, v.mod)

model <- as.character(v.mod$model)
model.new <- sostituisci(model)
psill.stat <- v.stat.fit$psill
range.stat <- v.stat.fit$range
psill.nostat <- v.nostat.fit$psill
range.nostat <- v.nostat.fit$range
yrange <- c(0,sum(psill.stat))


# Compare the variogram of the data and the variogram of the residuals
# red: stationary; blue: non stationary
plot(v.stat$dist,v.stat$gamma, xlab='distance', ylab='semivariance',pch=19,col='red',ylim=yrange)
curve(v.f.est(x, C0=psill.stat[1]+psill.stat[2], 
              cov.pars=rbind(c(psill.stat[1], range.stat[1]),
                             c(psill.stat[2], range.stat[2])), 
              cov.model = model.new), from = 0.0001, # to = 1600,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "Variogram model", add=TRUE,col='red', lwd=2) # ,ylim=c(0,110)
points(v.nostat$dist,v.nostat$gamma,xlab='distance',ylab='semivariance',pch=19,col='steelblue',ylim=yrange)
curve(v.f.est(x, C0=psill.nostat[1]+psill.nostat[2], 
              cov.pars=rbind(c(psill.nostat[1], range.nostat[1]), 
                             c(psill.nostat[2], range.nostat[2])), 
              cov.model = model.new), from = 0.0001, #to = 1600,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "Variogram model",add=TRUE,col='steelblue',lwd=2) # ,ylim=c(0,110))




