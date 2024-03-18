# package loading
library(tidyverse)
library(doRNG)
library(foreach)
library(doParallel)
library(parallel)
library(geepack)
library(lavaan)
library(Rfast)
library(WeightIt)
library(cobalt)

# lag1 wave2: --------------------------------------------------------------
#simple longitudinal setting with 1-period lagged effect and 2 waves of data
#estimating the cross-lagged effects under different data generation scenarios

fun_lag1_wave2 <- function(bxy, byx, n, gL, heterogeneity){
  # data generation
  n=n
  if(heterogeneity == FALSE){
    C0 = rnorm(n, 0, 1)
    L0 = rnorm(n, 0, 1)+0.1*C0
    gL0= do.call(gL, list(L0))
    X0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
    Y0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
    
    L1 = 0.1*C0 + 0.5*L0 + rnorm(n, 0, 1)
    gL1= do.call(gL, list(L1))
    X1 = 0.1*C0 + 0.5*X0 + byx*Y0 + 0.2*gL1 + rnorm(n, 0, 1)
    Y1 = 0.1*C0 + 0.5*Y0 + bxy*X0 + 0.2*gL1 + rnorm(n, 0, 1)
    
    L2 = 0.1*C0 + 0.5*L1 + rnorm(n, 0, 1)
    gL2= do.call(gL, list(L2))
    X2 = 0.1*C0 + 0.5*X1 + byx*Y1 + 0.2*gL2 + rnorm(n, 0, 1)
    Y2 = 0.1*C0 + 0.5*Y1 + bxy*X1 + 0.2*gL2 + rnorm(n, 0, 1)
    
  } else{
    
    C0 = rnorm(n, 0, 1)
    L0 = rnorm(n, 0, 1)+0.1*C0
    gL0= do.call(gL, list(L0))
    X0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
    Y0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
    
    L1 = 0.1*C0 + 0.5*L0 + rnorm(n, 0, 1)
    gL1= do.call(gL, list(L1))
    X1 = 0.1*C0 + 0.5*X0 + byx*Y0 + 0.1*gL1 + 0.3*Y0*gL1 + rnorm(n, 0, 1)
    Y1 = 0.1*C0 + 0.5*Y0 + bxy*X0 + 0.1*gL1 + 0.3*X0*gL1 + rnorm(n, 0, 1)
    
    L2 = 0.1*C0 + 0.5*L1 + rnorm(n, 0, 1)
    gL2= do.call(gL, list(L2))
    X2 = 0.1*C0 + 0.5*X1 + byx*Y1 + 0.1*gL2 + 0.3*Y1*gL2 + rnorm(n, 0, 1)
    Y2 = 0.1*C0 + 0.5*Y1 + bxy*X1 + 0.1*gL2 + 0.3*X1*gL2 + rnorm(n, 0, 1)
  }
  
  data <- data.frame(C0,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2) 
  cor(data)
  
  # CLPM
  res_data <- data.frame(X1=residuals(lm(X1 ~ C0 +L1)),
                         Y1=residuals(lm(Y1 ~ C0 +L1)),
                         X2=residuals(lm(X2 ~ C0 +L2)),
                         Y2=residuals(lm(Y2 ~ C0 +L2))) 
  model <- 
    "Y2 ~ a1*Y1 + b1*X1
     X2 ~ a2*X1 + b2*Y1
     X1 ~~r1*Y1
     X2 ~~0*Y2
     X1 ~~X1
     X2 ~~X2
     Y1 ~~Y1
     Y2 ~~Y2"
  fit <- sem(model, res_data)
  CLPM_result <- summary(fit, fit.measures = T, 
                         standardized = TRUE,
                         ci = TRUE,
                         rsquare = TRUE)
  CLPM <- CLPM_result[["pe"]][c(2,4), c("est", "se")]
  
  # G-CLPM
  SL.library <- c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.earth")

    fit_XY <- weightit(X1 ~ Y1 + L1 + C0, method = "super", 
                       SL.library = SL.library)
    fit_YX <- weightit(Y1 ~ X1 + L1 + C0, method = "super", 
                       SL.library = SL.library)
  
  w.cv_XY <- sd(fit_XY[["weights"]], na.rm = TRUE)/mean(fit_XY[["weights"]], na.rm = TRUE)
  if(!is.finite(w.cv_XY) | w.cv_XY > 4){
    fit_XY <- trim(fit_XY, at=0.1)
  }
  w.cv_YX <- sd(fit_YX[["weights"]], na.rm = TRUE)/mean(fit_YX[["weights"]], na.rm = TRUE)
  if(!is.finite(w.cv_YX) | w.cv_YX > 4){
    fit_YX <- trim(fit_YX, at=0.1)
  }
  
  bal.tab(fit_XY)
  bal.tab(fit_YX)
  
  G_CLPM_XY <- summary(geeglm(Y2 ~ X1, weights = fit_XY[["weights"]], id= rownames(data)))
  G_CLPM_YX <- summary(geeglm(X2 ~ Y1, weights = fit_YX[["weights"]], id= rownames(data)))
  
  
  G_CLPM <- rbind(G_CLPM_XY$coefficients[2, c("Estimate","Std.err")], 
                  G_CLPM_YX$coefficients[2, c("Estimate","Std.err")])
  
  colnames(G_CLPM) <- c("est", "se")
  
  result <- matrix(c(as.matrix(rbind(CLPM, G_CLPM))), nrow = 1)
  colnames(result) <- c("CLPM_est_xy","CLPM_est_yx","GCLPM_est_xy","GCLPM_est_yx",
                        "CLPM_se_xy","CLPM_se_yx","GCLPM_se_xy","GCLPM_se_yx")
  return(result)
}

clus <- makeCluster(20)
registerDoParallel(clus)

# different sample size
n_list <- c(200, 400, 600, 800, 1000, 1200, 1600, 2000, 2500, 3000, 3500, 4000, 5000)

# construct lag1_wave2 to store the results
lag1_wave2 <- vector("list",2)
names(lag1_wave2) <- c("coef", "non_hete")

lag1_wave2$coef <- vector("list",3)
names(lag1_wave2$coef) <- c("bxy0.2 byx0.2", "bxy0 byx0.2", "bxy0.2 byx0")

lag1_wave2$coef$`bxy0.2 byx0.2` <- 
  lag1_wave2$coef$`bxy0 byx0.2` <- 
  lag1_wave2$coef$`bxy0.2 byx0` <- vector("list",length(n_list))
names(lag1_wave2$coef$`bxy0.2 byx0.2`) <- 
  names(lag1_wave2$coef$`bxy0 byx0.2`) <- 
  names(lag1_wave2$coef$`bxy0.2 byx0`) <- n_list

lag1_wave2$non_hete <- vector("list",3)
names(lag1_wave2$non_hete) <- c("non", "hete", "non_hete")
lag1_wave2$non_hete$non <- 
  lag1_wave2$non_hete$hete <- 
  lag1_wave2$non_hete$non_hete <- vector("list",length(n_list))
names(lag1_wave2$non_hete$non) <- 
  names(lag1_wave2$non_hete$hete) <- 
  names(lag1_wave2$non_hete$non_hete) <- n_list

for (i in 1: length(n_list)) {
  # linear, homogeneous
  lag1_wave2$coef$`bxy0.2 byx0.2`[[i]] <- foreach(1:1000,.combine='rbind') %dorng%
    fun_lag1_wave2(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x)x, heterogeneity = F)
  lag1_wave2$coef$`bxy0 byx0.2`[[i]] <- foreach(1:1000,.combine='rbind') %dorng%
    fun_lag1_wave2(bxy=0, byx=0.2, n=n_list[i], gL=function(x)x, heterogeneity = F)
  lag1_wave2$coef$`bxy0.2 byx0`[[i]] <- foreach(1:1000,.combine='rbind') %dorng%
    fun_lag1_wave2(bxy=0.2, byx=0, n=n_list[i], gL=function(x)x, heterogeneity = F)
  # non-linear
  lag1_wave2$non_hete$non[[i]] <- foreach(1:1000,.combine='rbind', .errorhandling = "remove") %dorng%
    fun_lag1_wave2(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x){x+2*x^2+ifelse(x<0,1,0)}, heterogeneity = F)
  # heterogeneous
  lag1_wave2$non_hete$hete[[i]] <- foreach(1:1000,.combine='rbind', .errorhandling = "remove") %dorng%
    fun_lag1_wave2(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x)x, heterogeneity = T)
  # both non-linear and heterogeneous
  lag1_wave2$non_hete$non_hete[[i]] <- foreach(1:1000,.combine='rbind', .errorhandling = "remove") %dorng%
    fun_lag1_wave2(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x){x+2*x^2+ifelse(x<0,1,0)}, heterogeneity = T)
}

# calculate the true value of estimand
## When heterogeneity exist, the causal effect will not be 0.2, construct a 
## function to estimate the causal effect 
cal_estimand <- function(gL){
  bxy=0.2
  byx=0.2
  n = 1000000
  C0 = rnorm(n, 0, 1)
  L0 = rnorm(n, 0, 1)+0.1*C0
  gL0= do.call(gL, list(L0))
  X0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
  Y0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
  
  L1 = 0.1*C0 + 0.5*L0 + rnorm(n, 0, 1)
  gL1= do.call(gL, list(L1))
  X1 = 0.1*C0 + 0.5*X0 + byx*Y0 + 0.1*gL1 + 0.3*Y0*gL1 + rnorm(n, 0, 1)
  Y1 = 0.1*C0 + 0.5*Y0 + bxy*X0 + 0.1*gL1 + 0.3*X0*gL1 + rnorm(n, 0, 1)
  
  L2 = 0.1*C0 + 0.5*L1 + rnorm(n, 0, 1)
  gL2= do.call(gL, list(L2))
  X2 = 0.1*C0 + 0.5*X1 + byx*Y1 + 0.1*gL2 + 0.3*Y1*gL2 + rnorm(n, 0, 1)
  Y2 = 0.1*C0 + 0.5*Y1 + bxy*X1 + 0.1*gL2 + 0.3*X1*gL2 + rnorm(n, 0, 1)
  
  cor(data.frame(C0,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2))
  summary(data.frame(C0,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2))
  
  X1= rep(1, n)
  Y2_potential_X1 = mean(0.1*C0 + 0.5*Y1 + bxy*X1 + 0.1*gL2 + 0.3*X1*gL2)
  X1= rep(0, n)
  Y2_potential_X0 = mean(0.1*C0 + 0.5*Y1 + bxy*X1 + 0.1*gL2 + 0.3*X1*gL2)
  estimand = Y2_potential_X1-Y2_potential_X0
  return(estimand)
}

# heterogeneous
estimand1 <- cal_estimand(gL=function(x){x})
# both non-linear and heterogeneous
estimand2 <- cal_estimand(gL=function(x){x+2*x^2+ifelse(x<0,1,0)})

# performance evaluation
lag1_wave_true <-lag1_wave2_mean <-lag1_wave2_bias <- lag1_wave2_biasP <- 
  lag1_wave2_var <- lag1_wave2_mse <- lag1_wave2_CIcove <- lag1_wave2_CI <- lag1_wave2

for (i in 1:length(lag1_wave2)) {
  for (j in 1:length(lag1_wave2[[i]])) {
    for (k in 1:length(lag1_wave2[[i]][[j]])) {
      lag1_wave_true[[i]][[j]][[k]] <- rep(0.2 ,4)
    }
  }
}  

for (i in 1:length(n_list)){
  lag1_wave_true$coef$`bxy0 byx0.2`[[i]] <- c(0, 0.2, 0, 0.2)
  lag1_wave_true$coef$`bxy0.2 byx0`[[i]] <- c(0.2, 0, 0.2, 0)
  lag1_wave_true$non_hete$hete[[i]] <- rep(estimand1, 4)
  lag1_wave_true$non_hete$non_hete[[i]] <- rep(estimand2, 4)
}

## calculate bias\variance\mse
for (i in 1:length(lag1_wave2)) {
  for (j in 1:length(lag1_wave2[[i]])) {
    for (k in 1:length(lag1_wave2[[i]][[j]])) {
      lag1_wave2_mean[[i]][[j]][[k]] <- colMeans(lag1_wave2[[i]][[j]][[k]][ ,1:(ncol(lag1_wave2[[i]][[j]][[k]])/2)], na.rm = T)
      lag1_wave2_bias[[i]][[j]][[k]] <- lag1_wave2_mean[[i]][[j]][[k]]-lag1_wave_true[[i]][[j]][[k]]
      lag1_wave2_var[[i]][[j]][[k]] <- colVars(lag1_wave2[[i]][[j]][[k]][ ,1:(ncol(lag1_wave2[[i]][[j]][[k]])/2)], na.rm = T)
      lag1_wave2_mse[[i]][[j]][[k]] <- (lag1_wave2_bias[[i]][[j]][[k]])^2 + lag1_wave2_var[[i]][[j]][[k]]
    }
  }
}

## calculate CI coverage rate
for (i in 1:length(lag1_wave2)) {
  for (j in 1:length(lag1_wave2[[i]])) {
    for (k in 1:length(lag1_wave2[[i]][[j]])){
      for (l in 1:(ncol(lag1_wave2[[i]][[j]][[k]])/2)) {
        L <- lag1_wave2[[i]][[j]][[k]][ ,l] - 1.96*lag1_wave2[[i]][[j]][[k]][ ,ncol(lag1_wave2[[i]][[j]][[k]])/2+l]
        H <- lag1_wave2[[i]][[j]][[k]][ ,l] + 1.96*lag1_wave2[[i]][[j]][[k]][ ,ncol(lag1_wave2[[i]][[j]][[k]])/2+l]
        lag1_wave2_CIcove[[i]][[j]][[k]][ ,l] <- lag1_wave_true[[i]][[j]][[k]][l]>=L & lag1_wave_true[[i]][[j]][[k]][l]<=H
      }
      lag1_wave2_CIcove[[i]][[j]][[k]] <- lag1_wave2_CIcove[[i]][[j]][[k]][ ,1:(ncol(lag1_wave2[[i]][[j]][[k]])/2)]
    }
  }
}

for (i in 1:length(lag1_wave2)) {
  for (j in 1:length(lag1_wave2[[i]])) {
    for (k in 1:length(lag1_wave2[[i]][[j]])) {
      lag1_wave2_CI[[i]][[j]][[k]] <- colMeans(lag1_wave2_CIcove[[i]][[j]][[k]])
    }
  }
}

# lag1 wave3 -------------------------------------------------------------
#simple longitudinal setting with 1-period lagged effect and 3 waves of data
#to demonstrate the impact of “GEE bias” under different correlation structures

fun_lag1_wave3 <- function(bxy, byx, gL, n, corstr){
  # data generation
  n=n
  C0 = rnorm(n, 0, 1)
  L0 = rnorm(n, 0, 1)+0.1*C0
  gL0= do.call(gL, list(L0))
  X0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)
  Y0 = 0.1*C0 + 0.2*gL0 + rnorm(n, 0, 1)

  L1 = 0.1*C0 + 0.5*L0 + rnorm(n, 0, 1)
  gL1= do.call(gL, list(L1))
  X1 = 0.1*C0 + 0.5*X0 + byx*Y0 + 0.2*gL1 + rnorm(n, 0, 1)
  Y1 = 0.1*C0 + 0.5*Y0 + bxy*X0 + 0.2*gL1 + rnorm(n, 0, 1)
  
  L2 = 0.1*C0 + 0.5*L1 + rnorm(n, 0, 1)
  gL2= do.call(gL, list(L2))
  X2 = 0.1*C0 + 0.5*X1 + byx*Y1 + 0.2*gL2 + rnorm(n, 0, 1)
  Y2 = 0.1*C0 + 0.5*Y1 + bxy*X1 + 0.2*gL2 + rnorm(n, 0, 1)
  
  L3 = 0.1*C0 + 0.5*L2 + rnorm(n, 0, 1)
  gL3= do.call(gL, list(L3))
  X3 = 0.1*C0 + 0.5*X2 + byx*Y2 + 0.2*gL3 + rnorm(n, 0, 1)
  Y3 = 0.1*C0 + 0.5*Y2 + bxy*X2 + 0.2*gL3 + rnorm(n, 0, 1)
  
  data <- data.frame(C0,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2)  
  cor(data)
  
  # CLPM
  res_data <- data.frame(X1=residuals(lm(X1 ~ C0 +L1)),
                         Y1=residuals(lm(Y1 ~ C0 +L1)),
                         X2=residuals(lm(X2 ~ C0 +L2)),
                         Y2=residuals(lm(Y2 ~ C0 +L2)),
                         X3=residuals(lm(X3 ~ C0 +L3)),
                         Y3=residuals(lm(Y3 ~ C0 +L3))) 
  model <- 
    "Y3 ~ a1*Y2 + b1*X2
     Y2 ~ a1*Y1 + b1*X1
     X3 ~ a2*X2 + b2*Y2
     X2 ~ a2*X1 + b2*Y1
     X1 ~~r1*Y1
     X2 ~~0*Y2
     X3 ~~0*Y3
     X1 ~~X1
     X2 ~~X2
     X3 ~~X3
     Y1 ~~Y1
     Y2 ~~Y2
     Y3 ~~Y3"
  fit <- sem(model, res_data)
  CLPM_result <- summary(fit, fit.measures = T, 
                         standardized = TRUE,
                         ci = TRUE,
                         rsquare = TRUE)
  CLPM <- CLPM_result[["pe"]][c(2,6), c("est", "se")]
  
  # G-CLPM
  fit_XY1 <- weightit(X1 ~ Y1 + L1 + C0, method = Wmethod, 
                        SL.library = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.earth"))
  fit_XY2 <- weightit(X2 ~ Y2 + L2 + C0, method = Wmethod, 
                        SL.library = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.earth"))
  fit_YX1 <- weightit(Y1 ~ X1 + L1 + C0, method = Wmethod, 
                        SL.library = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.earth"))
  fit_YX2 <- weightit(Y2 ~ X2 + L2 + C0, method = Wmethod, 
                        SL.library = c("SL.glm", "SL.gam", "SL.glm.interaction", "SL.earth"))
  
  dataXY1 <- data.frame(id = 1:n, X=X1, Y=Y2, weight= fit_XY1[["weights"]], wave= rep(1, n))
  dataXY2 <- data.frame(id = 1:n, X=X2, Y=Y3, weight= fit_XY2[["weights"]], wave= rep(2, n))
  dataXY <- rbind(dataXY1, dataXY2)
  dataXY <- dataXY[order(dataXY$id, dataXY$wave), ]
  
  dataYX1 <- data.frame(id = 1:n, Y=Y1, X=X2, weight= fit_YX1[["weights"]], wave= rep(1, n))
  dataYX2 <- data.frame(id = 1:n, Y=Y2, X=X3, weight= fit_YX2[["weights"]], wave= rep(2, n))
  dataYX <- rbind(dataYX1, dataYX2)
  dataYX <- dataYX[order(dataYX$id, dataYX$wave), ]
  
  G_CLPM_XY <- summary(geeglm(Y ~ X, id = id, family = gaussian, data = dataXY,
                              weights = weight, corstr = corstr))
  
  G_CLPM_YX <- summary(geeglm(X ~ Y, id = id, family = gaussian, data = dataYX,
                              weights = weight, corstr = corstr))
  
  G_CLPM <- rbind(G_CLPM_XY$coefficients[2, c("Estimate","Std.err")], 
                  G_CLPM_YX$coefficients[2, c("Estimate","Std.err")])
  colnames(G_CLPM) <- c("est", "se")
  
  result <- matrix(c(as.matrix(rbind(CLPM, G_CLPM))), nrow = 1)
  colnames(result) <- c("CLPM_est_xy","CLPM_est_yx","GCLPM_est_xy","GCLPM_est_yx",
                        "CLPM_se_xy","CLPM_se_yx","GCLPM_se_xy","GCLPM_se_yx")
  
  return(result)
}

clus <- makeCluster(20)
registerDoParallel(clus)

# different sample size
n_list <- c(200, 400, 600, 800, 1000, 1200, 1600, 2000, 2500, 3000, 3500, 4000, 5000)

# construct lag1_wave2 to store the results
lag1_wave3 <- vector("list",3)
names(lag1_wave3) <- c("independence", "exchangeable", "ar1")
lag1_wave3$independence <- 
  lag1_wave3$exchangeable <- 
  lag1_wave3$ar1 <- vector("list",length(n_list))
names(lag1_wave3$independence) <- 
  names(lag1_wave3$exchangeable) <- 
  names(lag1_wave3$ar1) <- n_list

for (i in 1:length(n_list)) {
  # different correlation structure
  lag1_wave3$independence[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag1_wave3(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x)x, corstr="independence")
  lag1_wave3$exchangeable[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag1_wave3(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x)x, corstr="exchangeable")
  lag1_wave3$ar1[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag1_wave3(bxy=0.2, byx=0.2, n=n_list[i], gL=function(x)x, corstr="ar1")
}

# performance evaluation
lag1_wave3_mean <- lag1_wave3_bias <- lag1_wave3_biasP <- 
  lag1_wave3_var <- lag1_wave3_mse <- lag1_wave3_CIcove <- lag1_wave3_CI <- lag1_wave3

## calculate bias\variance\MSE
for (i in 1:length(lag1_wave3)) {
  for (j in 1:length(lag1_wave3[[i]])) {
    lag1_wave3_mean[[i]][[j]] <- colMeans(lag1_wave3[[i]][[j]][ ,1:(ncol(lag1_wave3[[i]][[j]])/2)])
    lag1_wave3_bias[[i]][[j]] <- lag1_wave3_mean[[i]][[j]]-0.2
    lag1_wave3_var[[i]][[j]] <- colVars(lag1_wave3[[i]][[j]][ ,1:(ncol(lag1_wave3[[i]][[j]])/2)])
    lag1_wave3_mse[[i]][[j]] <- (lag1_wave3_bias[[i]][[j]])^2 + lag1_wave3_var[[i]][[j]]
  }
}

## calculate CI coverage
for (i in 1:length(lag1_wave3)) {
  for (j in 1:length(lag1_wave3[[i]])) {
    for (l in 1:(ncol(lag1_wave3[[i]][[j]])/2)) {
      L <- lag1_wave3[[i]][[j]][ ,l] - 1.96*lag1_wave3[[i]][[j]][ ,ncol(lag1_wave3[[i]][[j]])/2+l]
      H <- lag1_wave3[[i]][[j]][ ,l] + 1.96*lag1_wave3[[i]][[j]][ ,ncol(lag1_wave3[[i]][[j]])/2+l]
      lag1_wave3_CIcove[[i]][[j]][,l] <- 0.2>=L & 0.2<= H
    }
    lag1_wave3_CIcove[[i]][[j]] <- lag1_wave3_CIcove[[i]][[j]][ ,1:(ncol(lag1_wave3[[i]][[j]])/2)]
  }
}

for (i in 1:length(lag1_wave3)) {
  for (j in 1:length(lag1_wave3[[i]])) {
    lag1_wave3_CI[[i]][[j]] <- colMeans(lag1_wave3_CIcove[[i]][[j]])
  }
}

# lag2_wave4 --------------------------------------------------------------
#complex longitudinal setting with 2-period lagged effect and 4 waves of data
#identifying the cross-lagged effects under linear and homogeneous data generation scenarios

fun_lag2_wave4 <- function(bxy1, byx1, bxy2, byx2, n, gL, Wmethod, min){
  # data generation
  n=n
  C0 = rnorm(n, 0, 1)
  L_1 = rnorm(n, 0, 1)+0.1*C0
  gL_1= do.call(gL, list(L_1))
  X_1 = 0.1*C0 + 0.2*gL_1 + rnorm(n, 0, 1)
  Y_1 = 0.1*C0 + 0.2*gL_1 + rnorm(n, 0, 1)
  
  L0 = 0.1*C0 + 0.3*L_1 + 0.2*X_1 + 0.2*Y_1 +rnorm(n, 0, 1)
  gL0= do.call(gL, list(L0))
  X0 = 0.1*C0 + 0.3*X_1 + byx1*Y_1 + 0.2*gL0 + rnorm(n, 0, 1)
  Y0 = 0.1*C0 + 0.3*Y_1 + bxy1*X_1 + 0.2*gL0 + rnorm(n, 0, 1)
  
  L1 = 0.1*C0 + 0.3*L0 + 0.2*X0 + 0.2*Y0 + rnorm(n, 0, 1)
  gL1= do.call(gL, list(L1))
  X1 = 0.1*C0 + 0.3*X0 + 0.1*X_1 + byx1*Y0 + byx2*Y_1 + 0.2*gL1 + rnorm(n, 0, 1)
  Y1 = 0.1*C0 + 0.3*Y0 + 0.1*Y_1 + bxy1*X0 + bxy2*X_1 + 0.2*gL1 + rnorm(n, 0, 1)
  
  L2 = 0.1*C0 + 0.3*L1 + 0.2*X1 + 0.2*Y1 + rnorm(n, 0, 1)
  gL2= do.call(gL, list(L2))
  X2 = 0.1*C0 + 0.3*X1 + 0.1*X0 + byx1*Y1 + byx2*Y0 + 0.2*gL2 + rnorm(n, 0, 1)
  Y2 = 0.1*C0 + 0.3*Y1 + 0.1*Y0 + bxy1*X1 + bxy2*X0 + 0.2*gL2 + rnorm(n, 0, 1)
  
  L3 = 0.1*C0 + 0.3*L2 + 0.2*X2 + 0.2*Y2 + rnorm(n, 0, 1)
  gL3= do.call(gL, list(L3))
  X3 = 0.1*C0 + 0.3*X2 + 0.1*X1 + byx1*Y2 + byx2*Y1 + 0.2*gL3 + rnorm(n, 0, 1)
  Y3 = 0.1*C0 + 0.3*Y2 + 0.1*Y1 + bxy1*X2 + bxy2*X1 + 0.2*gL3 + rnorm(n, 0, 1)
  
  data <- data.frame(C0,L_1,X_1,Y_1,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2,L3,X3,Y3)
  cor(data)
  
  # CLPM
  res_data <- data.frame(X1=residuals(lm(X1 ~ C0 +L1)),
                         Y1=residuals(lm(Y1 ~ C0 +L1)),
                         X2=residuals(lm(X2 ~ C0 +L2)),
                         Y2=residuals(lm(Y2 ~ C0 +L2)),
                         X3=residuals(lm(X3 ~ C0 +L3)),
                         Y3=residuals(lm(Y3 ~ C0 +L3))) 
  model <- 
    "Y3 ~ a1*Y2 + b1*X2 + c1*Y1 + d1*X1
     Y2 ~ a1*Y1 + b1*X1
     X3 ~ a2*X2 + b2*Y2 + c2*X1 + d2*Y1
     X2 ~ a2*X1 + b2*Y1
     X1 ~~r1*Y1
     X2 ~~r2*Y2
     X3 ~~0*Y3
     X1 ~~X1
     X2 ~~X2
     X3 ~~X3
     Y1 ~~Y1
     Y2 ~~Y2
     Y3 ~~Y3"
  fit <- sem(model, res_data)
  CLPM_result <- summary(fit, fit.measures = T, 
                         standardized = TRUE,
                         ci = TRUE,
                         rsquare = TRUE)
  # b1 d1 b2 d2
  CLPM <- CLPM_result[["pe"]][c(2,4,8,10), c("est", "se")]
  
  # G-CLPM
  if(min ==T){
    fit_XY <- weightitMSM(list(X1 ~ X0 + Y1 + Y0 + L1 + C0,
                               X2 ~ X1 + Y2 + Y1 + L2 + C0), method = Wmethod, stabilize = T, data = data)
    fit_YX <- weightitMSM(list(Y1 ~ Y0 + X1 + X0 + L1 + C0,
                               Y2 ~ Y1 + X2 + X1 + L2 + C0), method = Wmethod, stabilize = T, data = data)
  }else{
    fit_XY <- weightitMSM(list(X1 ~ X0 + X_1 + Y0 + Y_1 + L1 + C0,
                               X2 ~ X1 + X0 + Y1 + Y0 + L2 + C0), method = Wmethod, stabilize = T, data = data)
    fit_YX <- weightitMSM(list(Y1 ~ Y0 + Y_1 + X0 + X_1 + L1 + C0,
                               Y2 ~ Y1 + Y0 + X1 + X0 + L2 + C0), method = Wmethod, stabilize = T, data = data)
  }
  
  bal.tab(fit_XY)
  bal.tab(fit_YX)
  
  G_CLPM_XY <- summary(geeglm(Y3 ~ X2 + X1, weights = fit_XY[["weights"]], id= rownames(data)))
  G_CLPM_YX <- summary(geeglm(X3 ~ Y2 + Y1, weights = fit_YX[["weights"]], id= rownames(data)))
  
  G_CLPM <- rbind(G_CLPM_XY$coefficients[c(2,3), c("Estimate","Std.err")], 
                  G_CLPM_YX$coefficients[c(2,3), c("Estimate","Std.err")])
  
  colnames(G_CLPM) <- c("est", "se")
  
  result <- matrix(c(as.matrix(rbind(CLPM, G_CLPM))), nrow = 1)
  
  colnames(result) <- c("CLPM_est_xy1","CLPM_est_xy2","CLPM_est_yx1","CLPM_est_yx2",
                        "GCLPM_est_xy1","GCLPM_est_xy2","GCLPM_est_yx1","GCLPM_est_yx2",
                        "CLPM_se_xy1","CLPM_se_xy2","CLPM_se_yx1","CLPM_se_yx2",
                        "GCLPM_se_xy1","GCLPM_se_xy2","GCLPM_se_yx1","GCLPM_se_yx")
  
  return(result)
}

clus <- makeCluster(20)
registerDoParallel(clus)

#### different sample size
n_list <- c(200, 400, 600, 800, 1000, 1200, 1600, 2000, 2500, 3000, 3500, 4000, 5000)

#### construct lag2_wave4 to store the results
lag2_wave4 <- vector("list",3)
names(lag2_wave4) <- c("bxy1_0.2 byx1_0.2", "bxy1_0 byx1_0.2", "bxy1_0.2 byx1_0")
lag2_wave4$`bxy1_0.2 byx1_0.2` <- 
  lag2_wave4$`bxy1_0 byx1_0.2` <- 
  lag2_wave4$`bxy1_0.2 byx1_0` <- vector("list",length(n_list))
names(lag2_wave4$`bxy1_0.2 byx1_0.2`) <- 
  names(lag2_wave4$`bxy1_0 byx1_0.2`) <- 
  names(lag2_wave4$`bxy1_0.2 byx1_0`) <- n_list

for (i in 1: length(n_list)){
  lag2_wave4$`bxy1_0.2 byx1_0.2`[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag2_wave4(bxy1=0.2, byx1=0.2, bxy2=0.1, byx2=0.1, n=n_list[i], gL=function(x)x, Wmethod="glm", min = T)
  lag2_wave4$`bxy1_0 byx1_0.2`[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag2_wave4(bxy1=0, byx1=0.2, bxy2=0, byx2=0.1, n=n_list[i], gL=function(x)x, Wmethod="glm", min = T)
  lag2_wave4$`bxy1_0.2 byx1_0`[[i]] <- foreach(1:1000,.combine='rbind') %dorng% 
    fun_lag2_wave4(bxy1=0.2, byx1=0, bxy2=0.1, byx2=0, n=n_list[i], gL=function(x)x, Wmethod="glm", min = T)
}

# calculate the true value of estimand
cal_estimand <- function(bxy1, byx1, bxy2, byx2){
  n = 1000000
  gL=function(x)x
  Wmethod="glm"
  C0 = rnorm(n, 0, 1)
  L_1 = rnorm(n, 0, 1)+0.1*C0
  gL_1= do.call(gL, list(L_1))
  X_1 = 0.1*C0 + 0.2*gL_1 + rnorm(n, 0, 1)
  Y_1 = 0.1*C0 + 0.2*gL_1 + rnorm(n, 0, 1)
  
  L0 = 0.1*C0 + 0.3*L_1 + 0.2*X_1 + 0.2*Y_1 +rnorm(n, 0, 1)
  gL0= do.call(gL, list(L0))
  X0 = 0.1*C0 + 0.3*X_1 + byx1*Y_1 + 0.2*gL0 + rnorm(n, 0, 1)
  Y0 = 0.1*C0 + 0.3*Y_1 + bxy1*X_1 + 0.2*gL0 + rnorm(n, 0, 1)
  
  L1 = 0.1*C0 + 0.3*L0 + 0.2*X0 + 0.2*Y0 + rnorm(n, 0, 1)
  gL1= do.call(gL, list(L1))
  X1 = 0.1*C0 + 0.3*X0 + 0.1*X_1 + byx1*Y0 + byx2*Y_1 + 0.2*gL1 + rnorm(n, 0, 1)
  Y1 = 0.1*C0 + 0.3*Y0 + 0.1*Y_1 + bxy1*X0 + bxy2*X_1 + 0.2*gL1 + rnorm(n, 0, 1)
  
  L2 = 0.1*C0 + 0.3*L1 + 0.2*X1 + 0.2*Y1 + rnorm(n, 0, 1)
  gL2= do.call(gL, list(L2))
  X2 = 0.1*C0 + 0.3*X1 + 0.1*X0 + byx1*Y1 + byx2*Y0 + 0.2*gL2 + rnorm(n, 0, 1)
  Y2 = 0.1*C0 + 0.3*Y1 + 0.1*Y0 + bxy1*X1 + bxy2*X0 + 0.2*gL2 + rnorm(n, 0, 1)
  
  L3 = 0.1*C0 + 0.3*L2 + 0.2*X2 + 0.2*Y2 + rnorm(n, 0, 1)
  gL3= do.call(gL, list(L3))
  X3 = 0.1*C0 + 0.3*X2 + 0.1*X1 + byx1*Y2 + byx2*Y1 + 0.2*gL3 + rnorm(n, 0, 1)
  Y3 = 0.1*C0 + 0.3*Y2 + 0.1*Y1 + bxy1*X2 + bxy2*X1 + 0.2*gL3 + rnorm(n, 0, 1)
  
  cor(data.frame(C0,L_1,X_1,Y_1,L0,X0,Y0,L1,X1,Y1,L2,X2,Y2,L3,X3,Y3))
  
  # X1=1, X2=0
  X1 = rep(1, n)
  X2 = rep(0, n)
  PL2 = 0.1*C0 + 0.3*L1 + 0.2*X1 + 0.2*Y1
  PY2 = 0.1*C0 + 0.3*Y1 + 0.1*Y0 + bxy1*X1 + bxy2*X0 + 0.2*PL2
  PL3 = 0.1*C0 + 0.3*PL2 + 0.2*X2 + 0.2*PY2
  PY3 = 0.1*C0 + 0.3*PY2 + 0.1*Y1 + bxy1*X2 + bxy2*X1 + 0.2*PL3
  Y3_potential_10 = mean(PY3)
  
  # X1=0, X2=1
  X1 = rep(0, n)
  X2 = rep(1, n)
  PL2 = 0.1*C0 + 0.3*L1 + 0.2*X1 + 0.2*Y1
  PY2 = 0.1*C0 + 0.3*Y1 + 0.1*Y0 + bxy1*X1 + bxy2*X0 + 0.2*PL2
  PL3 = 0.1*C0 + 0.3*PL2 + 0.2*X2 + 0.2*PY2
  PY3 = 0.1*C0 + 0.3*PY2 + 0.1*Y1 + bxy1*X2 + bxy2*X1 + 0.2*PL3
  Y3_potential_01 = mean(PY3)
  
  # X1=0, X2=0
  X1 = rep(0, n)
  X2 = rep(0, n)
  PL2 = 0.1*C0 + 0.3*L1 + 0.2*X1 + 0.2*Y1
  PY2 = 0.1*C0 + 0.3*Y1 + 0.1*Y0 + bxy1*X1 + bxy2*X0 + 0.2*PL2
  PL3 = 0.1*C0 + 0.3*PL2 + 0.2*X2 + 0.2*PY2
  PY3 = 0.1*C0 + 0.3*PY2 + 0.1*Y1 + bxy1*X2 + bxy2*X1 + 0.2*PL3
  Y3_potential_00 = mean(PY3)
  
  estimand1 = Y3_potential_01-Y3_potential_00
  estimand2 = Y3_potential_10-Y3_potential_00
  
  return(c(estimand1, estimand2))
}
estimand_0.2_0.1 <- cal_estimand(bxy1=0.2,byx1=0.2,bxy2=0.1,byx2=0.1)
estimand_0_0 <- cal_estimand(bxy1=0, byx1=0.2, bxy2=0, byx2=0.1)

for (i in 1:length(n_list)) {
  lag2_wave4_true$`bxy1_0.2 byx1_0.2`[[i]] <- c(estimand_0.2_0.1, estimand_0.2_0.1,
                                                estimand_0.2_0.1, estimand_0.2_0.1)
  
  lag2_wave4_true$`bxy1_0 byx1_0.2`[[i]] <- c(estimand_0_0, estimand_0.2_0.1,
                                              estimand_0_0, estimand_0.2_0.1)
  
  lag2_wave4_true$`bxy1_0.2 byx1_0`[[i]] <- c(estimand_0.2_0.1, estimand_0_0,
                                              estimand_0.2_0.1, estimand_0_0)
}

# performance evaluation
lag2_wave4_true <- lag2_wave4_mean <- lag2_wave4_bias <- lag2_wave4_biasP <- 
  lag2_wave4_var <- lag2_wave4_mse <- lag2_wave4_CIcove <- lag2_wave4_CI <- lag2_wave4

## calculate bias\variance\mse
for (i in 1:length(lag2_wave4)) {
  for (j in 1:length(lag2_wave4[[i]])) {
    lag2_wave4_mean[[i]][[j]] <- colMeans(lag2_wave4[[i]][[j]][ ,1: (ncol(lag2_wave4[[i]][[j]])/2)])
    lag2_wave4_bias[[i]][[j]] <- lag2_wave4_mean[[i]][[j]]-lag2_wave4_true[[i]][[j]]
    lag2_wave4_var[[i]][[j]] <- colVars(lag2_wave4[[i]][[j]][ ,1: (ncol(lag2_wave4[[i]][[j]])/2)])
    lag2_wave4_mse[[i]][[j]] <- (lag2_wave4_bias[[i]][[j]])^2 + lag2_wave4_var[[i]][[j]]
  }
}


## calculate CI coverage rate
for (i in 1:length(lag2_wave4)) {
  for (j in 1:length(lag2_wave4[[i]])) {
    for (l in 1:(ncol(lag2_wave4[[i]][[j]])/2)) {
      L <- lag2_wave4[[i]][[j]][ ,l] - 1.96*lag2_wave4[[i]][[j]][ ,ncol(lag2_wave4[[i]][[j]])/2+l]
      H <- lag2_wave4[[i]][[j]][ ,l] + 1.96*lag2_wave4[[i]][[j]][ ,ncol(lag2_wave4[[i]][[j]])/2+l]
      lag2_wave4_CIcove[[i]][[j]][,l] <- lag2_wave4_true[[i]][[j]][l]>=L & lag2_wave4_true[[i]][[j]][l]<= H
    }
    lag2_wave4_CIcove[[i]][[j]] <- lag2_wave4_CIcove[[i]][[j]][ ,1:(ncol(lag2_wave4[[i]][[j]])/2)]
  }
}

for (i in 1:length(lag2_wave4)) {
  for (j in 1:length(lag2_wave4[[i]])) {
    lag2_wave4_CI[[i]][[j]] <- colMeans(lag2_wave4_CIcove[[i]][[j]])
  }
}