##########################################################################
#### Cox-based Restricted Prob. test function ############################
##########################################################################

# Require packages 
require(survival)
require(timereg)


# Create function
CoxBasedResProb <- function(time_obs, 
                            cens_ind,
                            trt,
                            dataset, 
                            tau, 
                            delta, 
                            alpha = 0.05) {
  
  
  # Change names in data set
  dataset$time_obs <- time_obs
  dataset$cens_ind <- cens_ind
  dataset$trt <- trt
  
  # Fit model 
  fitCA <- cox.aalen(Surv(time_obs, cens_ind) ~ prop(trt), data = dataset,  
                     resample.iid = T, 
                     max.timepoint.sim = NULL)
  
  
  
  # Finding tau again
  # note that cum is cumulative coefficients- but same structure
  tauplac <- which(fitCA$cum[,"time"] == 
                     max(fitCA$cum[, "time"][fitCA$cum[, "time"] <= tau]) )
  
  ## Non-parametric part of model resampled : fitCA$B.iid
  ## So if I want the tau time point, I need to look across the 1000 simulations.
  cumhaziid <- unname(sapply(1:length(fitCA$B.iid), function(i) fitCA$B.iid[[i]][tauplac,]))
  
  
  ## Parametric part of model resampled : fitCA$gamma.iid
  betaiid <- fitCA$gamma.iid
  
  ### So, we need to compute an estimate of V(eps_1, eps_2)
  ## We know that we can find an consistent estimator of this (see LaTeX)
  Veps1 <- sum(betaiid^2)
  Veps2 <- sum(cumhaziid^2)
  Veps1eps2 <- sum(betaiid * cumhaziid)
  
  VarMat <- matrix(c(Veps1, Veps1eps2, 
                     Veps1eps2, Veps2), byrow = T, ncol = 2)
  ####### This is the variance-covariance matrix of (beta, Lambda0(tau))
  
  
  
  ################### TRANSFORMATIONS ##################### 
  ### Step 1: Finding an estimate of Omega
  
  JacobianFunc <- function(beta, cumhaz){
    matrix(c(exp(beta), 0, 
             exp(beta)^(-1/(exp(beta)-1)) * ((exp(beta) * beta - exp(beta) + 1) / ((exp(beta)-1)^2)), 
             0, 
             0, -exp(-cumhaz)), 
           byrow = T, ncol = 2)
  }
  
  Sigma <- VarMat 
  
  beta <- unname(fitCA$gamma)
  cumhaz <- unname(fitCA$cum[tauplac,][2])
  
  J_f <- JacobianFunc(beta = beta, cumhaz = cumhaz)
  
  Omega <- J_f %*% Sigma %*% t(J_f)
  Omega
  ####### This is the variance-covariance matrix of (theta, theta^{-1/(theta -1)}, S0(tau))
  
  
  ##### STEP 2
  
  # Transformations
  theta <- exp(beta)
  S0tau <- exp(-cumhaz)
  thetatrans <- theta^{-1/(theta-1)}
  
  J_g1_func <- function(x, y, z){
    matrix(c(1, 0, 0, 
             0, 1, 0, 
             0, 0, 1, 
             0, 2*(y-z), -2*(y-z)), 
           byrow = T, ncol = 3)
  }
  
  J_g1 <- J_g1_func(x = theta, y = thetatrans, z = S0tau)
  
  up <- J_g1 %*% Omega %*% t(J_g1)
  up
  
  ##### STEP 3
  J_g2_func <- function(x, y, z, w){
    matrix(c(1, 0, 0, 0, 
             0, 1, 0, 0, 
             0, 0, 1, 0, 
             0, 0, 0, 1 / (2 * sqrt(w))), byrow = T, ncol = 4)
  }
  
  J_g2 <- J_g2_func(x = theta, y = thetatrans, z = S0tau, w = (thetatrans - S0tau)^2)
  
  
  up2 <- J_g2 %*% up %*% t(J_g2)
  up2
  
  
  ##### STEP 4
  J_g3_func <- function(){
    matrix(c(1, 0, 0, 0, 
             0, 1/2, 1/2, 1/2), byrow = T, ncol = 4)
  }
  J_g3 <- J_g3_func()
  
  up3 <- J_g3 %*% up2 %*% t(J_g3)
  up3
  
  
  #### STEP 5
  J_k1_func <- function(x, y) {
    J_k <- c(- y^x * log(y), 1 - x * y^(x-1))
    t(as.matrix(J_k))
  }
  
  xi <- max(thetatrans, S0tau)
  # Same as: 1/2 * (thetatrans + S0tau + abs(thetatrans - S0tau))
  
  J_k1 <- J_k1_func(x = theta, y = xi)
  
  up4 <- J_k1 %*% up3 %*% t(J_k1)
  
  
  # This is the asymp SD
  asymSD <- sqrt(up4)
  
  
  # Collecting results
  larg_diff  <- abs(xi - xi^theta)
  
  # Make the test 
  c_a <- qnorm(alpha)
  TestStat <- (larg_diff - delta) / asymSD
  TestRes <- TestStat < c_a
  
  # Display result
  if (TestRes == FALSE) {Result = "Cannot reject H0 (treatment effects are non-equivalent)"}
  if (TestRes == TRUE) {Result = "Can reject H0 (treatment effects non-equivalent)"}
  
  list(asymSD = asymSD, 
       larg_diff = larg_diff, 
       TestStat = TestStat, 
       TestRes = TestRes, 
       delta = delta, 
       quantile = c_a,
       Result = Result)
}




# Check on gastric cancer data
setwd("~/Desktop/Speciale/R article")
gasdata <- read.table("gasdata.txt")
head(gasdata)

# Run the function
res <- CoxBasedResProb(time_obs = gasdata$time, 
                       cens_ind = gasdata$status, 
                       trt = gasdata$trt, 
                       dataset = gasdata, tau = 2000, delta = 0.2, alpha = 0.05)
res