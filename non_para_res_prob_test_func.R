##########################################################################
#### Non-parametric Restricted Prob. test function #######################
##########################################################################


# Require packages 
require(survival)
require(timereg)


# Create function 
nonParamResProb <- function(time_obs, 
                            cens_ind,
                            trt,
                            dataset, 
                            tau, 
                            delta, 
                            alpha = 0.05){
  
  # Change names in data set
  dataset$time_obs <- time_obs
  dataset$cens_ind <- cens_ind
  dataset$trt <- trt
  
  #Fit the model
  fit_aalen <- aalen(Surv(time_obs, cens_ind) ~ trt, data = dataset, 
                     resample.iid = 1)
  
  # Extracting time, and cumhaz, also wrt. tau
  time <- fit_aalen$cum[,1]
  time_wh <- time < tau
  
  trt0_cumhaz <- fit_aalen$cum[,2]
  trt11_cumhaz <- fit_aalen$cum[,3]
  trt1_cumhaz <- fit_aalen$cum[,2] + fit_aalen$cum[,3]
  
  trt0_cumhaz_tau <- fit_aalen$cum[,2][time_wh]
  trt11_cumhaz_tau <- fit_aalen$cum[,3][time_wh]
  trt1_cumhaz_tau <- fit_aalen$cum[,2][time_wh] + fit_aalen$cum[,3][time_wh]
  
  # Survival estimates
  trt0_bresl <- exp(-trt0_cumhaz)
  trt1_bresl <- exp(-trt1_cumhaz)
  
  trt0_bresl_tau <- exp(-trt0_cumhaz_tau)
  trt1_bresl_tau <- exp(-trt1_cumhaz_tau)
  
  # Largest differences, which t
  larg_diff <- max(abs(trt0_bresl_tau - trt1_bresl_tau))
  larg_diff_no <- which(abs(trt0_bresl_tau - trt1_bresl_tau) == larg_diff)[1]
  larg_diff_t <- time[larg_diff_no]
  
  
  # Estimated difference
  estS_that <- trt0_bresl[larg_diff_no] - trt1_bresl[larg_diff_no]
  
  #################################################
  ######### Finding variance estimate #############
  #################################################
  
  
  cumhaziid <- unname(sapply(1:length(fit_aalen$B.iid), function(i) fit_aalen$B.iid[[i]][larg_diff_no,]))

  
  ### So, we need to compute an estimate of V(eps_intercept, eps_trt) at time larg_diff_no
  Veps_int <- sum(cumhaziid[1,]^2)
  Veps_trt <- sum(cumhaziid[2,]^2)
  Veps_int_trt <- sum(cumhaziid[1,] * cumhaziid[2,])
  
  # Collect
  VarMat <- matrix(c(Veps_int, Veps_int_trt,
                     Veps_int_trt, Veps_trt), byrow = T, ncol = 2)
  
  #### FIRST TRANSFORMATION - SURVIVAL FUNCTIONS ######
  
  J_h <- function(x, y){
    matrix(c(-exp(-x), 0,
             -exp(-(x+y)), -exp(-(x+y))), byrow = T, ncol = 2)
  }
  J_h_ins <- J_h(x = trt0_cumhaz[larg_diff_no], y = trt11_cumhaz[larg_diff_no])
  J_h_ins
  
  VarUpt <- J_h_ins %*% VarMat %*% t(J_h_ins)
  VarUpt
  
  ### SECOND TRANSFORMATION - DIFFERENCE in survival functions ###
  # Depends on location of extremal point
  if (trt0_bresl[larg_diff_no] > trt1_bresl[larg_diff_no]){
    J_f <- matrix(c(1,-1), ncol = 2)
  }
  if (trt0_bresl[larg_diff_no] < trt1_bresl[larg_diff_no]){
    J_f <- matrix(c(-1, 1), ncol = 2)
  }
  #
  VarUpt2 <- J_f %*% VarUpt %*% t(J_f)
  VarUpt2
  
  
  # This is the asymSD
  asymSD <- sqrt(VarUpt2)
  
  
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
