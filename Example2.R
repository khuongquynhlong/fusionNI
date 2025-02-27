library(tidyverse)
library(survey)
library(doParallel)
library(foreach)



#===============================================================================
#------------------------- Example 2: Simulation study ------------------------- 
#===============================================================================

#----- Generate data function
#===============================================================================
sim_data <- function(n1, n2, seed = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  N <- n1 + n2
  X1 <- rnorm(N, 0, 0.5)
  X2 <- rbinom(N, 1, 0.4)
  
  #--- Generate S: 1 = local trial, 0 = distal trial
  # P(S = 1 | X) = plogis(a_est - X1 - X2) such that
  # the expected number of S == 1 equals n1.
  f_a <- function(a) {
    sum(plogis(a - X1 - X2)) - n1
  }
  
  # Find the root
  a_est <- uniroot(f_a, interval = c(-10, 10), tol = 1e-10, maxiter = 10000)$root
  S <- rbinom(N, size = 1, prob = plogis(a_est - X1 - X2))
  
  #--- Generate A: P(A = 1| S = 1)
  # For individuals with S == 1, A is drawn to be either 1 or 2.
  # For individuals with S == 0, A is drawn to be either 0 or 1.
  A <- rep(NA, N)
  A[S == 1] <- rbinom(sum(S == 1), size = 1, prob = 0.5) + 1  # results in {1, 2}
  A[S == 0] <- rbinom(sum(S == 0), size = 1, prob = 0.5)      # results in {0, 1}
  
  #---- Combine dataset
  gen_data <- data.frame(X1 = X1, X2 = X2, S = S, A = A)
  
  #--- Generate potential outcome: Y0, Y1, Y2
  # Y ~ expit(X1 + X2 + 3*I(A = 0) + I(A = 1) + 2*I(A = 2))
  lin_pred0 <- -3 + gen_data$X1 + gen_data$X2 + 2    # for A = 0
  lin_pred1 <- -3 + gen_data$X1 + gen_data$X2 + 1    # for A = 1
  lin_pred2 <- -3 + gen_data$X1 + gen_data$X2 + 1.2  # for A = 2
  
  gen_data$Y0 <- rbinom(n = N, size = 1, prob = plogis(lin_pred0))
  gen_data$Y1 <- rbinom(n = N, size = 1, prob = plogis(lin_pred1))
  gen_data$Y2 <- rbinom(n = N, size = 1, prob = plogis(lin_pred2))
  
  # Observed Y
  # Y = I{A = 0} * Y0 + I{A = 1} * Y1 + I{A = 2} * Y2
  gen_data$Y <- ifelse(gen_data$A == 0, gen_data$Y0,
                       ifelse(gen_data$A == 1, gen_data$Y1, gen_data$Y2))
  
  return(gen_data)
}


#----- Weight generating function
#===============================================================================
gen_weight <- function(data) {
  # W_a = P(A = a|X, S=d)
  distal_data <- subset(data, S == 0)
  a_mod       <- glm(A ~ X1 + X2, family = "binomial", data = distal_data)
  a_pred      <- predict(a_mod, newdata = data, type = "response")
  w_a         <- 1/a_pred

  # W_trans
  s_mod       <- glm(S ~ X1 + X2, family = "binomial", data = data)
  s_pred      <- predict(s_mod, type = "response")
  w_trans     <- s_pred/(1-s_pred)
  w_trial     <-  (1 - mean(data$S))/mean(data$S)
  w_2t        <- w_trans*w_trial

  # All components
  w_all        <- w_a*w_2t
  data$w_trans <- w_trans
  data$w_2t    <- w_2t
  data$w_all   <- w_all

  return(data)
}


# Simplified weights: this is sufficient
gen_weight_s <- function(data) {
  s_mod        <- glm(S ~ X1 + X2, family = "binomial", data = data)
  s_pred       <- predict(s_mod, type = "response")
  w_trial      <- s_pred/(1-s_pred)
  data$w_trans <- w_trans
  
  return(data)
}



#----- Fusion estimator
#===============================================================================

#----- Estimates from local and distal trials separately
trial_est <- function(data){
  dataS0 <- subset(data, S == 0)
  dataS1 <- subset(data, S == 1)
  
  # psi_local
  psi_local1 <- mean(dataS1[dataS1$A == 1, ]$Y)
  psi_local2 <- mean(dataS1[dataS1$A == 2, ]$Y)
  psi_local  <- psi_local1 - psi_local2
  
  # psi_distal
  psi_distal1 <- mean(dataS0[dataS0$A == 1, ]$Y) 
  psi_distal0 <- mean(dataS0[dataS0$A == 0, ]$Y)
  psi_distal  <- psi_distal1 - psi_distal0
  
  est_list <- list(psi_local1  = psi_local1, 
                   psi_local2  = psi_local2, 
                   psi_local   = psi_local,
                   psi_distal0 = psi_distal0,
                   psi_distal1 = psi_distal1,
                   psi_distal  = psi_distal)
  return(est_list)
}



#----- Estimates by transporting distal to local (psi_trans)
trans_est <- function(data){
  dataS0     <- subset(data, S == 0)
  # Applying the weight
  Y1w_mod    <- svyglm(Y ~ A, family = "quasibinomial", 
                    design = svydesign(~ 1, weights = ~dataS0$w_trans, data = dataS0))
  pY1        <- predict(Y1w_mod, newdata = data.frame(A = 1), type = "response") 
  pY0        <- predict(Y1w_mod, newdata = data.frame(A = 0), type = "response") 
  
  psi_trans1 <- mean(pY1)
  psi_trans0 <- mean(pY0)
  psi_trans  <- psi_trans1 - psi_trans0
  
  est_list   <- list(psi_trans0 = psi_trans0,
                     psi_trans1 = psi_trans1, 
                     psi_trans = psi_trans)
  return(est_list)
}




#----- Combine function for simulation example
#===============================================================================
sim_combine <- function(data){
  
  dataS1      <- subset(data, S == 1)
  
  # True values
  true_Y2     <- mean(dataS1$Y2)
  true_Y0     <- mean(dataS1$Y0)
  true_fusion <- true_Y2 - true_Y0
  
  # Estimates
  trial_list  <- trial_est(data)
  trans_list  <- trans_est(data)
  
  # Fusion estimate
  psi_fusion  <- trans_list$psi_trans - trial_list$psi_local
 
  # Crude estimate 
  psi_crude   <- trial_list$psi_local2 - trial_list$psi_distal0
  
  # Difference between distal and trans for active control (A = 1)
  C_diff      <- trans_list$psi_trans1 - trial_list$psi_distal1
  
  # Bias: absolute
  bias_fusion <- true_fusion - psi_fusion
  bias_crude  <- true_fusion - psi_crude
  
  est_list    <- data.frame(true_fusion = true_fusion,
                            psi_fusion  = psi_fusion, 
                            psi_crude   = psi_crude,
                            C_diff      = C_diff,
                            bias_fusion = bias_fusion,
                            bias_crude  = bias_crude)
  return(est_list)
}




num_cores <- parallel::detectCores() - 4

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define simulation grid parameters
n1_vals    <- c(200, 500, 1000, 2000)
n2_vals    <- c(200, 500, 1000, 2000)
iterations <- 100000

# Create a grid of parameter combinations
param_grid <- expand.grid(n1 = n1_vals, n2 = n2_vals)

# List to store simulation results
results_list <- vector("list", nrow(param_grid))


start_time <- Sys.time()

set.seed(12345) 

# Loop over each combination of n1 and n2
for (j in 1:nrow(param_grid)) {
  n1 <- param_grid$n1[j]
  n2 <- param_grid$n2[j]
  
  # Parallel loop for iterations for current parameter combination
  sim_results <- foreach(i = 1:iterations, .combine = rbind, 
                         .packages = c("survey")) %dopar% {
                           df <- sim_data(n1 = n1, n2 = n2) |> gen_weight()
                           # Compute estimates 
                           sim_result <- sim_combine(df)
                           # Append simulation parameters
                           sim_result$n1   <- n1
                           sim_result$n2   <- n2
                           sim_result$iter <- i
                           # Return
                           sim_result 
                         }
  
  results_list[[j]] <- sim_results
}

result_tab <- do.call(rbind, results_list)

stopCluster(cl)

# Runtime
end_time <- Sys.time()
runtime <- end_time - start_time
runtime


# saveRDS(result_tab, "Example2.RDS")


result_tab <- readRDS("Example2.RDS")






# Probability scale
result_sum <- result_tab |>
  group_by(n1, n2) |>
  summarise(
    true_fusion   = mean(true_fusion),
    psi_fusion    = mean(psi_fusion),
    psi_crude     = mean(psi_crude),
    C_diff        = mean(C_diff),
    mbias_fusion  = mean(bias_fusion),
    mbias_crude   = mean(bias_crude),
    rmse_fusion   = sqrt(mean(bias_fusion^2)),
    rmse_crude    = sqrt(mean(bias_crude^2)),
    fusion_mce    = sd(bias_fusion),
    crude_mce     = sd(bias_crude)
  ) |>
  mutate_at(c("mbias_fusion", "mbias_crude", "rmse_fusion",
              "rmse_crude", "fusion_mce", "crude_mce"), 
            function(x){x = round(x, 4)})

writexl::write_xlsx(result_sum, "result_sum.xlsx")



# Percentage scale
result_per <- result_tab |>
  mutate(
    bias_fusion   = bias_fusion*100,
    bias_crude    = bias_crude*100) |>
  group_by(n1, n2) |>
  summarise(
    true_fusion   = mean(true_fusion*100),
    psi_fusion    = mean(psi_fusion*100),
    psi_crude     = mean(psi_crude*100),
    C_diff        = mean(C_diff*100),
    mbias_fusion  = mean(bias_fusion),
    mbias_crude   = mean(bias_crude),
    rmse_fusion   = sqrt(mean(bias_fusion^2)),
    rmse_crude    = sqrt(mean(bias_crude^2)),
    fusion_mce    = sd(bias_fusion),
    crude_mce     = sd(bias_crude)
  ) |>
  mutate_at(c("mbias_fusion", "mbias_crude", "rmse_fusion",
              "rmse_crude", "fusion_mce", "crude_mce"), 
            function(x){x = round(x, 4)})






