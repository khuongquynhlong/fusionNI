library(tidyverse)
# library(survey)
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
gen_weight     <- function(data) {
  s_mod        <- glm(S ~ X1 + X2, family = "binomial", data = data)
  s_pred       <- predict(s_mod, type = "response")
  w_trans      <- s_pred/(1-s_pred)
  data$w_trans <- w_trans
  
  return(data)
}



#----- Fusion estimator
#===============================================================================

#----- Estimates from local and distal trials
trial_est <- function(data){
  Y12 <- data[data$S == 1 & data$A == 2, ]$Y
  Y11 <- data[data$S == 1 & data$A == 1, ]$Y
  Y01 <- data[data$S == 0 & data$A == 1, ]$Y
  Y00 <- data[data$S == 0 & data$A == 0, ]$Y
  
  # psi_local
  psi_local1 <- mean(Y11)
  psi_local2 <- mean(Y12)
  psi_local  <- psi_local1 - psi_local2
  
  # Binomial se for psi
  se_psi_local1 <- sqrt((psi_local1 * (1 - psi_local1))/length(Y11))
  se_psi_local2 <- sqrt((psi_local2 * (1 - psi_local2))/length(Y12))
  se_psi_local  <- sqrt(se_psi_local1^2 + se_psi_local2^2)
  
  # psi_distal
  psi_distal1 <- mean(Y01) 
  psi_distal0 <- mean(Y00)
  psi_distal  <- psi_distal1 - psi_distal0
  
  # Binomial se for psi
  se_psi_distal1 <- sqrt((psi_distal1 * (1 - psi_distal1))/length(Y01))
  se_psi_distal0 <- sqrt((psi_distal0 * (1 - psi_distal0))/length(Y00))
  se_psi_distal  <- sqrt(se_psi_distal1^2 + se_psi_distal0^2)
  
  
  est_list <- list(psi_local1     = psi_local1, 
                   psi_local2     = psi_local2, 
                   psi_local      = psi_local,
                   psi_distal0    = psi_distal0,
                   psi_distal1    = psi_distal1,
                   psi_distal     = psi_distal,
                   se_psi_local1  = se_psi_local1,
                   se_psi_local2  = se_psi_local2,
                   se_psi_local   = se_psi_local,
                   se_psi_distal1 = se_psi_distal1,
                   se_psi_distal0 = se_psi_distal0,
                   se_psi_distal  = se_psi_distal)
  return(est_list)
}



#----- Estimates by transporting distal to local (psi_trans)
trans_est <- function(data){
  dataS0_A0  <- subset(data, S == 0 & A == 0)
  dataS0_A1  <- subset(data, S == 0 & A == 1)
  
  # Estimates: Hájek estimator
  psi_trans0 <- sum(dataS0_A0$w_trans*dataS0_A0$Y)/sum(dataS0_A0$w_trans)
  psi_trans1 <- sum(dataS0_A1$w_trans*dataS0_A1$Y)/sum(dataS0_A1$w_trans)
  
  psi_trans  <- psi_trans1 - psi_trans0
  
  # Sandwich variance
  var_psi_trans0 <- sum((dataS0_A0$w_tran*(dataS0_A0$Y - psi_trans0))^2)/(sum(dataS0_A0$w_tran))^2
  var_psi_trans1 <- sum((dataS0_A1$w_tran*(dataS0_A1$Y - psi_trans1))^2)/(sum(dataS0_A1$w_tran))^2
  
  se_psi_trans   <- sqrt(var_psi_trans0 + var_psi_trans1)
  
  est_list   <- list(psi_trans0    = psi_trans0,
                     psi_trans1    = psi_trans1, 
                     psi_trans     = psi_trans,
                     se_psi_trans1 = sqrt(var_psi_trans1),
                     se_psi_trans0 = sqrt(var_psi_trans0),
                     se_psi_trans  = se_psi_trans)
  return(est_list)
}




#----- Combine function for simulation example
#===============================================================================
sim_combine <- function(data){
  
  dataS1        <- subset(data, S == 1)
  
  # True values
  true_Y2       <- mean(dataS1$Y2)
  true_Y0       <- mean(dataS1$Y0)
  true_fusion   <- true_Y2 - true_Y0
  
  # Estimates
  trial_list    <- trial_est(data)
  trans_list    <- trans_est(data)
  
  # Fusion estimate
  psi_fusion    <- trans_list$psi_trans - trial_list$psi_local
  se_psi_fusion <- sqrt(trans_list$se_psi_trans^2 + trial_list$se_psi_local^2)
  
  # Crude estimate 
  psi_crude     <- trial_list$psi_local2 - trial_list$psi_distal0
  se_psi_crude  <- sqrt(trial_list$se_psi_local2^2 + trial_list$se_psi_distal0^2)
  
  # Difference between distal and trans for active control (A = 1)
  C_diff        <- trans_list$psi_trans1 - trial_list$psi_distal1
  
  # Bias: absolute
  bias_fusion   <- true_fusion - psi_fusion
  bias_crude    <- true_fusion - psi_crude
  
  est_list      <- data.frame(true_fusion   = true_fusion,
                              psi_fusion    = psi_fusion, 
                              psi_crude     = psi_crude,
                              C_diff        = C_diff,
                              bias_fusion   = bias_fusion,
                              bias_crude    = bias_crude,
                              se_psi_fusion = se_psi_fusion,
                              se_psi_crude  = se_psi_crude)
  return(est_list)
}






#----- Bootstrap
#===============================================================================
sim_boot <- function(B, data) {
  # Subset the data by S and A groups
  dataS1A2 <- subset(data, S == 1 & A == 2)
  dataS1A1 <- subset(data, S == 1 & A == 1)
  dataS0A1 <- subset(data, S == 0 & A == 1)
  dataS0A0 <- subset(data, S == 0 & A == 0)
  
  # Get the number of observations
  n12 <- nrow(dataS1A2)
  n11 <- nrow(dataS1A1)
  n01 <- nrow(dataS0A1)
  n00 <- nrow(dataS0A0)
  
  # List to store bootstrap results
  boot_results <- vector("list", B)
  
  for (i in 1:B) {
    # Resample with replacement from each group
    dataS1A2_b <- slice_sample(dataS1A2, n = n12, replace = TRUE)
    dataS1A1_b <- slice_sample(dataS1A1, n = n11, replace = TRUE)
    dataS0A1_b <- slice_sample(dataS0A1, n = n01, replace = TRUE)
    dataS0A0_b <- slice_sample(dataS0A0, n = n00, replace = TRUE)
    
    # Combine the bootstrap samples
    data_boot <- rbind(dataS1A2_b, dataS1A1_b, dataS0A1_b, dataS0A0_b)
    
    # Calculate estimates on the bootstrap sample
    temp_result <- sim_combine(data_boot)
    
    # Store the result
    boot_results[[i]] <- temp_result
  }
  
  # Combine all bootstrap results into one data frame
  boot_df <- do.call(rbind, boot_results)
  return(boot_df)
}







#===============================================================================
#----------------------- Run simulation: Estimates + SE ------------------------ 
#===============================================================================

num_cores <- parallel::detectCores() - 4

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define simulation grid parameters
n1_vals    <- c(200, 500, 1000, 2000)
n2_vals    <- c(200, 500, 1000, 2000)
iterations <- 10000

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
  sim_results <- foreach(i = 1:iterations, .combine = rbind) %dopar% {
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
    se_ave_fusion = mean(se_psi_fusion),
    se_ave_crude  = mean(se_psi_crude),
    rmse_fusion   = sqrt(mean(bias_fusion^2)),
    rmse_crude    = sqrt(mean(bias_crude^2)),
    fusion_mce    = sd(bias_fusion),
    crude_mce     = sd(bias_crude)
  ) |> mutate_all(function(x){x = round(x, 4)})







#===============================================================================
#------------------------------ Inference: 95% CI ------------------------------ 
#===============================================================================

#----- Finding coverage of Asymptotic 95%CI
#===============================================================================
result_tab <- readRDS("Example2_est.RDS")

# Probability scale
param_grid <- result_tab |>
  group_by(n1, n2) |>
  summarise(
    true_fusion   = mean(true_fusion),
    psi_fusion    = mean(psi_fusion)) |>
  ungroup()

# Number of Monte Carlo replicates
MC <- 10000

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

start_time <- Sys.time()

set.seed(12345)

# Parallel loop
coverage_results <- foreach(i = 1:nrow(param_grid), .combine = rbind, 
                            .packages = c("dplyr")) %dopar% {
                              
                              n1 <- param_grid$n1[i]
                              n2 <- param_grid$n2[i]
                              true_fusion <- param_grid$true_fusion[i]
                              
                              cover_fusion  <- 0
                              cover_crude   <- 0
                              
                              # Run MC replicates for this combination
                              for (j in 1:MC) {
                                # Simulate data and generate weights
                                data_sim <- sim_data(n1 = n1, n2 = n2) |> gen_weight()
                                
                                sim_re <- sim_combine(data_sim)
                                
                                # Calculate the asymptotic 95% CI for the fusion estimator
                                lower_ci <- sim_re$psi_fusion - 1.96 * sim_re$se_psi_fusion
                                upper_ci <- sim_re$psi_fusion + 1.96 * sim_re$se_psi_fusion
                                
                                lower_ci_cr <- sim_re$psi_crude - 1.96 * sim_re$se_psi_crude
                                upper_ci_cr <- sim_re$psi_crude + 1.96 * sim_re$se_psi_crude
                                
                                # Check whether the true fusion effect falls inside the CI
                                
                                if (true_fusion >= lower_ci && true_fusion <= upper_ci) {
                                  cover_fusion <- cover_fusion + 1
                                }
                                
                                # Crude 
                                
                                if (true_fusion >= lower_ci_cr && true_fusion <= upper_ci_cr) {
                                  cover_crude <- cover_crude + 1
                                }
                                
                              }
                              
                              # Compute the coverage probability
                              cover_fusion <- cover_fusion / MC
                              
                              cover_crude <- cover_crude / MC
                              
                              data.frame(n1 = n1, n2 = n2, 
                                         cover_fusion = cover_fusion,
                                         cover_crude = cover_crude)
                            }

stopCluster(cl)

# Runtime
end_time <- Sys.time()
runtime <- end_time - start_time
runtime

coverage_results <- coverage_results |> mutate_all(function(x){x = round(x, 4)})







#----- Bootstrap SE: averaged over 10,000 runs
#===============================================================================
result_tab <- readRDS("Example2_est.RDS")

# Probability scale
param_grid <- result_tab |>
  group_by(n1, n2) |>
  summarise(
    true_fusion   = mean(true_fusion),
    psi_fusion    = mean(psi_fusion)) |>
  ungroup()

# Number of Monte Carlo replicates
MC <- 200

cl <- makeCluster(detectCores() - 4)
registerDoParallel(cl)

start_time <- Sys.time()

# set.seed(12345)

results_boot_se <- foreach(i = 1:nrow(param_grid), .combine = rbind,
                           .packages = c("survey", "dplyr")) %dopar% {
                             n1 <- param_grid$n1[i]
                             n2 <- param_grid$n2[i]
                             
                             boot_se_fusion <- numeric(MC)
                             boot_se_crude  <- numeric(MC)
                             
                             for (j in 1:MC) {
                               # Simulate data and generate weights
                               data_sim <- sim_data(n1 = n1, n2 = n2) |> gen_weight()
                               
                               # Run the bootstrap procedure 
                               boot_re <- sim_boot(B = 2000, data = data_sim)
                               
                               boot_se_fusion[j] <- sd(boot_re$psi_fusion)
                               boot_se_crude[j]  <- sd(boot_re$psi_crude)
                             }
                             
                             # Compute the average bootstrap SE for the current (n1, n2) combination
                             boot_se_fusion_ave <- mean(boot_se_fusion)
                             boot_se_crude_ave <- mean(boot_se_crude)
                             
                             # Return a data frame row with the sample sizes and average bootstrap SE
                             data.frame(n1 = n1, n2 = n2, 
                                        boot_se_fusion = boot_se_fusion,
                                        boot_se_crude = boot_se_crude)
                           }


stopCluster(cl)

# Runtime
end_time <- Sys.time()
runtime <- end_time - start_time
runtime

results_boot_se <- results_boot_se |> group_by(n1, n2) |>
  summarise(boot_se_fusion_ave = mean(boot_se_fusion)) |> mutate_all(function(x){x = round(x, 4)})






#----- Compare 95%CI: Supplemental Table 1
#===============================================================================
n1_values <- c(200, 500, 1000, 2000)
n2_values <- c(200, 500, 1000, 2000)
param_grid <- expand.grid(n1 = n1_values, n2 = n2_values)

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

start_time <- Sys.time()

set.seed(12345)
# Parallel loop over each parameter combination
boot_results <- foreach(i = 1:nrow(param_grid), .combine = rbind, 
                        .packages = c("survey", "dplyr")) %dopar% {
                          n1 <- param_grid$n1[i]
                          n2 <- param_grid$n2[i]
                          
                          data_sim <- sim_data(n1 = n1, n2 = n2) |> gen_weight()
                          
                          # Obtain the combined estimates
                          sim_re <- sim_combine(data_sim)
                          
                          # Perform the bootstrap
                          boot_re <- sim_boot(B = 5000, data = data_sim)
                          
                          # Asymptotic 95%CI
                          fusion_asymp_lower <- sim_re$psi_fusion - 1.96 * sim_re$se_psi_fusion
                          fusion_asymp_upper <- sim_re$psi_fusion + 1.96 * sim_re$se_psi_fusion
                          
                          crude_asymp_lower <- sim_re$psi_crude - 1.96 * sim_re$se_psi_crude
                          crude_asymp_upper <- sim_re$psi_crude + 1.96 * sim_re$se_psi_crude
                          
                          # Bootstrap 95%CI
                          fusion_boot_CI <- quantile(boot_re$psi_fusion, probs = c(0.025, 0.975))
                          crude_boot_CI  <- quantile(boot_re$psi_crude, probs = c(0.025, 0.975))
                          
                          re_list <- data.frame(
                            n1 = n1,
                            n2 = n2,
                            psi_fusion         = sim_re$psi_fusion,
                            fusion_asymp_se    = sim_re$se_psi_fusion,
                            fusion_asymp_lower = fusion_asymp_lower,
                            fusion_asymp_upper = fusion_asymp_upper,
                            fusion_boot_se     = sd(boot_re$psi_fusion),
                            fusion_boot_lower  = fusion_boot_CI[1],
                            fusion_boot_upper  = fusion_boot_CI[2],
                            psi_crude          = sim_re$psi_crude,
                            crude_asymp_se     = sim_re$se_psi_crude,
                            crude_asymp_lower  = crude_asymp_lower,
                            crude_asymp_upper  = crude_asymp_upper,
                            crude_boot_se      = sd(boot_re$psi_crude),
                            crude_boot_lower   = crude_boot_CI[1],
                            crude_boot_upper   = crude_boot_CI[2]
                          )
                          
                          return(re_list)
                        }

stopCluster(cl)

# Runtime
end_time <- Sys.time()
runtime <- end_time - start_time
runtime

boot_results <- boot_results |> mutate_all(function(x){x = round(x, 4)})

