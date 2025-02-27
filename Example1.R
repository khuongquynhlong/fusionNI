library(tidyverse)
library(tableone)
library(gmodels)

#===============================================================================
#------------ Example 1: Simple case with single binary covariate X ------------ 
#===============================================================================

#----- Generate data
#===============================================================================
sim_data1 <- function(n1, n2, seed = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  N <- n1 + n2
  X1 <- rbinom(N, 1, 0.4)
  
  #--- Generate S
  # P(S = 1 | X) = plogis(a_est - X1 - X2) such that
  # the expected number of S == 1 equals n1.
  f_a <- function(a) {
    sum(plogis(a + X1)) - n1
  }
  
  # Find the root
  a_est <- uniroot(f_a, interval = c(-10, 10), tol = 1e-10, maxiter = 10000)$root
  S <- rbinom(N, size = 1, prob = plogis(a_est + X1))
  
  #--- Generate A: P(A = 1| S = 1)
  # For individuals with S == 1, A is drawn to be either 1 or 2.
  # For individuals with S == 0, A is drawn to be either 0 or 1.
  A <- rep(NA, N)
  A[S == 1] <- rbinom(sum(S == 1), size = 1, prob = 0.5) + 1  # results in {1, 2}
  A[S == 0] <- rbinom(sum(S == 0), size = 1, prob = 0.5)      # results in {0, 1}
  
  #---- Combine dataset
  gen_data <- data.frame(X1 = X1, S = S, A = A)
  
  #--- Generate potential outcome: Y0, Y1, Y2
  # Y ~ expit(X1 + X2 + 3*I(A = 0) + I(A = 1) + 2*I(A = 2))
  lin_pred0 <- -3 + gen_data$X1 + 2    # for A = 0
  lin_pred1 <- -3 + gen_data$X1 + 1    # for A = 1
  lin_pred2 <- -3 + gen_data$X1 + 1.2  # for A = 2
  
  gen_data$Y0 <- rbinom(n = N, size = 1, prob = plogis(lin_pred0))
  gen_data$Y1 <- rbinom(n = N, size = 1, prob = plogis(lin_pred1))
  gen_data$Y2 <- rbinom(n = N, size = 1, prob = plogis(lin_pred2))
  
  # Observed Y
  # Y = I{A = 0} * Y0 + I{A = 1} * Y1 + I{A = 2} * Y2
  gen_data$Y <- ifelse(gen_data$A == 0, gen_data$Y0,
                       ifelse(gen_data$A == 1, gen_data$Y1, gen_data$Y2))
  
  return(gen_data)
}

n1 = 1000; n2 = 1000
N = n1 + n2

gen_data1 <- sim_data1(n1 = n1, n2 = n2, seed = 12345)

gen_data1 <- gen_data1 |>
  mutate(X1 = as.factor(X1),
         Y = as.factor(Y))

# Create the table by stratifying on A
CreateTableOne(vars = c("X1", "Y"), strata = "A", data = subset(gen_data1, S == 1))
CreateTableOne(vars = c("X1", "Y"), strata = "A", data = subset(gen_data1, S == 0))


# Proportion of Y
tab1 <- gen_data1 |> group_by(S, X1) |>
  summarise(n = round(n()/100))

# Proportion of Y
tab1_1 <- gen_data1 |> group_by(S, A) |>
  summarise(probY = mean(Y == "1"))

# Proportion of Y, by X
tab2 <- gen_data1 |> group_by(S, A, X1) |>
  summarise(probY = mean(Y == "1"),
            n = n())


# Sample
n11 <- tab1$n[4]
n10 <- tab1$n[3]
n01 <- tab1$n[2]
n00 <- tab1$n[1]

# P(S = l | X = 1)
p11 <- n11/ (n11 + n01)
p10 <- n10/ (n10 + n00)

w_trans1 <- p11/(1-p11)*(n2/N)
w_trans0 <- p10/(1-p10)*(n2/N)

psi_trans1 = (tab2$probY[4]*tab2$n[4]*w_trans1 + tab2$probY[3]*tab2$n[3]*w_trans0)/
  (tab2$n[4]*w_trans1 + tab2$n[3]*w_trans0)


psi_trans1 - tab1_1$probY[3]

