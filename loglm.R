# Packages and other settings
library(tidyverse)
library(magrittr)
library(truncnorm)

# Functions: model evaluation metrics
`%mfe%` <- function(a, b) {
  if(length(a) != length(b)) stop("Two vectors have different lengths!")
  if(length(a) == 0)         stop("Vector is length-0!")
  if(any((a + b) == 0))      stop("Value zero is used as a denominator!")
  mean(abs(a - b)/((a + b)/2))
}

`%mfb%` <- function(a, b) {
  if(length(a) != length(b)) stop("Two vectors have different lengths!")
  if(length(a) == 0)         stop("Vector is length-0!")
  if(any((a + b) == 0))      stop("Value zero is used as a denominator!")
  mean((a - b)/((a + b)/2))
}

# Function: original log-transformed linear model
loglm_original <- function(formula, .data, 
                           LOWER = rep(0.1, str_count(Reduce(paste, deparse(formula)), "\\+") + 1),
                           UPPER = rep(10, str_count(Reduce(paste, deparse(formula)), "\\+") + 1)){
  
  lm(formula, .data) -> ols 
  ols %>% coef() %>% pmax(0.00001) -> init
  ols$model[-1] %>% as.matrix() -> X
  
  if(length(init) != ncol(X)) stop("Did you forget -1 in your regression formula?")
  if(length(LOWER) != length(UPPER)) stop("Lower boundary and upper bounday do NOT match!")
  if(length(LOWER) != ncol(ols$model[-1])) stop("Constrain boundary and number of variables do NOT match!")
  
  func0 <- function(para){
    log(as.vector(X %*% para)) - log(ols$model[[1]]) -> difference
    sum(difference^2)
  }
  
  optim(init, fn = func0, method = "L-BFGS-B", lower = LOWER, upper = UPPER)$par -> optimized
  
  optimized
}

# Function: original log-transformed linear model for Monte Carlo Simulation
loglm_original_trial <- function(formula, .data,
                                 std = rep(0.4^2, str_count(Reduce(paste, deparse(formula)), "\\+") + 2),
                                 LOWER = rep(0.1, str_count(Reduce(paste, deparse(formula)), "\\+") + 1),
                                 UPPER = rep(10, str_count(Reduce(paste, deparse(formula)), "\\+") + 1)){
  
  lm(formula, .data) -> ols
  ols$model -> .data2
  
  if(ncol(.data2) != length(std)) stop("Number of uncertainties provided do not match the column number of data!")
  
  modify2(.data2, std, ~.x * rtruncnorm(length(.x), 0, Inf, 1, .y)) -> .data2
  loglm(formula, .data2, LOWER, UPPER)
}


# Function: regularized log-transformed linear model
loglm_reg <- function(formula, .data, lnL = 0,
                      LOWER = rep(0.1, str_count(Reduce(paste, deparse(formula)), "\\+") + 1),
                      UPPER = rep(10, str_count(Reduce(paste, deparse(formula)), "\\+") + 1)){
  
  lm(formula, .data) -> ols 
  ols %>% coef() %>% pmax(0.00001) -> init
  ols$model[-1] %>% as.matrix() -> X
  
  if(length(init) != ncol(X)) stop("Did you forget -1 in your regression formula?")
  if(length(LOWER) != length(UPPER)) stop("Lower boundary and upper bounday do NOT match!")
  if(length(LOWER) != ncol(ols$model[-1])) stop("Constrain boundary and number of variables do NOT match!")
  
  func0 <- function(para){
    log(as.vector(X %*% para)) - log(ols$model[[1]]) -> difference
    sum(difference^2) + exp(lnL) * sum((log(para))^2)
  }
  
  optim(init, fn = func0, method = "L-BFGS-B", lower = LOWER, upper = UPPER)$par -> optimized
  
  sum((log(as.vector(X %*% optimized)) - log(ols$model[[1]]))^2) -> Q0
  sum((log(optimized))^2) -> Qp
  c(Q0 = Q0, Qp = Qp, optimized)
}

# Function: regularized log-transformed linear model for Monte Carlo Simulation
loglm_reg_trial <- function(formula, .data, lnL,
                            std = rep(0.4^2, str_count(Reduce(paste, deparse(formula)), "\\+") + 2),
                            LOWER = rep(0.1, str_count(Reduce(paste, deparse(formula)), "\\+") + 1),
                            UPPER = rep(10, str_count(Reduce(paste, deparse(formula)), "\\+") + 1)){
  
  lm(formula, .data) -> ols
  ols$model -> .data2
  
  if(ncol(.data2) != length(std)) stop("Number of uncertainties provided do not match the column number of data!")
  
  modify2(.data2, std, ~.x * rtruncnorm(length(.x), 0, Inf, 1, .y)) -> .data2
  loglm_reg(formula, .data2, lnL, LOWER, UPPER)
}
