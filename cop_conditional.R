library(quantmod)
library(copula)
library(urca)


copula_conditional <- function(stocks, train_start, train_end, test_start, test_end, p1, p2, fixed_cost, percentage_cost){
  # training data
  
  adjusted_col1 <- paste0(stocks[1], ".Close")
  adjusted_col2 <- paste0(stocks[2], ".Close")
  
  df_train_1 <- as.data.frame(getSymbols(stocks[1], src='yahoo', auto.assign=FALSE, from=train_start, to=train_end))[adjusted_col1]
  df_train_2 <- as.data.frame(getSymbols(stocks[2], src='yahoo', auto.assign=FALSE, from=train_start, to=train_end))[adjusted_col2]
  
  df_train <- data.frame(s1 = df_train_1[[adjusted_col1]], s2 = df_train_2[[adjusted_col2]])
  
  #Compare copulas (older code)
  
  get_stock_data<-function(stocks, train_start, train_end){
    vec<-list()
    for (i in 1:length(stocks)){
      data <- as.data.frame(getSymbols(stocks[i], src = "yahoo", from = train_start, to = train_end, auto.assign = FALSE))
      adjusted_col <- paste0(stocks[i], ".Adjusted")
      adjusted_prices <- data[[adjusted_col]]
      
      vec[[i]]<-adjusted_prices
    }
    df<-data.frame(vec)
    colnames(df)<-stocks
    return (df)
  }
  
  df<-(get_stock_data(stocks, train_start, train_end))
  
  convert_to_pobs <- function(data) {
    # Convert data to matrix
    data_matrix <- as.matrix(data)
    
    # Convert data to pseudo-observations
    pobs_data <- pobs(data_matrix)
    
    # Check if pseudo-observations contain NA values
    if (any(is.na(pobs_data))) {
      stop("Pseudo-observations contain NA values.")
    }
    
    return(pobs_data)
  }
  
  fit_and_evaluate_copula <- function(data, copula_type) {
    # Convert data to pseudo-observations
    pobs_data <- convert_to_pobs(data)
    
    # Fit the copula model and handle errors
    tryCatch({
      if (copula_type == "Student-t") {
        copula_fit <- fitCopula(tCopula(dim=2), pobs_data, method="ml")
      } else if (copula_type == "Clayton") {
        copula_fit <- fitCopula(claytonCopula(dim=2), pobs_data, method="ml")
      } else if (copula_type == "Gumbel") {
        copula_fit <- fitCopula(gumbelCopula(dim=2), pobs_data, method="ml")
      } else {
        stop(paste("Unknown copula type:", copula_type))
      }
      
      # Calculate Log-Likelihood, AIC, and BIC
      log_likelihood <- logLik(copula_fit)
      n <- nrow(data)
      p <- length(copula_fit@estimate)
      aic <- -2 * log_likelihood + 2 * p
      bic <- -2 * log_likelihood + log(n) * p
      
      return(c(LogLikelihood = log_likelihood, AIC = aic, BIC = bic))
      
    }, error = function(e) {
      warning(paste("Error fitting copula of type", copula_type, ":", e$message))
      return(c(LogLikelihood = NA, AIC = NA, BIC = NA))
    })
  }
  
  compare_copulas <- function(stocks, train_start, train_end) {
    # Get the stock data
    data <- get_stock_data(stocks, train_start, train_end)
    
    # Check data before fitting copulas
    if (ncol(data) < 2) {
      stop("Not enough columns in data for copula fitting.")
    }
    
    # Define copula types
    copulas <- c("Student-t", "Clayton", "Gumbel")
    
    # Fit and evaluate each copula
    results <- sapply(copulas, function(copula) {
      fit_and_evaluate_copula(data, copula)
    }, simplify = "data.frame")
    
    return(results)
  }
  
  # Run comparison
  copulas <- c("Student-t", "Clayton", "Gumbel")
  results <- compare_copulas(stocks, train_start, train_end)
  mini <- 2147483647
  
  if(results[2]<mini){
    mini <- results[2]
    copu <- "Student-t"
  }
  if(results[5]<mini){
    mini <- results[5]
    copu <- "Clayton"
  }
  if(results[8]<mini){
    mini <- results[8]
    copu <- "Gumbel"
  }
  
  #
  
  # testing data
  
  df_test_1 = as.data.frame(getSymbols(stocks[1], src='yahoo', auto.assign=FALSE, from=test_start, to=test_end))[adjusted_col1]
  df_test_2 = as.data.frame(getSymbols(stocks[2], src='yahoo', auto.assign=FALSE, from=test_start, to=test_end))[adjusted_col2]
  
  df_test <- data.frame(s1 = df_test_1[[adjusted_col1]], s2 = df_test_2[[adjusted_col2]])
  
  # Johansen cointegration test on training data
  jtest_train <- ca.jo(df_train, type="trace", K=2, ecdet="none", spec="longrun")
  
  # Extract cointegration vectors and calculate hedge ratio
  cointegration_vectors_train <- jtest_train@V
  hedge_ratio <- cointegration_vectors_train[1, 1] / cointegration_vectors_train[2, 1]
  
  # Calculate empirical CDF values for training data
  u_s1_train <- ecdf(df_train$s1)(df_train$s1)
  u_s2_train <- ecdf(df_train$s2)(df_train$s2)
  
  # Prepare matrix for copula fitting on training data
  mat_train <- cbind(u_s1_train, u_s2_train)
  
  if(copu=="Student-t"){
    fitted <- tCopula(dim=2)
  }else if(copu=="Clayton"){
    fitted <- claytonCopula(dim=2)
  }else{
    fitted <- gumbelCopula(dim=2)
  }
  
  # Fit Clayton copula on training data (best fitted according to AIC)
  cop <- fitCopula(fitted, mat_train, method="itau")
  
  # Calculate empirical CDF values for testing data
  u_s1_test <- ecdf(df_test$s1)(df_test$s1)
  u_s2_test <- ecdf(df_test$s2)(df_test$s2)
  
  # Prepare matrix for copula fitting on testing data
  mat_test <- cbind(u_s1_test, u_s2_test)
  
  # Initialize variables for trading strategy
  initial_capital <- 100  # Starting capital
  money <- initial_capital
  position <- 0
  equity <- c()  # Vector to store returns for each trading day
  
  # Trading logic on testing data
  for (i in 1:nrow(df_test)) {
    p_copula_current <- pCopula(c(mat_test[i, 1], mat_test[i, 2]), cop@copula)
    p_copula_s1 <- pCopula(c(1, mat_test[i, 2]), cop@copula)
    p_copula_s2 <- pCopula(c(mat_test[i, 1], 1), cop@copula)
    trade_cost <- 0
    
    if (((p_copula_current / p_copula_s1) <= p1) | ((p_copula_current / p_copula_s2) >= p2) & (position == 0)) {
      money <- money - (hedge_ratio * df_test$s1[i]) + (df_test$s2[i])
      trade_value <- abs(hedge_ratio * df_test[,1][i]) + abs(df_test[,2][i])
      trade_cost <- fixed_cost + (percentage_cost * trade_value)
      money <- money - trade_cost
      if(money!=0){
        equity <- c(equity, money)
      }
      position <- 1
      
    } else if (((p_copula_current / p_copula_s1) >= p2) & ((p_copula_current / p_copula_s2) <= p1) & (position == 1)) {
      money <- money + (hedge_ratio * df_test$s1[i]) - (df_test$s2[i])
      trade_value <- abs(hedge_ratio * df_test[,1][i]) + abs(df_test[,2][i])
      trade_cost <- fixed_cost + (percentage_cost * trade_value)
      money <- money - trade_cost
      if(money!=0){
        equity <- c(equity, money)
      }
      position <- 0
    }
    
  }# Print final result and Sharpe Ratio
  eq <- c()
  for (i in 1:length(equity)){
    if(equity[i]!=0){
      eq <- c(eq, equity[i])
    }
  }
  
  ret <- c()
  for (i in 2:length(eq)){
    ret <- c(ret, ((eq[i]-eq[i-1])/eq[i-1]))
  }
  
  volatility <- sd(ret)
  
  return (list(money = money, volatility=volatility))
}


stocks <- c('AAPL', 'GOOG')
train_start <- '2015-01-01'
train_end <-'2020-01-01'
test_start <- '2021-01-01'
test_end <- '2023-12-31'
p1 <- 0.05
p2 <- 0.95
fixed_cost <- 10
percentage_cost <- 0.001

res <- copula_conditional(stocks, train_start, train_end, test_start, test_end, p1, p2, fixed_cost, percentage_cost)
print(res$money)
print(res$volatility)