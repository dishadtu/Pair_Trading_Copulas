library(quantmod)
library(copula)
library(urca)
library(KFAS)
library(pracma)
library(ggplot2)

# Transaction Fee Parameters
fixed_fee <- 1  # Fixed fee per transaction
percentage_fee <- 0.001  # Percentage of the value of the trade

# Risk-Free Rate Parameters (e.g., 3% annual rate)
risk_free_annual <- 0.03

extract_df <- function(stocks, start, end, type) {
  print(paste("Fetching data for stocks: ", stocks, "from", start, "to", end))
  if (type == "raw") {
    data1 <- getSymbols(stocks[1], src = 'yahoo', auto.assign = FALSE, from = start, to = end)
    data2 <- getSymbols(stocks[2], src = 'yahoo', auto.assign = FALSE, from = start, to = end)
  } else {
    data1 <- getSymbols(stocks[1], src = 'yahoo', auto.assign = FALSE, from = start, to = end)
    data2 <- getSymbols(stocks[2], src = 'yahoo', auto.assign = FALSE, from = start, to = end)
    
    data1 <- Ad(data1)
    data2 <- Ad(data2)
    
    data1 <- diff(log(data1))
    data2 <- diff(log(data2))
  }
  print("Data fetched successfully")
  return(data.frame(s1=data1, s2=data2))
}

# Function to compute transaction fee for a trade
calculate_transaction_fee <- function(trade_value) {
  fee <- fixed_fee + (trade_value * percentage_fee)
  return(fee)
}

# Function to fit copulas and return the best-fit copula
compare_copulas <- function(df) {
  data_matrix <- as.matrix(df)
  u <- pobs(data_matrix)
  
  best_copula <- NULL
  best_loglik <- -Inf
  best_copula_name <- NULL
  
  copula_models <- list(
    Gaussian = normalCopula(),
    Gumbel = gumbelCopula(),
    Clayton = claytonCopula()
  )
  
  for (copula_name in names(copula_models)) {
    copula_model <- copula_models[[copula_name]]
    
    tryCatch({
      fitted_copula <- fitCopula(copula_model, u, method = "ml")
      loglik <- logLik(fitted_copula)
      
      if (loglik > best_loglik) {
        best_loglik <- loglik
        best_copula <- fitted_copula
        best_copula_name <- copula_name
      }
      
    }, error = function(e) {
      warning(paste("Error fitting", copula_name, "copula:", e$message))
    })
  }
  
  if (is.null(best_copula)) {
    stop("Failed to fit any copula to the data.")
  }
  
  cat("Best fitted copula based on log-likelihood:", best_copula_name, "\n")
  return(best_copula)
}

# Rolling function with transaction fee, volatility, and Sharpe ratio including risk-free rate
rolling <- function(stocks, total_df, test_df, p1, p2, window, multiple, last_value) {
  
  cointegration_test <- function(df) {
    jtest <- ca.jo(df, type = "trace", ecdet = "const", K = 2)
    summary(jtest)
    hedge_ratio <- -jtest@V[, 1][2] / jtest@V[, 1][1]
    return(hedge_ratio)
  }
  
  money <- 10000
  position <- 0
  stocks_of_s1 <- 0
  stocks_of_s2 <- 0
  equity <- numeric(nrow(test_df) + 1)
  equity[1] <- money
  transaction_fees <- c()
  return_matrix <- matrix(, nrow=nrow(test_df), ncol=5)
  conditional_df <- data.frame(c1=numeric(), c2=numeric())
  
  u <- pobs(tail(total_df, n=window))
  count <- 0
  time_fitted <- 1
  
  cop <- compare_copulas(tail(total_df, n=window))@copula
  conditional_copula <- cCopula(u, cop)
  new_matrix <- matrix(, nrow=nrow(conditional_copula), ncol=2)
  new_matrix[, 1] <- conditional_copula[, 1]
  new_matrix[, 2] <- conditional_copula[, 2]
  
  df_c <- as.data.frame(new_matrix)
  
  conditional_df <- rbind(conditional_df, df_c)
  
  count <- 0
  time_fitted <- 1
  
  for (i in 1:nrow(test_df)) {
    # Update portfolio value
    current_portfolio_value <- stocks_of_s1 * test_df[i, 1] + stocks_of_s2 * test_df[i, 2]
    
    # Re-fit copula every 'window' periods
    if (count >= window) {
      cop <- compare_copulas(tail(total_df, n=window))@copula
      conditional_copula <- cCopula(u, cop)
      new_matrix <- matrix(, nrow=nrow(conditional_copula), ncol=2)
      new_matrix[, 1] <- conditional_copula[, 1]
      new_matrix[, 2] <- conditional_copula[, 2]
      df_c <- as.data.frame(new_matrix)
      conditional_df <- rbind(conditional_df, df_c)
      count <- 0
      time_fitted <- time_fitted + 1
    } else {
      count <- count + 1
    }
    
    # Update returns in total_df
    if (i == 1) {
      current_return1 <- (test_df[i, 1] / last_value[1, 1]) - 1
      current_return2 <- (test_df[i, 2] / last_value[1, 2]) - 1
    } else {
      current_return1 <- (test_df[i, 1] / test_df[i-1, 1]) - 1
      current_return2 <- (test_df[i, 2] / test_df[i-1, 2]) - 1
    }
    new_row <- data.frame(s1=current_return1, s2=current_return2)
    colnames(new_row) <- colnames(total_df)
    total_df <- rbind(total_df, new_row)
    
    hedge_ratio <- cointegration_test(tail(total_df, n=window))
    mat <- pobs(tail(total_df, n=window))
    
    condition1 <- conditional_df[i, 1]
    condition2 <- conditional_df[i, 2]
    
    return_matrix[i, 1] <- condition1
    return_matrix[i, 2] <- condition2
    return_matrix[i, 3] <- hedge_ratio
    
    # Trading logic
    # Closing positions
    if ((condition1 > 0.5 && condition2 > 0.5) && position != 0) {
      trade_value <- abs(stocks_of_s1 * test_df[i, 1]) + abs(stocks_of_s2 * test_df[i, 2])
      fee <- calculate_transaction_fee(trade_value)
      transaction_fees <- c(transaction_fees, fee)
      
      money <- money + (stocks_of_s1 * test_df[i, 1]) + (stocks_of_s2 * test_df[i, 2]) - fee
      position <- 0
      stocks_of_s1 <- 0
      stocks_of_s2 <- 0
    } 
    # Opening long position in s1, short in s2
    else if ((condition1 <= p1 && condition2 >= p2) && position == 0) {
      trade_value <- (test_df[i, 1] * multiple) + (hedge_ratio * multiple * test_df[i, 2])
      fee <- calculate_transaction_fee(trade_value)
      transaction_fees <- c(transaction_fees, fee)
      
      stocks_of_s1 <- stocks_of_s1 + (1 * multiple)
      stocks_of_s2 <- stocks_of_s2 - (hedge_ratio * multiple)
      money <- money - (test_df[i, 1] * multiple) + (hedge_ratio * multiple * test_df[i, 2]) - fee
      position <- 1
    }
    # Opening short position in s1, long in s2
    else if ((condition1 >= p2 && condition2 <= p1) && position == 0) {
      trade_value <- (test_df[i, 1] * multiple) + (hedge_ratio * multiple * test_df[i, 2])
      fee <- calculate_transaction_fee(trade_value)
      transaction_fees <- c(transaction_fees, fee)
      
      stocks_of_s1 <- stocks_of_s1 - (1 * multiple)
      stocks_of_s2 <- stocks_of_s2 + (hedge_ratio * multiple)
      money <- money + (test_df[i, 1] * multiple) - (hedge_ratio * multiple * test_df[i, 2]) - fee
      position <- -1
    }
    
    # Update equity after trading
    equity[i+1] <- money + (stocks_of_s1 * test_df[i,1]) + (stocks_of_s2 * test_df[i,2])
    
    return_matrix[i, 4] <- stocks_of_s1
    return_matrix[i, 5] <- stocks_of_s2
  }
  
  # Calculate performance metrics
  log_returns <- diff(log(equity))
  volatility_daily <- sd(log_returns, na.rm=TRUE)
  volatility_annual <- volatility_daily * sqrt(252)
  
  initial_value <- equity[1]
  final_value <- equity[length(equity)]
  years <- length(equity) / 252  # Approximate number of years
  CAGR <- ((final_value / initial_value)^(1 / years)) - 1
  
  daily_risk_free <- (1 + risk_free_annual)^(1/252) - 1
  daily_returns <- diff(log(equity))
  excess_returns <- daily_returns - daily_risk_free
  sharpe_ratio <- mean(excess_returns, na.rm=TRUE) / sd(daily_returns, na.rm=TRUE) * sqrt(252)
  
  print("Equity values:")
  print(equity)
  
  print("Daily returns:")
  print(daily_returns)
  
  return(list(
    returns = money,
    transaction_fees = transaction_fees,
    volatility_daily = volatility_daily,
    volatility_annualized = volatility_annual,
    cagr = CAGR,
    equity = equity,
    times = time_fitted,
    sharpe = sharpe_ratio,
    daily = daily_returns,
    mat = return_matrix
  ))
}



# Function to find optimal multiple for the highest Sharpe ratio
find_sharpe_plot <- function(pair, total_df, test_df, p1, p2, window, last_value, multiple_range = seq(1, 100, by = 0.5)) {
  # Initialize a vector to store Sharpe ratios for each multiple value
  sharpe_ratios <- numeric(length(multiple_range))
  
  # Loop through all values of multiple
  for (i in seq_along(multiple_range)) {
    multiple <- multiple_range[i]
    
    # Run the rolling function with the current multiple value
    results <- rolling(pair, total_df, test_df, p1, p2, window, multiple, last_value)
    
    # Store the Sharpe ratio for the current multiple
    sharpe_ratios[i] <- results$sharpe
  }
  
  # Create a data frame for plotting
  plot_data <- data.frame(multiple = multiple_range, sharpe_ratio = sharpe_ratios)
  
  # Plot the Sharpe ratio as a function of the multiple
  ggplot(plot_data, aes(x = multiple, y = sharpe_ratio)) +
    geom_line() +
    geom_point(color = 'blue') +
    labs(title = 'Sharpe Ratio vs. Multiple', x = 'Multiple', y = 'Sharpe Ratio') +
    theme_minimal()
  
  # Find the multiple that maximizes the Sharpe ratio
  optimal_multiple <- multiple_range[which.max(sharpe_ratios)]
  
  # Print and return the results
  cat("Optimal multiple:", optimal_multiple, "\n")
  cat("Maximum Sharpe ratio:", max(sharpe_ratios), "\n")
  
  return(optimal_multiple)
}

# Main execution
#selected_pairs <- list(c("AMZN", "ADBE"))
test_start <- '2021-01-01'
test_end <- '2023-01-01'
buffer_start <- '2019-01-01'
buffer_end <- '2020-12-31'
p1 <- 0.33
p2 <- 0.875
window <- 252
multiple <- 100

# total_df <- extract_df(pair, buffer_start, buffer_end, "returns")
# test_df <- extract_df(pair, test_start, test_end, "raw")
# last_value <- extract_df(pair, buffer_start, buffer_end, "raw")
# last_value <- tail(last_value, 1)
# 
# results <- rolling(pair, total_df, test_df, p1, p2, window, multiple, last_value)
# print(results$sharpe)

# Uncomment to run multiple pairs
final_df <- data.frame()
# selected_pairs<-list(
#   c("AMZN", "ADBE"), c("V", "ACN"), c("V", "ADBE"), c("MA", "LIN"),
#   c("V", "PYPL"), c("HON", "LIN"), c("MA", "DHR"), c("NEE", "LIN"),
#   c("AMZN", "NFLX"), c("HD", "TXN"), c("MA", "COST"), c("WMT", "LIN"),
#   c("V", "COST"), c("AMZN", "PYPL"), c("WMT", "ACN"), c("GOOGL", "ADBE"),
#   c("HD", "LIN"), c("HD", "MA"), c("HD", "ADBE"), c("NFLX", "ADBE"),
#   c("HON", "TXN"), c("MCD", "LIN"), c("UNH", "ADBE"), c("MCD", "NEE"),
#   c("AMZN", "GOOGL"), c("NKE", "COST"), c("ACN", "COST"), c("HON", "DHR"),
#   c("AAPL", "MSFT"), c("MCD", "ACN"), c("UNH", "NFLX"), c("WMT", "COST"),
#   c("CSCO", "ACN"), c("HD", "COST"), c("CSCO", "LIN"), c("MCD", "DHR"),
#   c("CSCO", "HON"), c("NKE", "NEE"), c("MCD", "UNP"), c("CSCO", "DHR"),
#   c("WMT", "TXN"), c("MRK", "LLY"), c("HON", "COST"), c("NKE", "ACN"),
#   c("GOOGL", "NFLX"), c("NKE", "LIN")
# )

selected_pairs<-list(
  c("NKE", "COST"),c("MA","COST") ,c("AMZN", "ADBE"), c("UNH", "ADBE"), c("HON", "DHR"),
  c("UNH", "NFLX"), c("MCD", "ACN"), c("AMZN", "NFLX"), c("NEE", "LIN"),
  c("NFLX", "ADBE"), c("NKE", "NEE"), c("NKE", "ACN"), c("GOOGL", "ADBE"),
  c("CSCO", "HON"),c("AMZN","PYPL"),c("CSCO","ACN"), c("WMT", "COST"), c("CSCO", "LIN"), c("WMT", "LIN"),
  c("NKE", "LIN")
)
# selected_pairs<-list(
#      c("NKE", "COST")
# )
for(pair in selected_pairs){
  total_df <- extract_df(pair, buffer_start, buffer_end, "returns")
  test_df <- extract_df(pair, test_start, test_end, "raw")
  last_value <- extract_df(pair, buffer_start, buffer_end, "raw")
  last_value <- tail(last_value, 1)
  
  # Find optimal multiple and plot the Sharpe ratio for the current pair
  #optimal_multiple <- find_sharpe_plot(pair, total_df, test_df, p1, p2, window, last_value)
  #optimal_multiple <- 1640

  # Run the rolling function and get results
  results <- rolling(pair, total_df, test_df, p1, p2, window, multiple, last_value)

  temp_df <- data.frame(
    stock1=pair[1],
    stock2=pair[2],
    daily_volatility = results$volatility_daily,
    annualized_volatility = results$volatility_annualized,
    sharpe_ratio=results$sharpe,
    cagr=results$cagr
  )

  final_df <- rbind(final_df, temp_df)
}

print(final_df)
write.csv(final_df, file = "Strategy Result_For_Selected_Singular.csv", row.names = FALSE)