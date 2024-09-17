# Pair_Trading_Copulas
This repository features a copula-based trading strategy in R that models asset pair dependencies using conditional probability. It selects the best copula model (Student-t, Clayton, or Gumbel) based on AIC, BIC, and Log-Likelihood, and uses the Johansen cointegration test for hedge ratio calculation.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Parameters](#parameters)

## Overview

The script performs the following steps:

1. *Data Collection*: Retrieves historical adjusted closing prices for two specified stocks over training and testing periods.
2. *Copula Selection*: Compares Student-t, Clayton, and Gumbel copulas based on AIC to model the dependency structure between the stocks.
3. *Cointegration Analysis*: Performs the Johansen cointegration test on the training data to calculate the hedge ratio.
4. *Trading Strategy*: Implements a conditional trading strategy based on copula probabilities and hedge ratios.
5. *Performance Evaluation*: Calculates the final capital (money) and the volatility of returns.

## Requirements

- R (version >= 3.6.0)
- Packages:
  - quantmod
  - copula
  - urca

Install the required packages using:

R
install.packages(c("quantmod", "copula", "urca"))


## Usage

The main function is copula_conditional, which takes the following parameters:

### Parameters

- stocks: A vector of two stock ticker symbols (e.g., c('AAPL', 'GOOG')).
- train_start, train_end: The start and end dates for the training period (format: 'YYYY-MM-DD').
- test_start, test_end: The start and end dates for the testing period (format: 'YYYY-MM-DD').
- p1, p2: Threshold probabilities for the trading strategy (values between 0 and 1).
- fixed_cost: Fixed transaction cost per trade.
- percentage_cost: Variable transaction cost as a percentage of the trade value.

