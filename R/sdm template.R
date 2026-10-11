## ---------------------------
## Script name: 
## Purpose of script: 
## Project name: 
## Project reference number: 

## Author: ADD YOUR NAME
## Date created: ADD DATE
## Date last updated: 
## Copyright (c) Kent Public Health Observatory, Kent County Council, ADD YEAR
## Email: kpho@kent.gov.uk

## ---------------------------
## Notes:
## init gives the starting stock values; state is the current stock-value 
## vector passed into the function; and dA through dD are rates, not stock
## values.
## ode() repeatedly calls sd_mod() with the current time, stock values, and
## parameters, then integrates the returned rates to estimate updated stocks.
## With method = "lsoda", internal step sizes are chosen by the solver and
## need not match times; times specifies when results are returned.
##
## ---------------------------



## ---- Packages ---- ##
library(deSolve)
library(tidyr)
library(ggplot2)

## ---- Settings ---- ##

# Filepaths
# Will be a data sub folder ...

# Colours
col_1 <- "#2D7C82"
col_2 <- "#A9A9A9"
my_linewidth <- 1.5


## ---- Get data ---- ##




## ---- Model setup ---- ##

# Parameters
# All parameters used in the model must be defined here
# Use lower case short names
parameters <- c(
  param1 = 0.9,   # description parameter 1
  param2 = 1/5,  # description parameter 2
  param3 = 1/70,  # description parameter 3
  unused = 100,  # doesn't matter if params defined but not used
  n = 10000  # total population
)

# What unit does the model measure?
value_unit <- "People"

# Initial states of Stocks
# Reserve capital first letters for stock names to make it easier to identify them in the function
# These names must match those used in the sd_mode function below
# And must be in the same order
init <- c(
  A = parameters[["n"]] - 1, 
  B = 0, 
  C = 1, 
  D = 0
)

# Time frame
# Must be numeric and in increasing order
times <- seq(0, 160, by = 1)
time_unit <- "days"

# Validation
# times is finite and increasing
if (!is.numeric(times) || any(!is.finite(times)) || !all(diff(times) > 0)) {
  stop("'times' must be a non-negative numeric vector of finite values in increasing order.")
}
# check that stock and parameter names are unique
if (length(unique(names(init))) != length(init)) {
  stop("Stock names in 'init' must be unique.")
}
if (length(unique(names(parameters))) != length(parameters)) {
  stop("Parameter names in 'parameters' must be unique.")
}
# check that names are disjoint between stocks and parameters
if (length(intersect(names(init), names(parameters))) > 0) {
  stop("Stock names and parameter names must be disjoint (no overlap).")
}

# Define a function that returns the rates of change for each state variable
# The function must follow this format
sd_mod <- function(time, state, parameters) {
  # Make the named stocks and parameters available by name in the equations.
  with(as.list(c(state, parameters)), {
    # Each d-variable is the rate of change of its corresponding stock.
    dA <- -param1 * A * C / n
    dB <- param1 * A * C / n - param2 * B
    dC <- param2 * B - param3 * C
    dD <- param3 * C

    # Return rates in the same order as the initial stocks in init.
    # deSolve integrates these rates to calculate the stock values.
    return(list(c(dA, dB, dC, dD)))
  })
}

## ---- Processing ---- ##

# The ode() function from deSolve numerically integrates the system of equations.
output <- deSolve::ode(
    y = init,          # Starting values for the stocks; updated as the model runs
    times = times,     # Times at which to return results, not necessarily every solver calculation
    func = sd_mod,     # Function that calculates the rates of change
    parms = parameters,# Parameters supplied to sd_mod on each call
    method = "lsoda"  # default method, suitable for most problems
)

output <- as.data.frame(output)

## ---- Outputs ---- ##

# Convert to long format
output_long <- tidyr::pivot_longer(
  output,
  cols = -time,
  names_to = "Stock", 
  values_to = "Value"
)

# Plot using ggplot2
p <- ggplot2::ggplot(output_long, aes(x = time, y = Value, color = Stock)) +
  geom_line(linewidth = my_linewidth) +
  labs(title = "System dynamic model",
       x = paste0("Time (", time_unit, ")"),
       y = value_unit) +
  theme_minimal()

print(p)
