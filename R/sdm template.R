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
##
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
init <- c(A = parameters[["n"]] - 1, B = 0, C = 1, D = 0)

# Time frame
# Must be numeric and in increasing order
times <- seq(0, 160, by = 1)
time_unit <- "days"

# Define a function that returns the rates of change for each state variable
# The function must follow this format
sd_mod <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Rates of change or net flows of people
    # For rates of change, d prefixes the variable names to indicate they are derivatives (rates of change) of the respective compartments.
    # These are not flows but rather net changes in stocks
    # at each timestep
    # Replace capital letters with your stock names
    # define how parameters and stocks interact to determine rate of change for stocks
    dA <- -param1 * A * C / n  # Note on change in stock A
    dB <- param1 * A * C / n - param2 * B # Note on change in stock B
    dC <- param2 * B - param3 * C # Note on change in stock C
    dD <- param3 * C # Note on change in stock D
    return(list(c(dA, dB, dC, dD)))  # must return a list
  })
}

## ---- Processing ---- ##

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

# The ode() function from deSolve numerically integrates the system of equations.
output <- deSolve::ode(
    y = init, 
    times = times, 
    func = sd_mod, 
    parms = parameters,  # parameters passed to the sd_mod function
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
