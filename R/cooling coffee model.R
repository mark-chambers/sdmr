## ---------------------------
## Script name: 
## Purpose of script: SD model of X ...
## Project name: 
## Project reference number: 
##
## Author: chambm06
## Date created: 2023-11-12
## Date last updated: 
## Copyright (c) Kent Public Health Observatory, Kent County Council, 2023
## Email: kpho@kent.gov.uk
##
## ---------------------------
## Notes:
##   This R script contains a simple example of an sd model slightly modified
##   from the example in the deSolve package. The stocks, parameters and model
##   function will need to be modified to match the project.
##
## ---------------------------

## ---- Settings ---- ##

# Setup simulation times and time step
start <- 0
finish <- 60
time_step <- 1

# Create time vector for simulation time
sim_time <- seq(start, finish, by = time_step)

## ---- Packages ---- ##

library(deSolve)
library(ggplot2)
library(tidyr)


## ---- Get data ---- ##



## ---- Processing ---- ##

# Define stocks for system dynamic model
stocks <- c(coffeeTemp = 80)

# Rates and constants for sd model
pars  <- c(cMass                    = 100,   # mass in grams
           cRoomTemp                = 20,    # temperature of room in degrees
           cSpecificHeat            = 1,     # energy heat 1g water by 1 degree
           cSurfaceAreaCup          = 0.03,  # external surface area of cup
           cHc                      = 50,  # convective heat transfer coeff watts per m2
           cWattPerHrCalories       = 859.845,  # convert watts/hr to calories
           rAbilityCupToReleaseHeat = 0.1)   # zero = no heat loss, 1 = max 

# Create sd model as a function
sd_mod <- function(time, stocks, pars) {
  with(as.list(c(stocks, pars)), {
    
    # Calcs
    calories          <- coffeeTemp * cSpecificHeat * cMass
    tempGap           <- coffeeTemp - cRoomTemp
    watts             <- cSurfaceAreaCup * cHc * tempGap
    caloriesPerMinute <- (watts * cWattPerHrCalories) / 60
    
    # Flows
    heatLossPerMinute <- rAbilityCupToReleaseHeat * caloriesPerMinute / (cMass * cSpecificHeat)

    
    # Changes in stocks (d = delta or difference)
    dCoffeeTemp       <- (- heatLossPerMinute)
    
    
    # Return changes (in same order as stocks vector)
    return(list(c(dCoffeeTemp)))
  })
}
# 
# Run model using ordinary differential equations function from deSolve package
# Note that this creates a matrix.
out <- ode(stocks, sim_time, sd_mod, pars)

## ---- Outputs ---- ##

# See the output
View(out)

# Plot all the stocks over time
plot(out)

# To use in ggplot, must be converted to dataframe.
out_df <- data.frame(out)
out_df <- pivot_longer(out_df, !time, names_to = "state", values_to = "count")
ggplot(data = out_df, mapping = aes(x = time, y = count, colour = state, linetype = state)) + 
  geom_line()
