
library(deSolve)

# Define a function that returns the rates of change for each state variable
# The function must follow this format
SEIR <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Rates of change or net flows of people
    dS <- -beta * S * I / N  # Susceptible to Exposed
    dE <- beta * S * I / N - sigma * E # Exposed to Infectious minus Infectious to Recovered
    dI <- sigma * E - gamma * I  # Infectious to Recovered minus recovered individuals accumulating
    dR <- gamma * I  # Recovered individuals accumulating
    return(list(c(dS, dE, dI, dR)))  # must return a list
  })
}

# SEIR Explainer: 
# time: The current time point (not used directly in this function but required by ODE solvers).
# state: A named vector of current values for S, E, I, and R.
# parameters: A named vector of model parameters: beta, sigma, gamma, and N.
# The with(as.list) part combines state and parameters into a single list and unpacks them so you can refer to S, E, I, R, beta, sigma, gamma, and N directly.
# For rates of change, d prefixes the variable names to indicate they are derivatives (rates of change) of the respective compartments.

# Parameters
# All parameters used in the model must be defined here
parameters <- c(
  beta = 0.9,   # transmission rate
  sigma = 1/5,  # incubation rate (1/incubation period)
  gamma = 1/70,   # recovery rate (1/infectious period)
  N = 10000  # total population
)


# Initial states
# These names must match those used in the SEIR function
# And must be in the same order
init <- c(S = unname(parameters["N"]) - 1, E = 0, I = 1, R = 0)


# Time frame
# Must be numeric and in increasing order
times <- seq(0, 160, by = 1)

# The ode() function from deSolve numerically integrates the system of equations.
output <- ode(
    y = init, 
    times = times, 
    func = SEIR, 
    parms = parameters,  # parameters passed to the SEIR function
    method = "lsoda"  # default method, suitable for most problems
)



# Convert output to a data frame for easier handling
output <- as.data.frame(output)


matplot(output$time, output[,2:5], type = "l", lty = 1, col = 1:4,
        xlab = "Time (days)", ylab = "Population",
        main = "SEIR Model Dynamics")
legend("left", legend = c("Susceptible", "Exposed", "Infectious", "Recovered"),
       col = 1:4, lty = 1)




