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
## Template for SMALL models: a handful of stocks, no arrays/dimensions.
## (For models split by age, sex, locality etc. use the array template.)
## It runs as it stands (a simple SIR). Work down the sections and replace each
## with your own model. Look for the <<< EDIT markers.
##
## Naming conventions:
## - Stocks: capital letters (S, I, R). Parameters: short lower-case names.
## - Flows: name them (infection, recovery), then build the d-variables from them.
## - d-variables (dS, dI, ...) are rates of change, not stock values.
##
## How ode() works:
## 1. y = init gives ode() the starting values for all stocks.
## 2. ode() calls sd_mod() with the current time, current stock values, and
##    parms = parameters. The argument name 'state' is user-chosen.
## 3. sd_mod() calculates each stock's net flow (inflows minus outflows),
##    in stock units per time.
## 4. sd_mod() returns the rates in a list, in the same order as the stocks in
##    init. These rates are not the updated stock values.
## 5. The solver integrates the rates to estimate updated stock values. With
##    method = "lsoda" internal step sizes are chosen by the solver and need not
##    match `times`; `times` specifies when results are returned.
##
## Scenarios: each scenario is a set of overrides to `parameters` and/or `init`.
## Every scenario is run and plotted together. One scenario ("base") is fine.
## ---------------------------


## ---- 1. Packages ---- ##
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)


## ---- 2. Settings ---- ##

# Filepaths
# dir_data <- "data"                      # <<< EDIT: data sub folder

# Units (used to label the plot)
value_unit <- "People"                    # <<< EDIT: what does the model measure?
time_unit  <- "days"                      # <<< EDIT: every rate below is per this unit

# Can the stocks be negative? TRUE for counts of people or money in a pot.
# FALSE for quantities like temperature, which can legitimately be below zero.
non_negative <- TRUE                      # <<< EDIT

# Colours (one per scenario, cycled) and line width
base_cols    <- c("#2D7C82", "#A9A9A9", "#D98E04", "#8C4A8C", "#B5413B")
my_linewidth <- 1.5
plot_title   <- "System dynamic model"    # <<< EDIT


## ---- 3. Get data ---- ##
# Read real inputs here and use them in section 4.
# dat <- read.csv(file.path(dir_data, "inputs.csv"))                     # <<< EDIT


## ---- 4. Model setup ---- ##

# Parameters (the Stella constants). Short lower-case names, with units.
parameters <- c(                                                 # <<< EDIT
  beta  = 0.4,    # transmission rate, per day
  gamma = 0.2     # recovery rate, per day (1 / infectious period)
)

# Initial states of the stocks. Names must match those used in sd_mod below
# and must be in the same order as the d-variables it returns.
init <- c(                                                       # <<< EDIT
  S = 9990,
  I = 10,
  R = 0
)

# Time frame: numeric, finite, increasing
times <- seq(0, 160, by = 1)                                     # <<< EDIT

# Scenarios: a named list. Each element can override `parameters` and/or `init`
# by name; anything not mentioned keeps its base value above.
scenarios <- list(                                               # <<< EDIT
  "base" = list(),
  "lower transmission" = list(parameters = c(beta = 0.3)),
  "bigger seed"        = list(init = c(S = 9900, I = 100))
)


## ---- 5. Validation ---- ##
if (!is.numeric(times) || any(!is.finite(times)) || !all(diff(times) > 0)) {
  stop("'times' must be a numeric vector of finite values in increasing order.")
}
if (length(unique(names(init))) != length(init) || is.null(names(init))) {
  stop("Stock names in 'init' must be present and unique.")
}
if (length(unique(names(parameters))) != length(parameters) || is.null(names(parameters))) {
  stop("Parameter names in 'parameters' must be present and unique.")
}
if (length(intersect(names(init), names(parameters))) > 0) {
  stop("Stock names and parameter names must be disjoint (no overlap).")
}
if (any(!is.finite(init)) || (non_negative && any(init < 0))) {
  stop("Initial stock values must be finite", if (non_negative) " and not negative.", ".")
}
if (length(scenarios) < 1 || is.null(names(scenarios)) ||
    any(names(scenarios) == "") || anyDuplicated(names(scenarios))) {
  stop("'scenarios' must be a list with a unique name for every element.")
}
# overrides must name things that exist (catches typos)
for (nm in names(scenarios)) {
  bad <- c(setdiff(names(scenarios[[nm]]$parameters), names(parameters)),
           setdiff(names(scenarios[[nm]]$init), names(init)))
  if (length(bad) > 0) {
    stop("Scenario '", nm, "' overrides unknown name(s): ", paste(bad, collapse = ", "))
  }
}


## ---- 6. The model: auxiliaries, flows, then derivatives ---- ##
# deSolve calls this at every time step. Keep to the same order each time:
#   a) auxiliaries   b) flows   c) derivatives
# A stock's derivative = inflows - outflows. Return them in the same order as `init`.
sd_mod <- function(time, state, parameters) {
  # Make the named stocks and parameters available by name in the equations.
  with(as.list(c(state, parameters)), {

    # a) auxiliaries (variables worked out from the stocks)           <<< EDIT
    N <- S + I + R

    # b) flows (rates, per time unit)
    infection <- beta * S * I / N
    recovery  <- gamma * I

    # c) derivatives: net flow into each stock. Positive values increase the
    # stock, negative values decrease it.
    dS <- -infection
    dI <-  infection - recovery
    dR <-  recovery

    # The first element is the derivatives. Anything after is an extra output
    # column (handy for plotting an auxiliary), e.g. list(c(dS, dI, dR), N = N)
    return(list(c(dS, dI, dR)))
  })
}


## ---- 7. Run every scenario ---- ##
override <- function(x, new) {
  if (length(new) > 0) x[names(new)] <- new
  x
}

run_scenario <- function(sc, method = "lsoda") {
  deSolve::ode(
    y      = override(init, sc$init),               # starting values for the stocks
    times  = times,                                 # times at which to return results
    func   = sd_mod,                                # calculates the rates of change
    parms  = override(parameters, sc$parameters),   # supplied to sd_mod on each call
    method = method                                 # "lsoda": adaptive, handles stiff systems
  )
}
# method = "euler" uses fixed steps at the spacing of `times` (to match Stella's
# Euler, use the same DT as the Stella model).

output <- lapply(names(scenarios), function(nm) {
  as.data.frame(run_scenario(scenarios[[nm]])) |> mutate(scenario = nm)
}) |>
  bind_rows() |>
  mutate(scenario = factor(scenario, levels = names(scenarios)))


## ---- 8. Checks (do these every time) ---- ##
stock_names <- names(init)

if (non_negative) {
  cat("smallest stock value:", signif(min(output[stock_names]), 3), "\n")
}

# Optional: if the model is closed (nothing enters or leaves), the total across
# stocks should not change. Delete this if your stocks are not additive
# (temperature, price, etc.), because the "drift" would be meaningless.
check_conserved <- TRUE                                          # <<< EDIT
if (check_conserved) {
  total <- rowSums(output[stock_names])
  drift <- tapply(total, output$scenario, function(x) max(abs(x - x[1])))
  cat("max drift in total, by scenario:\n"); print(signif(drift, 3))
}

# Better still: if your model has an exact solution, compare against it.


## ---- 9. Outputs ---- ##

# Long format: one row per time x scenario x variable (stocks and any extra outputs)
output_long <- output |>
  pivot_longer(-c(time, scenario), names_to = "Variable", values_to = "Value") |>
  mutate(Variable = factor(Variable, levels = unique(Variable)))

# One panel per variable, one line per scenario
n_scen       <- nlevels(output_long$scenario)
scen_palette <- setNames(rep_len(base_cols, n_scen), levels(output_long$scenario))

p <- ggplot(output_long, aes(x = time, y = Value, colour = scenario)) +
  geom_line(linewidth = my_linewidth) +
  facet_wrap(~ Variable, scales = "free_y") +
  scale_colour_manual(values = scen_palette) +
  labs(title = plot_title,
       x = paste0("Time (", time_unit, ")"),
       y = value_unit,
       colour = NULL) +
  theme_minimal()

if (n_scen == 1) p <- p + theme(legend.position = "none")

print(p)
