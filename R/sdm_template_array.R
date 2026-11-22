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
## Template for models whose stocks are split by one or more dimensions
## (age, sex, locality, deprivation quintile ...). For a handful of stocks with no
## dimensions, use sdm_template_simple.R instead.
## It runs as it stands (a simple SIR split by age group). Work down the sections
## and replace each with your own model. Look for the <<< EDIT markers.
##
## Naming conventions:
## - Stocks: capital letters (S, I, R). They are stored as arrays with one cell
##   per combination of the dimensions defined in the Dimensions section.
## - Parameters: short lower-case names, all held in the list `parameters`.
##   Each is a single value, or an array shaped like the dimensions.
## - d-variables (dS, dI, ...) are rates of change, not stock values.
##
## How ode() works:
## 1. y = y0 gives ode() the starting values for all stocks (as one flat vector).
## 2. ode() calls sd_mod() with the current time, current stock values, and
##    parms = parameters.
## 3. sd_mod() unpacks the stocks into arrays, works out the auxiliaries and
##    flows, and calculates each stock's net flow (inflows minus outflows),
##    in stock units per time.
## 4. sd_mod() returns the rates in a list, in the same order as the stocks in
##    `init`. These rates are not the updated stock values.
## 5. The solver integrates the rates to estimate updated stock values. With
##    method = "lsoda" internal step sizes are chosen by the solver and need not
##    match `times`; `times` specifies when results are returned.
##
## Final steps: converting to a data frame and reshaping for plotting changes the
## output format, not the model calculation.
## ---------------------------


## ---- 1. Packages ---- ##
library(deSolve)
library(dplyr)
library(tidyr)      # needs version 1.3.0 or later for separate_wider_delim()
library(ggplot2)


## ---- 2. Settings ---- ##

# Filepaths
dir_data <- "data"                       # <<< EDIT: data sub folder

# Units (used to label the plots)
value_unit <- "People"                   # <<< EDIT: what does the model measure?
time_unit  <- "days"                     # <<< EDIT: every rate below is per this unit

# Can the stocks be negative? TRUE for counts of people.
# FALSE for quantities like temperature, which can legitimately be below zero.
non_negative <- TRUE                     # <<< EDIT

# Is the model closed (nothing enters or leaves)? If so the total across all
# stocks should not change. Set FALSE if you have births, deaths, migration in or
# out, events that change stocks, or stocks that cannot be added together.
check_conserved <- TRUE                  # <<< EDIT

# Colours (cycled across stocks, so give at least as many as you have stocks)
base_cols    <- c("#2D7C82", "#A9A9A9", "#D98E04", "#8C4A8C", "#B5413B")
my_linewidth <- 1.5
plot_title   <- "System dynamic model"   # <<< EDIT


## ---- 3. Get data ---- ##
# Read real inputs here (population, rates, etc.) and use them in the Model setup.
# pop <- read.csv(file.path(dir_data, "population.csv"))                 # <<< EDIT


## ---- 4. Dimensions (the Stella "dimensions") ---- ##
# One entry per dimension, in a fixed order. Add more for more dimensions, e.g.
#   dn <- list(age = ..., sex = c("F","M"), loc = c("A","B"), dep = paste0("Q",1:5))
# Level names must not contain the character "|" (it is used to label cells).
dn <- list(                                                      # <<< EDIT
  age = c("0-14", "15-64", "65+")
)

dims    <- unname(lengths(dn))     # size of each dimension
n_cells <- prod(dims)              # number of cells in each stock

# "15-64|F|A|Q3"-style label for every cell, in the same order as as.vector()
cell_labels <- function() {
  do.call(paste, c(expand.grid(dn, stringsAsFactors = FALSE), sep = "|"))
}

# TRUE if x is an array shaped like the dimensions
is_cell_array <- function(x) {
  !is.null(dim(x)) && length(dim(x)) == length(dims) && all(dim(x) == dims)
}

# Spread a vector along dimension k to fill a whole array.
# Use it to build arrays of parameters or populations that vary by one dimension.
#   along(c(1, 0.8, 1.2), 1)   # a parameter that varies by the first dimension
along <- function(v, k) {
  array(v[slice.index(array(0, dims), k)], dim = dims, dimnames = dn)
}

# Division that returns 0 where the denominator is 0 (e.g. a cell with nobody in it).
# Without it, a per-cell I / N gives NaN, which then spreads through the model.
safe_div <- function(a, b) ifelse(b > 0, a / b, 0)

# Full stock|cell names for the cells you describe, for use in `events`.
#   cells_where("S", age = "15-64")                    # every S cell in that age group
#   cells_where("S", age = "65+", loc = c("A", "B"))   # filter on several dimensions
cells_where <- function(stock, ...) {
  filt    <- list(...)
  unknown <- setdiff(names(filt), names(dn))
  if (length(unknown) > 0) stop("Unknown dimension(s): ", paste(unknown, collapse = ", "))
  grid <- expand.grid(dn, stringsAsFactors = FALSE)
  keep <- rep(TRUE, nrow(grid))
  for (d in names(filt)) keep <- keep & grid[[d]] %in% filt[[d]]
  paste(stock, cell_labels()[keep], sep = "|")
}

# Structured mixing, one contact matrix per dimension (optional; see the model below).
# mix_dim() applies matrix M along dimension k of array x; mix_all() applies one matrix
# per dimension in turn. The result has the same shape as the input. Assumes each
# dimension mixes independently of the others.
mix_dim <- function(x, M, k) {
  d    <- dim(x)
  perm <- c(k, setdiff(seq_along(d), k))          # bring dimension k to the front
  xm   <- M %*% matrix(aperm(x, perm), nrow = d[k])
  aperm(array(xm, dim = d[perm]), order(perm))    # put the dimensions back in order
}
mix_all <- function(P, mats) {
  for (k in seq_along(mats)) P <- mix_dim(P, mats[[k]], k)
  P
}


## ---- 5. Model setup ---- ##

# Parameters (the Stella constants)
# All parameters used in the model must be defined here. Use lower case short
# names. Each is a single value, or an array built with along() (or array()) that
# is shaped like the dimensions. It doesn't matter if one is defined but not used.
parameters <- list(                                              # <<< EDIT
  beta  = 0.4,                          # transmission rate, per day
  gamma = 0.2                           # recovery rate, per day (1 / infectious period)
  # sus = along(c(1.0, 0.8, 1.2), 1)    # example: relative susceptibility by age
)

# Contact matrices for structured mixing: one per dimension, in the same order as dn
# (rows and columns follow that dimension's levels). Used by mix_all() in the model.
# contact <- list(matrix(c(8, 3, 1,
#                          3, 6, 2,
#                          1, 2, 3), 3, 3, byrow = TRUE))            # <<< EDIT

# Lookup tables (the Stella graphical functions)
# rule = 2 holds the end values flat outside the range, as Stella does.
# effect_fn <- approxfun(x = c(0, 1, 2), y = c(0, 1, 1.3), rule = 2)   # <<< EDIT

# Initial states of the stocks. Each stock is an array shaped like the dimensions.
# Names must match those used in sd_mod below, and be in the same order.
N0 <- along(c(20000, 60000, 20000), 1)       # starting population  <<< EDIT
I0 <- 0 * N0
I0["15-64"] <- 10                            # seed cases by name, e.g. I0["15-64", , "A", ] <- 5

init <- list(
  S = N0 - I0,
  I = I0,
  R = 0 * N0
)

# Flat state vector for deSolve, with a label on every element
y0 <- unlist(lapply(init, as.vector))
names(y0) <- paste(rep(names(init), each = n_cells),
                   rep(cell_labels(), times = length(init)), sep = "|")

# Turn the flat vector back into one array per stock (used inside sd_mod)
unpack <- function(y) {
  stocks <- names(init)
  setNames(lapply(seq_along(stocks), function(i)
    array(y[(i - 1) * n_cells + seq_len(n_cells)], dim = dims, dimnames = dn)),
    stocks)
}

# Time frame
# Must be numeric and in increasing order
times <- seq(0, 365, by = 1)                                     # <<< EDIT


## ---- 6. Validation ---- ##
# times is finite and increasing
if (!is.numeric(times) || any(!is.finite(times)) || !all(diff(times) > 0)) {
  stop("'times' must be a numeric vector of finite values in increasing order.")
}
# dimension names are unique, every dimension has levels, and no level contains "|"
if (length(unique(names(dn))) != length(dn) || is.null(names(dn)) || any(lengths(dn) < 1)) {
  stop("Dimension names in 'dn' must be present and unique, and each dimension needs at least one level.")
}
if (any(grepl("|", unlist(dn), fixed = TRUE))) {
  stop("Dimension level names must not contain the character '|'.")
}
# parameters is a fully named list, with unique names
if (!is.list(parameters) || is.null(names(parameters)) || any(names(parameters) == "")) {
  stop("'parameters' must be a list with a name for every element.")
}
if (length(unique(names(parameters))) != length(parameters)) {
  stop("Parameter names in 'parameters' must be unique.")
}
# each parameter is a single value or an array shaped like the dimensions
# (a wrongly shaped one would otherwise be silently recycled)
bad_shape <- !vapply(parameters, function(p) is.list(p) || length(p) == 1 || is_cell_array(p),
                     logical(1))
if (any(bad_shape)) {
  stop("Parameter(s) with the wrong shape (need one value, or an array shaped like the dimensions): ",
       paste(names(parameters)[bad_shape], collapse = ", "))
}
# stock names are unique
if (length(unique(names(init))) != length(init) || is.null(names(init))) {
  stop("Stock names in 'init' must be present and unique.")
}
# names are disjoint between stocks and parameters
if (length(intersect(names(init), names(parameters))) > 0) {
  stop("Stock names and parameter names must be disjoint (no overlap).")
}
# every stock is an array shaped like the dimensions, with finite values
if (!all(vapply(init, is_cell_array, logical(1)))) {
  stop("Every stock in 'init' must be an array shaped like the dimensions (build it with along() or array()).")
}
if (any(!is.finite(y0)) || (non_negative && any(y0 < 0))) {
  stop("Initial stock values must be finite", if (non_negative) " and not negative.", ".")
}


## ---- 7. The model: auxiliaries, flows, then derivatives ---- ##
# deSolve calls this at every time step. Keep to the same order each time:
#   a) unpack the stocks   b) auxiliaries   c) flows   d) derivatives
# A stock's derivative = inflows - outflows. Return them in the same order as `init`.
sd_mod <- function(time, state, parameters) {
  with(parameters, {

    # a) stocks
    st <- unpack(state)
    S <- st$S; I <- st$I; R <- st$R

    # b) auxiliaries (variables worked out from the stocks)             <<< EDIT
    N          <- S + I + R
    prevalence <- sum(I) / sum(N)           # everyone mixes equally
    lambda     <- beta * prevalence         # force of infection (per person per day)
    # Structured mixing instead (define `contact` in the Model setup):
    # lambda <- beta * mix_all(safe_div(I, N), contact)    # one value per cell
    # beta then no longer equals "transmissions per contact" for a given R0, so it
    # needs scaling to hit a target R0. See sd_array_models.R (scale_beta).

    # c) flows (rates, per time unit)                                   <<< EDIT
    infection  <- S * lambda
    recovery   <- gamma * I

    # d) derivatives, in the same order as `init`                       <<< EDIT
    # Each is the net flow into the stock: positive values increase it,
    # negative values decrease it.
    dS <- -infection
    dI <-  infection - recovery
    dR <-  recovery

    # The first element is the derivatives; anything after is an extra output
    # column (handy for plotting an auxiliary). Named, scalar values only.
    return(list(c(dS, dI, dR), prevalence = prevalence))
  })
}


## ---- 8. Processing ---- ##

# Interventions at set times ("events"): change stocks at a time.
# The var names must match names(y0): one row per cell, so use cells_where().
# method = "multiply" scales the stock, "replace" sets it, "add" adds to it.
# Leave as NULL if you have none.
events <- NULL
# events <- list(data = data.frame(var = cells_where("S", age = "15-64"),
#                                  time = 100, value = 0.9, method = "multiply"))
if (!is.null(events) && !all(events$data$var %in% names(y0))) {
  stop("Some events refer to cells that do not exist: ",
       paste(setdiff(events$data$var, names(y0)), collapse = ", "))
}

# A function, so scenarios and sensitivity runs are one-liners
run_model <- function(parameters, method = "lsoda") {
  deSolve::ode(
    y      = y0,           # Starting values for the stocks; updated as the model runs
    times  = times,        # Times at which to return results
    func   = sd_mod,       # Function that calculates the rates of change
    parms  = parameters,   # Parameters supplied to sd_mod on each call
    method = method,       # "lsoda" is adaptive and handles stiff systems
    events = events
  )
}
# method = "euler" uses fixed steps at the spacing of `times`. To match Stella's
# Euler, use a `times` grid with the same DT and add a non-negative guard.

output <- run_model(parameters)


## ---- 9. Checks (do these every time) ---- ##
state_cols <- grep("\\|", colnames(output))   # the stock columns (names contain "|")

# Tiny negative values (around 1e-8 or smaller) are solver noise, not a model fault.
if (non_negative) {
  cat("smallest stock value:", signif(min(output[, state_cols]), 3), "\n")
}

# For a closed model the total across stocks should stay at its starting value.
# (Not meaningful if events change stocks, or if stocks can't be added together.)
if (check_conserved) {
  total <- rowSums(output[, state_cols])
  cat("max drift in total:", signif(max(abs(total - total[1])), 3), "\n")
}

# Better still: if a simpler version of your model has an exact solution, compare
# against it (for example by collapsing all the dimensions to one group).


## ---- 10. Outputs ---- ##

# Convert to long format: one row per time x stock x cell, one column per dimension
output_long <- as.data.frame(output) |>
  select(time, contains("|")) |>
  pivot_longer(-time, names_to = "cell", values_to = "Value") |>
  separate_wider_delim(cell, delim = "|", names = c("Stock", names(dn))) |>
  mutate(Stock = factor(Stock, levels = names(init)),
         across(all_of(names(dn)), ~ factor(.x, levels = dn[[cur_column()]])))

# ARRAYSUM equivalent: total over every dimension except the ones you group by
totals <- output_long |> group_by(time, Stock) |> summarise(Value = sum(Value), .groups = "drop")

# Plot using ggplot2, one panel per level of one dimension
plot_by       <- names(dn)[1]                                    # <<< EDIT
stock_palette <- setNames(rep_len(base_cols, length(init)), names(init))

p <- output_long |>
  group_by(time, Stock, across(all_of(plot_by))) |>
  summarise(Value = sum(Value), .groups = "drop") |>
  ggplot(aes(x = time, y = Value, colour = Stock)) +
  geom_line(linewidth = my_linewidth) +
  facet_wrap(plot_by) +
  scale_colour_manual(values = stock_palette) +
  labs(title = plot_title,
       x = paste0("Time (", time_unit, ")"),
       y = value_unit) +
  theme_minimal()

print(p)


## ---- 11. Optional: scenarios and sensitivity ---- ##
# Because run_model() takes the parameter list, a scenario is a copy with one change:
# scenario <- modifyList(parameters, list(beta = 0.3))
# output_b <- run_model(scenario)
#
# Sensitivity over a range:
# betas <- seq(0.2, 0.6, by = 0.1)
# peaks <- sapply(betas, function(b) {
#   o <- run_model(modifyList(parameters, list(beta = b)))
#   max(rowSums(o[, grep("^I\\|", colnames(o))]))
# })
# plot(betas, peaks, type = "b", xlab = "beta", ylab = "Peak infectious")
