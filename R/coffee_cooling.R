## ---------------------------
## Script name: coffee_cooling.R
## Purpose of script: Simulate how fast a cup of coffee cools to room temperature
##                    from different starting temperatures (Newton's law of cooling)
## Project name:
## Project reference number:

## Author: ADD YOUR NAME
## Date created: ADD DATE
## Date last updated:
## Copyright (c) Kent Public Health Observatory, Kent County Council, ADD YEAR
## Email: kpho@kent.gov.uk

## ---------------------------
## Notes:
## One stock (C, coffee temperature). One flow (cooling), proportional to the
## gap between the coffee and the room:   dC/dt = -k * (C - t_room)
## This has an exact solution, so the script checks itself against it:
##   C(t) = t_room + (C0 - t_room) * exp(-k * t)
## Built from sdm_template_simple.R.
## ---------------------------


## ---- 1. Packages ---- ##
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)


## ---- 2. Settings ---- ##
value_unit <- "°C"                   # degrees Celsius (escape keeps the file ASCII)
time_unit  <- "minutes"

# Temperature can legitimately be negative, so don't apply the "not negative" rule
non_negative <- FALSE

base_cols    <- c("#2D7C82", "#A9A9A9", "#D98E04", "#8C4A8C", "#B5413B")
my_linewidth <- 1.5
plot_title   <- "Cooling coffee"


## ---- 4. Model setup ---- ##

parameters <- c(
  k      = 0.1,   # cooling rate, per minute (time constant 1/k = 10 minutes)
  t_room = 18     # room temperature, degrees C, assumed constant
)

# Starting temperatures are defined ONCE here, and the scenarios are built from them
start_temps <- c(100, 80, 60)
names(start_temps) <- paste0(start_temps, "°C")

init <- c(C = start_temps[[1]])           # base value; each scenario overrides it

times <- seq(0, 60, by = 1)               # about 6 time constants, so near room temperature at the end

scenarios <- lapply(start_temps, function(x) list(init = c(C = x)))


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
for (nm in names(scenarios)) {
  bad <- c(setdiff(names(scenarios[[nm]]$parameters), names(parameters)),
           setdiff(names(scenarios[[nm]]$init), names(init)))
  if (length(bad) > 0) {
    stop("Scenario '", nm, "' overrides unknown name(s): ", paste(bad, collapse = ", "))
  }
}


## ---- 6. The model ---- ##
sd_mod <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

    # a) auxiliaries
    discrepancy <- C - t_room       # how much hotter the coffee is than the room

    # b) flows
    cooling <- k * discrepancy      # degrees lost per minute. Negative if the
                                    # coffee is colder than the room (it warms up)

    # c) derivatives
    dC <- -cooling

    return(list(c(dC)))
  })
}


## ---- 7. Run every scenario ---- ##
override <- function(x, new) {
  if (length(new) > 0) x[names(new)] <- new
  x
}

run_scenario <- function(sc, method = "lsoda") {
  deSolve::ode(
    y      = override(init, sc$init),
    times  = times,
    func   = sd_mod,
    parms  = override(parameters, sc$parameters),
    method = method
  )
}

output <- lapply(names(scenarios), function(nm) {
  as.data.frame(run_scenario(scenarios[[nm]])) |> mutate(scenario = nm)
}) |>
  bind_rows() |>
  mutate(scenario = factor(scenario, levels = names(scenarios)))


## ---- 8. Checks ---- ##
# Temperatures are not conserved (the coffee loses heat), so there is no
# "total" to check. Instead, compare with the exact solution.
t_room <- parameters[["t_room"]]
k      <- parameters[["k"]]

c0    <- start_temps[as.character(output$scenario)]
exact <- t_room + (c0 - t_room) * exp(-k * output$time)

cat("max difference from exact solution:", signif(max(abs(output$C - exact)), 3),
    value_unit, "\n")
cat("lowest temperature reached:", signif(min(output$C), 4), value_unit,
    "(should stay above room temperature,", t_room, value_unit, ")\n")


## ---- Key result: minutes to cool to a drinkable temperature ---- ##
threshold <- 40                                                  # <<< EDIT (degrees C)

key <- output |>
  group_by(scenario) |>
  summarise(
    minutes_from_output = approx(C, time, xout = threshold)$y,   # interpolated; NA if never crossed
    minutes_exact       = log((first(C) - t_room) / (threshold - t_room)) / k,
    .groups = "drop"
  )
print(key)


## ---- 9. Outputs ---- ##
output_long <- output |>
  pivot_longer(-c(time, scenario), names_to = "Variable", values_to = "Value") |>
  mutate(Variable = factor(Variable, levels = unique(Variable)))

n_scen       <- nlevels(output_long$scenario)
scen_palette <- setNames(rep_len(base_cols, n_scen), levels(output_long$scenario))

p <- ggplot(output_long, aes(x = time, y = Value, colour = scenario)) +
  geom_hline(yintercept = t_room, linetype = "dashed", colour = "grey50") +
  annotate("text", x = max(times), y = t_room, label = "room temperature",
           hjust = 1, vjust = -0.6, colour = "grey40", size = 3.5) +
  geom_line(linewidth = my_linewidth) +
  facet_null() +                           # only one variable, so no panels needed
  scale_colour_manual(values = scen_palette) +
  labs(title = plot_title,
       x = paste0("Time (", time_unit, ")"),
       y = paste0("Coffee temperature (", value_unit, ")"),
       colour = "Starting\ntemperature") +
  theme_minimal()

print(p)
