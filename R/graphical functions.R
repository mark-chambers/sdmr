## ---------------------------
## Script name: graphical functions.R
## Purpose of script: Demo how to achieve Stella like converters
## Project name: 
## Project reference number: 

## Author: Mark Chambers
## Date created: 2025-09-25
## Date last updated: 
## Copyright (c) Kent Public Health Observatory, Kent County Council, 2025
## Email: kpho@kent.gov.uk

## ---------------------------
## Notes:
##
##
## ---------------------------


## ---- Settings ---- ##

# Filepaths
project_fp <- paste0(Sys.getenv("USERPROFILE"), "/OneDrive - KENT COUNTY COUNCIL/Shared Documents - Kent Analytics Team/Work/KPHO/[project folder]/data/")
main_fp <- "G:/ST Strategic Commissioning/Analytics (R)/Public Health Observatory (R)/data/datasets/"

# Colours
col_kent <- "#2D7C82"
col_eng <- "#A9A9A9"

my_linewidth <- 1.5

start <- 0
finish <- 160
time_step <- 1
timeseries <- seq(start, finish, by = time_step)

# Custom plotting function using base R
plot_gf <- function(x, y, xlab = "Time", ylab = "Value", title = "Line Plot") {
  plot(
    x = timeseries,
    y = y,
    ylim = c(0, max(y)),
    type = "l",
    col = col_kent,
    lwd = my_linewidth,
    xlab = xlab,
    ylab = ylab,
    main = title
  )
}

## ---- Packages ---- ##




## ---- Get data ---- ##




## ---- Processing ---- ##
# Use t as a time variable

# Linear function
linear <- function(t, slope, intercept, scale_output = TRUE) {
  y <- slope * t + intercept
  if(scale_output) {
    y <- y / max(y)  # Scale to max of 1
  }
  return(y)
}
linear_data <- linear(
    t = timeseries,
    slope = -0.5,
    intercept = 90
)
plot_gf(
    y = linear_data,
    title = "Linear function"
)

# Step function
step_change <- function(t, step_time, before, after, scale_output = TRUE) {
  y <- ifelse(t < step_time, before, after)
  if(scale_output) {
    y <- y / max(y)  # Scale to max of 1
  }
  return(y)
}  
step_change_data <- step_change(
    t = timeseries, 
    step_time = 45,
    before = 15,
    after = 50
)
plot_gf(
    y = step_change_data,
    title = "Step change function"
)

# Exponential function
exponential <- function(t, rate, initial, scale_output = TRUE) {
  y <- initial * exp(rate * t)
  if(scale_output) {
    y <- y / max(y)  # Scale to max of 1
  }
  return(y)
} 
exponential_data <- exponential(
    t = timeseries,
    rate = -0.05,  # negative for decay
    initial = 5,
    scale_output = TRUE
)  
plot_gf(
    y = exponential_data,
    title = "Exponential function"
)

# Logistic function
logistic <- function(t, carrying_capacity, growth_rate, midpoint, scale_output = TRUE) {
  y <- carrying_capacity / (1 + exp(-growth_rate * (t - midpoint)))
    if(scale_output) {
        y <- y / max(y)  # Scale to max of 1
    }
  return(y)
}  
logistic_data <- logistic(
    t = timeseries,
    carrying_capacity = 100,
    growth_rate = -0.15,  # negative for decrease
    midpoint = 50
)
plot_gf(
    y = logistic_data,
    title = "Logistic function"
)

# logarithmic function
logarithmic <- function(t, scale, shift, scale_output = TRUE) {
  y <- scale * log(t + shift)
    if(scale_output) {
        y <- y / max(y)  # Scale to max of 1
    }
  return(y)
}   
logarithmic_data <- logarithmic(
    t = timeseries,
    scale = 100,
    shift = 0
)
plot_gf(
    y = logarithmic_data,
    title = "Logarithmic function"
)

# gaussian function
gaussian <- function(t, mean, sd, amplitude) {
  y <- amplitude * exp(-0.5 * ((t - mean) / sd) ^ 2)
  return(y)
} 
gaussian_data <- gaussian(
    t = timeseries,
    mean = median(timeseries),  # center of the peak
    sd = 10,  # controls the width
    amplitude = 100  # peak height
)
plot_gf(
    y = gaussian_data,
    title = "Gaussian function"
)

# sine wave function
sine_wave <- function(t, amplitude, frequency, phase, vertical_shift) {
  y <- amplitude * sin(frequency * t + phase) + vertical_shift
  return(y)
}
sine_wave_data <- sine_wave(
    t = timeseries,
    amplitude = 15,  # height of the wave
    frequency = 12,  # controls how many cycles fit in the time frame
    phase = 0,  # left or right shift
    vertical_shift = 50  # center the wave vertically
)   
plot_gf(
    y = sine_wave_data,
    title = "Sine wave function"
)
sine_wave(t = 0:10, amplitude = 10, frequency = 1, phase = 0, vertical_shift = 50)

## ---- Outputs ---- ##

pdf("graphical functions.pdf", paper = "a4")
par(mfrow = c(1, 1))
plot_gf(y = linear_data, title = "Linear function")
plot_gf(y = step_change_data, title = "Step change function")
plot_gf(y = exponential_data, title = "Exponential function")
plot_gf(y = logistic_data, title = "Logistic function")
plot_gf(y = logarithmic_data, title = "Logarithmic function")
plot_gf(y = gaussian_data, title = "Gaussian function")
plot_gf(y = sine_wave_data, title = "Sine wave function")
dev.off()