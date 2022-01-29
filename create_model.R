library(devtools)
install_github("sizespectrum/mizerExperimental")
install_github("sizespectrum/mizer")
library(dplyr)
library(mizer)
library(mizerExperimental)

library(readxl)
Model_params <- read_excel("Model params.xlsx", range = "A1:I11")

sp <- Model_params
sp$species <- sp$'Species of interest'
sp$w_inf <- sp$W_infinity

# Improve these a and b values later
sp$a <- 0.006
sp$b <- 3

sp$w_mat <- sp$a * sp$L_mat ^ sp$b
sp$w_mat

sp$biomass_observed <- sp$'Observed Biomass kg/km^2'
sp$k_vb <- sp$'Growth coefficient K'

params <- newMultispeciesParams(sp)

gear_params(params)
gear_params(params)$catchability <- sp$`Fishing (external) mortality rate`
gear_params(params)$knife_edge_size <- as.numeric(sp$`Minimum landing size cm`)
gear_params(params)$knife_edge_size[7] <- 25 #change later

plotSpectra(params)

params <- steady(params)
plotSpectra(params)

params <- tuneParams(params, controls = c("abundance", "predation", "reproduction", "other", "interaction", "resource"))
#calibrate and match 
# run the system back to steady state

#saveParams(params, "params.rds")