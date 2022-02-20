library(devtools)
install_github("sizespectrum/mizerExperimental")
install_github("sizespectrum/mizer")
library(dplyr)
library(mizer)
library(mizerExperimental) #has some additional features to mizer

library(readxl)
Model_params <- read_excel("Model params.xlsx", range = "A1:I11")

sp <- Model_params #creates a data frame of the original model params

#Need to rename the columns in the data frame to fit the arguments of mizer (see cheat sheet names)
sp$species <- sp$'Species of interest'
sp$w_inf <- sp$W_infinity

# Improve these a and b values later
# mizer assumes a = 0.001 and b=3
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

Interaction.matrix <- read.delim("Interaction matrix.txt", header=FALSE)
inter <- as.matrix(Interaction.matrix)
rownames(inter) <- sp$`Species of interest`
colnames(inter) <- sp$`Species of interest`
params <- setInteraction(params, interaction = inter)

initial_effort(params) <- 1 #fishing effort 

params <- steady(params) #takes x years to find a steady state 

species_params(params)

plotSpectra(params)


params <- tuneParams(params, controls = c("abundance", "predation", "reproduction", "other", "interaction", "resource"))
#calibrate and match 
# run the system back to steady state

#saveParams(params, "params.rds")