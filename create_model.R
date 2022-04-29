#Load necessary R packages
library(devtools)
install_github("sizespectrum/mizerExperimental")
install_github("sizespectrum/mizer")
library(dplyr)
library(mizer)
library(mizerExperimental) #Has some additional features to mizer R package
library(ggplot2)
library(readxl)

#Define function that will be used to calculate Balanced Harvesting fishing mortality
productionFMort <- function(params, n, effort, e_growth, ...) {
  usual_f_mort <- colSums(mizerFMortGear(params, effort))
  production_density <- sweep(e_growth * n, 2, params@w, "*")
  production_density * usual_f_mort
}

#Define function that adds up yield of all species at each time step
TotalYield <- function(sim) {
  y <- getYield(sim)
  y_total <- rowSums(y)
  return(y_total)
}

#Define function to add up biomass of all species at each time step
TotalBiomass <- function(sim) {
  b <- getBiomass(sim)
  b_total <- rowSums(b)
  return(b_total)
}

#Define function the adds up abundance spectrum of all species at each time step
TotalAbundance <- function(sim) {
  n <- getN(sim)
  n_total <- rowSums(n)
  return(n_total)
}

Model_params <- read_excel("Model params.xlsx", range = "B2:J12") #Load compiled parameter list into R

sp <- Model_params #Creates a data frame of the original model parameters

#Need to rename the variables in the data frame to coincide with how mizer recognises them
sp$species <- sp$'Species of interest'
sp$w_inf <- sp$W_inf
sp$biomass_observed <- sp$'Observed Biomass (kg/km^2)'
sp$k_vb <- sp$'Growth coefficient K'

#Allometric parameters a and b
sp$a <- 0.006
sp$b <- 3
sp$w_mat <- sp$a * sp$L_mat ^ sp$b #Convert length at maturity to weight at maturity
sp$w_mat

#Adding predation parameters
Predation_Params <- read_excel("Predation Params.xlsx") #Load predation parameter list into R
sp$beta <- Predation_Params$beta
sp$sigma <- Predation_Params$sigma

params <- newMultispeciesParams(sp) #Create parameter object

#Adding fishing gear parameters
gear_params(params)
gear_params(params)$catchability <- sp$`Fishing (external) mortality rate`
gear_params(params)$knife_edge_size <- as.numeric(sp$`Minimum landing size (cm)`)
gear_params(params)$knife_edge_size[7] <- 29.6 #Choose a value for anglerfish because it has no min landing size

#Adding interaction matrix 
Interaction.matrix <- read.delim("Interaction matrix.txt", header=FALSE) #Load interaction matrix into R
inter <- as.matrix(Interaction.matrix)
rownames(inter) <- sp$`Species of interest`
colnames(inter) <- sp$`Species of interest`
params <- setInteraction(params, interaction = inter) #Pass interaction matrix to parameter object

###Calibrate model to current ecosystem with constant fishing effort at maturity size###
params_minlandfishing <- params

initial_effort(params_minlandfishing) <- 10 #fishing effort

params_minlandfishing <- steady(params_minlandfishing) #Running dynamics to steady state 

species_params(params) #Outputs table of all the species parameters

plotSpectra(params_minlandfishing, power=0, resource = FALSE, total = TRUE) #Plots abundance spectrum as a function of size
plotSpectra(params_minlandfishing, power=1, resource = FALSE, total = TRUE) #Plots biomass spectrum as a function of size or abundance spectrum as a function of log size
plotSpectra(params_minlandfishing, power=2, resource = FALSE, total = TRUE) #Plots biomass as a function of log size aka Sheldon spectrum

plotBiomassObservedVsModel(params_minlandfishing) #Plots model biomasses against observed biomasses
plotGrowthCurves(params_minlandfishing, species_panel = TRUE) #Plots model growth curves against vonB growth curves

#Opens R shiny gadget to tune growth curves 
params_minlandfishing <- tuneGrowth(params_minlandfishing)

#Save updated parameters
saveParams(params_minlandfishing, "params_minlandfishing.rds") #Can be accessed via Git repository

#To re-create following simulation & plots, load in updated parameters object
params_minlandfishing <- readRDS("params_minlandfishing.rds")

sim1 <- project(params_minlandfishing, effort=10, t_max = 75) #Run simulation for 75 years
plot(sim1)
plotBiomass(sim1, total = TRUE)
plotSpectra(sim1, power=1, resource = FALSE, total = TRUE)
slope1 <- getCommunitySlope(sim1, min_w = 0.01, max_w = 1000)
mean(slope1$slope) #-0.4367606

###No fishing###
params_nofishing <- params_minlandfishing
params_nofishing <- steady(params_nofishing) #Re-run dynamics to steady state
#zero fishing effort means no fishing in the system
sim0 <- project(params_nofishing, effort = 0, t_max = 75) #Run simulation for 75 years
plot(sim0)
plotBiomass(sim0, total = TRUE)
plotSpectra(sim0, power=1, resource = FALSE, total = TRUE)
slope0 <- getCommunitySlope(sim0, min_w = 0.01, max_w = 1000)
mean(slope0$slope) #-0.8117456

#Plot relative total abundances by size in the final time step
total_abund0 <- colSums(finalN(sim0)) 
total_abund1 <- colSums(finalN(sim1))
relative_abundance <- total_abund1/total_abund0
plot(x = w(params), y = relative_abundance, log = "xy", type = "n", 
     xlab = "Size [g]", ylab = "Relative abundance")
lines(x = w(params), y = relative_abundance)
lines(x = c(min(w(params)), max(w(params))), y = c(1, 1), lty = 2)

#Plot total predation mortality by size for each fishing scenario
PMort0 <- getPredMort(params_nofishing, finalN(sim0))[4,]
PMort1 <- getPredMort(params_minlandfishing, finalN(sim1))[4,]
plot(x = w(params), y = PMort0, log = "x", type = "n", 
     xlab = "Size [g]", ylab = "Predation Mortality [1/year]")
lines(x = w(params), y = PMort0, lty = 2)
lines(x = w(params), y = PMort1)
legend(x="topright", legend = c("no fishing", "fishing"), lty = c(2,1))

###Balanced harvesting###
params_BH <- params_minlandfishing
#Change size selectivity and gear to describe balanced harvesting
params_BH@gear_params$knife_edge_size <- params_BH@species_params$Minimum.landing.size..cm.
params_BH@gear_params$knife_edge_size[7] <- 29.6
params_BH@species_params$knife_edge_size <- params_BH@species_params$Minimum.landing.size..cm.
params_BH@species_params$knife_edge_size[7] <- 29.6
gear_params(params_BH)$gear <- "balanced"
species_params(params_BH)$gear <- "balanced"
params_BH <- setRateFunction(params_BH, "FMort", "productionFMort")
params_BH <- setFishing(params_BH, initial_effort = 10e4)
saveParams(params_BH, "params_BH_beforesteady.rds")
params_BH <- steady(params_BH) #Run dynamics to steady state
params_BH <- tuneGrowth(params_BH) #Tune parameters then re-run to steady state

#Save updated parameters
saveParams(params_BH, "params_BH.rds") #Can be accessed via Git repository

#To re-create following simulation & plots, load in updated parameters object
params_BH <- readRDS("params_BH.rds")

simBH <- project(params_BH, effort=10e4, t_max = 75) #Run simulation for 75 years
plot(simBH)
plotBiomass(simBH, total = TRUE)
plotSpectra(simBH, power=1, resource = FALSE, total = TRUE)
slopeBH <- getCommunitySlope(simBH, min_w = 0.01, max_w = 1000)
mean(slopeBH$slope) #-1.291445

#Simulate a BH model without running to steady state 
#To re-create following simulation & plots, load in updated parameters object
params_BH2 <- readRDS("params_BH_beforesteady.rds")
params_BH2 <- setFishing(params_BH2, initial_effort = 0.1)
params_BH2@gear_params$catchability <- 1
simBH2 <- project(params_BH2, effort=0.1, t_max = 75) #Run simulation for 75 years
plot(simBH2)
plotBiomass(simBH2, total = TRUE)
plotSpectra(simBH2, power=1, resource = FALSE, total = TRUE)
slopeBH2 <- getCommunitySlope(simBH2, min_w = 0.01, max_w = 1000)
mean(slopeBH2$slope) #-1.039897

#Plot total yield for each fishing method over time
TotYield1 <- TotalYield(sim1)
TotYieldBH <- TotalYield(simBH)
matplot(cbind(TotYield1, TotYieldBH), 
        log = "xy", type = "l", #main = "Total yield of all species over time",  
        xlab = "Year", ylab = "Yield", col = c(2,4), lty=c(1,1), 
        lwd = 3, bty = "l",cex.axis=1.2,cex.lab=1.2)
legend(x = "center",legend = c("Constant Fishing at Min Landing Size","Balanced Harvesting"),
       lty = c(1,1), col = c(2,4), lwd = 3, cex = 1.2, title = "Fishing method")

TotYieldBH2 <- TotalYield(simBH2)
matplot(cbind(TotYield1, TotYieldBH2), 
        log = "xy", type = "l", #main = "Total yield of all species over time",  
        xlab = "Year", ylab = "Yield", col = c(2,4), lty=c(1,1), 
        lwd = 3, bty = "l",cex.axis=1.2,cex.lab=1.2)
legend(x = "right",legend = c("Constant Fishing at Min Landing Size","Balanced Harvesting II"),
       lty = c(1,1), col = c(2,4), lwd = 3, cex = 1, title = "Fishing method")

#Plot total biomass for each fishing method over time
TotBiom1 <- TotalBiomass(sim1)
TotBiom0 <- TotalBiomass(sim0)
TotBiomBH <- TotalBiomass(simBH)
TotBiomBH2 <- TotalBiomass(simBH2)
matplot(cbind(TotBiom1, TotBiom0, TotBiomBH, TotAbunBH2), 
        log = "xy", type = "l", main = "Total Biomass",  
        xlab = "Year", ylab = "Biomass [kg km^-2 year^-1]", col = c(2,3,4,1), lty=c(1,1,1), 
        lwd = 3, cex.main=1.2, cex.axis=1.5, cex.lab=1.5)
legend(x = "topleft", legend = c("Constant Fishing at Min Landing Size","No Fishing","Balanced Harvesting","Balanced Harvesting II"),
       lty = c(1,1,1), col = c(2,3,4,1), lwd = 3, cex = 0.7, title = "Fishing method")
#Just plotting the two fishing systems
matplot(cbind(TotBiom1, TotBiomBH), 
        log = "xy", type = "l", main = "Total Biomass",  
        xlab = "Year", ylab = "Biomass [kg km^-2 year^-1]", col = c(2,4), lty=c(1,1), 
        lwd = 3, cex.main=1.2, cex.axis=1.5, cex.lab=1.5)
legend(x = "center",legend = c("Constant Fishing at Min Landing Size","Balanced Harvesting"),
       lty = c(1,1), col = c(2,4), lwd = 3, cex = 1.2, title = "Fishing method")

matplot(cbind(TotBiom1, TotBiomBH2), 
        log = "xy", type = "l", main = "Total Biomass",  
        xlab = "Year", ylab = "Biomass [kg km^-2 year^-1]", col = c(2,4), lty=c(1,1), 
        lwd = 3, cex.main=1.2, cex.axis=1.5, cex.lab=1.5)
legend(x = "right",legend = c("Constant Fishing at Min Landing Size","Balanced Harvesting II"),
       lty = c(1,1), col = c(2,4), lwd = 3, cex = 1, title = "Fishing method")

#Plot total abundance for each fishing method over time
TotAbun1 <- TotalAbundance(sim1)
TotAbun0 <- TotalAbundance(sim0)
TotAbunBH <- TotalAbundance(simBH)
TotAbunBH2 <- TotalAbundance(simBH2)
matplot(cbind(TotAbun1, TotAbun0, TotAbunBH, TotAbunBH2), 
        log = "xy", type = "l", main = "Total Abundance",  
        xlab = "Year", ylab = "Abundance", col = c(2,3,4,1), lty=c(1,1,1), 
        lwd = 3, cex.main=1.2, cex.axis=1.5, cex.lab=1.5)
legend(x = "topleft",legend = c("Constant Fishing at Min Landing Size","No Fishing","Balanced Harvesting", "Balanced Harvesting II"),
       lty = c(1,1,1), col = c(2,3,4,1), lwd = 3, cex = 0.7, title = "Fishing method")
mean(TotAbunBH) #47744.72
mean(TotAbunBH2) #49755.26
min(TotAbun0) #219.0246
max(TotAbun0) #10263464886
