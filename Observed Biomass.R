pkgs = c("tidyverse", "mapplots", "ggpubr", "RColorBrewer", "dplyr")
for(p in pkgs){
  if(!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
} 


load('South_west_survey_data.Rdata') #Load survey dataset
load('landings_data_processed.Rdata') #Loading landing dataset


#Choose our species of interest
spp_list = c('Atlantic cod', 'European plaice', 'Common sole', 'Haddock', 'Whiting', 'Anglerfishes nei', 'Saithe(=Pollock)', 'Ling', 'Megrims nei', 'European hake')
tax_list = unique(subset(land_df, English_name %in% spp_list)$new_taxa) %>% c('Lophius piscatorius', 'Lepidorhombus whiffiagonis') #Pulls out the taxonomic names which appear in both SW Survey data and Celtic sea landings data

#Average observed biomass across the years
survey.avg = surv_df %>%
  filter(taxa %in% tax_list,
         gear == 'beam') %>% 
  #Sum by haul first
  group_by(taxa, year, HaulID) %>%
  summarise(kgkm = sum(DensBiom_kg_Sqkm)) %>%
  filter(!is.na(kgkm)) %>%
  #Take mean over time
  group_by(taxa) %>%
  summarise(av_kgkm=mean(kgkm)) %>%
  ungroup()
survey.avg

#Average biomass over years
survey.trends = surv_df %>%
  filter(taxa %in% tax_list) %>% 
  #Sum by haul first
  group_by(HaulID,  year, taxa) %>%
  summarise(kgkm = sum(DensBiom_kg_Sqkm)) %>%
  #Take mean over time
  group_by(taxa, year) %>%
  summarise(av_kgkm=mean(kgkm)) %>%
  ungroup()
survey.trends

#Plot average biomass time series
p<- ggplot(survey.trends, aes(x=year, y=av_kgkm)) +
  geom_line(color="grey") +
  geom_point(color="blue") +
  facet_wrap(~taxa, scales="free") + 
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  labs(title = "Mean observed biomass",
       x = "Year",
       y = "Average biomass density (kg/km^2)")
p

### USE GEOMETRIC MEAN TO MINIMISE THE DOMINANCE OF OUTLIERS FOR THE ARITHMETIC MEAN ###

#Average observed biomass across the years
survey.geoavg = surv_df %>%
  filter(taxa %in% tax_list,
         gear == 'beam') %>% 
  #Sum by haul first
  group_by(taxa, year, HaulID) %>%
  summarise(kgkm = sum(DensBiom_kg_Sqkm)) %>%
  filter(!is.na(kgkm)) %>%
  #Take mean over time
  group_by(taxa) %>%
  summarise(av_kgkm=exp(mean(log(kgkm)))) %>%
  ungroup()
survey.geoavg

#Average biomass over years
survey.geotrends = surv_df %>%
  filter(taxa %in% tax_list) %>% 
  #Sum by haul first
  group_by(HaulID,  year, taxa) %>%
  summarise(kgkm = sum(DensBiom_kg_Sqkm)) %>%
  #Take mean over time
  group_by(taxa, year) %>%
  summarise(av_kgkm=exp(mean(log(kgkm)))) %>%
  ungroup()
survey.geotrends

#Plot average biomass time series
p2<- ggplot(survey.geotrends, aes(x=year, y=av_kgkm)) +
  geom_line(color="grey") +
  geom_point(color="blue") +
  facet_wrap(~taxa, scales="free") + 
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  labs(title = "Mean observed biomass",
       x = "Year",
       y = "(Geometric) Average biomass density (kg/km^2)")
p2
