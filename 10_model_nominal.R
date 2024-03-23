rm(list = ls())
# -------------------------------------------------------------------------
# Load libraries:
require(ggplot2)
require(viridis)
require(dplyr)
require(tidyr)
require(sf)
theme_set(theme_classic())
source('aux_functions.R')

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-LS'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS-LS_YFT/model_outputs'

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
predData = joinDF
st_geometry(predData) = NULL
PredTime = predData %>% group_by(year, quarter) %>% summarise(CPUE_est = mean(catch))
# Include area size in the CPUE-nominal calculation?
# If so, it will be influenced by the covered area by the fleet
# varData = predData %>% group_by(year, quarter, ID) %>% summarise(cpue_pred = mean(catch))
# PredGrid = left_join(MyGrid, varData, by = c('ID'))
# PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(CPUE_est = sum(cpue_pred*grid_area*portion_on_ocean))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4, 
                               CPUE_est = CPUE_est/mean(CPUE_est),
                               type = 'Nominal')
write.csv(PredTime, file = file.path(model_folder, 'CPUE_ts.csv'), row.names = FALSE)