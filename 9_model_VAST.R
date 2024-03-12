rm(list = ls())

# -------------------------------------------------------------------------
library(VAST)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
theme_set(theme_classic())

# -------------------------------------------------------------------------
data_folder = 'data'
plot_folder = 'plots/vast'
params = list(species = "SKJ", ORP = "IOTC")

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Map information for plotting:
limites = read.table(file.path(data_folder, "limites.csv"), header=TRUE,sep=",", na.strings="NA", dec=".", strip.white=TRUE)
limites = subset(limites, ORP == params$ORP)
xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
worldmap = map_data("world")
colnames(worldmap) = c("X", "Y", "PID", "POS", "region", "subregion")
yBreaks = seq(from = -20, to = 20, by = 20)
xBreaks = seq(from = 40, to = 100, by = 20)

# -------------------------------------------------------------------------
# Define dataset:
model_df = joinDF %>% mutate(presence = catch > 0, area_swept = 1,
                        time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                        time2 = (model_df$time - min(model_df$time))*4 + 1)
st_geometry(model_df) = NULL

settings = make_settings(n_x = 200, Region = 'User',
                          purpose = "index2", bias.correct=FALSE,
                          knot_method = 'grid')
settings$FieldConfig[2,] = 0 ## turn off temporal components

user_region = readRDS('user_region.rds')
vast_mod_1 = fit_model(settings = settings,
                       Lat_i = model_df$lat, Lon_i = model_df$lon,
                       t_i = model_df$time2, b_i = model_df$catch,
                       a_i = model_df$area_swept,
                       input_grid = user_region, run_model = TRUE)
plot_results(vast_mod_1, plot_set=3)




plot(vast_mod_1$spatial_list$MeshList$anisotropic_mesh)
plot(vast_mod_1$spatial_list$loc_i[,2], vast_mod_1$spatial_list$loc_i[,1])

