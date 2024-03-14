rm(list = ls())
# -------------------------------------------------------------------------
library(sdmTMB)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
theme_set(theme_classic())
source('aux_functions.R')


# -------------------------------------------------------------------------
data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-LS'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS-LS_YFT/images/sdmTMB'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS-LS_YFT/model_outputs/sdmTMB'
# create folders in case missing:
dir.create(plot_folder, showWarnings = FALSE)
dir.create(model_folder, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# Define dataset:
model_df = joinDF %>% 
              mutate(yyyyqq = paste(year, quarter, sep = '-'),
                     time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
# model_df = model_df %>% filter(time < 2015) 

# Make some changes before running model:
st_geometry(model_df) = NULL
model_df = sdmTMB::add_utm_columns(dat = model_df, ll_names = c('lon', 'lat'))

# Make mesh:
mesh = make_mesh(model_df, c("X", "Y"), cutoff = 60)
jpeg(filename = file.path(plot_folder, 'mesh.jpg'), width = 95, height = 100, units = 'mm', res = 500)
par(mar = c(1, 1, 1, 1))
plot(mesh)
dev.off()

# Spatial model:
my_st_model <- sdmTMB(
  data = model_df,
  formula = catch ~ 0 + as.factor(time),
  mesh = mesh, 
  family = delta_lognormal(),
  time = 'time',
  spatial = "on"
  #spatiotemporal = 'ar1'
)

# Check residuals (only fixed effects):

# TODO


# -------------------------------------------------------------------------
# Predictions:
varData = expand.grid(year = as.factor(levels(my_st_model$data$year)), 
                      quarter = as.factor(levels(my_st_model$data$quarter)), 
                      ID = unique(my_st_model$data$ID))
varData = varData %>% mutate(yyyyqq = paste(year, quarter, sep = '-'),
                             time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
predData = full_join(varData, MyGrid2[,c('ID', 'lon', 'lat')], by = 'ID')
predData$geometry = NULL
predData = sdmTMB::add_utm_columns(dat = predData, ll_names = c('lon', 'lat'))
predData = predData[complete.cases(predData), ] # check this later, not sure why this happens (some NAs)
predictions = predict(my_st_model, newdata = predData, return_tmb_object = TRUE)
predictions$data$cpue_pred = boot::inv.logit(predictions$data$est1)*exp(predictions$data$est2)
plot_data = predictions$data
PredGrid = left_join(MyGrid, plot_data, by = c('ID'))
plot_data = PredGrid %>% group_by(year, ID) %>% summarise(cpue_pred = mean(cpue_pred)) # aggregate by year and cluster, just to display in a map

p1 = ggplot() +  
  geom_sf(data = plot_data, aes(fill = cpue_pred, color = cpue_pred)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ year) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_predictions_gamst.jpg'), plot = p1, 
       width = 190, height = 180, units = 'mm', dpi = 500)


# Calculate index (weighted sum by area): or average?
PredTime = PredGrid
st_geometry(PredTime) = NULL
PredTime = PredTime[complete.cases(PredTime), ] # check this later, not sure why this happens (some NAs)
index = get_index(predictions, area = 4, bias_correct = TRUE)



# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
  geom_line() +
  labs(color = 'CPUE type') +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_inla_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

