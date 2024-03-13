rm(list = ls())
# -------------------------------------------------------------------------
# Load libraries:
require(statmod)
require(jtools)
require(interactions)
require(DHARMa)
require(ggplot2)
require(viridis)
require(dplyr)
require(tidyr)
require(mgcv)
require(sf)
require(mgcViz)
require(itsadug)
theme_set(theme_classic())
source('aux_functions.R')

# -------------------------------------------------------------------------

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-FS'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images/gamst'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/model_outputs/gamst'
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
# Hurdle model:

# Define dataset:
model_df = joinDF %>%  mutate(presence = catch > 0)
st_geometry(model_df) = NULL

# Run model 1:
gamst_mod_1 = bam(presence ~ year:quarter + s(lon, lat, bs = "re", by = year), data = model_df, family = binomial)
summary(gamst_mod_1)
save(gamst_mod_1, file = file.path(model_folder, 'my_model_1.RData'))

# Check residuals:
res_model = simulateResiduals(gamst_mod_1)
jpeg(filename = file.path(plot_folder, 'residuals_check_1.jpg'), width = 190, height = 100, units = 'mm', res = 500)
par(mar = c(4, 4, 1, 0.5))
plot(res_model, title = NULL)
dev.off()

# itsadug::pvisgam(gamst_mod_1, view=c("lon", "lat"), select=5)

# Run model 2:
gamst_mod_2 = bam(log(catch) ~ year:quarter + s(lon, lat, bs = "re", by = year), data = subset(model_df, catch > 0))
summary(gamst_mod_2)
save(gamst_mod_2, file = file.path(model_folder, 'my_model_2.RData'))

# Check residuals:
res_model = simulateResiduals(gamst_mod_2)
jpeg(filename = file.path(plot_folder, 'residuals_check_2.jpg'), width = 190, height = 100, units = 'mm', res = 500)
par(mar = c(4, 4, 1, 0.5))
plot(res_model, title = NULL)
dev.off()

# itsadug::pvisgam(gamst_mod_2, view=c("lon", "lat"), select=6)

# Predict CPUE (grid)
varData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), ID = unique(model_df$ID)) 
                      # h_sunrise = mean(model_df$h_sunrise), 
                      # den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow')
MyGrid2 = st_centroid(MyGrid) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
predData = full_join(varData, MyGrid2[,c('ID', 'lon', 'lat')], by = 'ID')
predData$geometry = NULL
predData$pred_1 = predict(gamst_mod_1, newdata = predData, type = 'response', allow.new.levels = TRUE)
predData$pred_2 = exp(predict(gamst_mod_2, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData$cpue_pred = predData$pred_1 * predData$pred_2
PredGrid = left_join(MyGrid, predData, by = c('ID'))
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

# Calculate index (weighted sum by grid)
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(value = sum(cpue_pred*grid_area*portion_on_ocean))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4, type = 'GAMst')
st_geometry(PredTime) = NULL
write.csv(PredTime, file = file.path(model_folder, 'CPUE_ts.csv'), row.names = FALSE)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value*1e-08)) +
  geom_line() +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_gamst.jpg'), width = 95, height = 70, units = 'mm', dpi = 500)

