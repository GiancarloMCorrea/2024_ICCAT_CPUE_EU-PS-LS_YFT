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
theme_set(theme_classic())
source('aux_functions.R')


# -------------------------------------------------------------------------

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images/gamst'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/model_outputs'

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# Model 1 (using only positive observations): 

# Define dataset:
model_df = joinDF %>% filter(catch > 0) %>%
  mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)

# Run model
gam_st_mod_1 = gam(log(catch) ~ year*quarter + te(lon, lat, time, k = c(12,12,6)) +
                   den_water2 + capacity + h_sunrise + follow_echo,
                 data = model_df)

# Plot and check
par(mfrow = c(2,3))
plot(gam_st_mod_1)
viz_mod = getViz(gam_st_mod_1)
check.gamViz(viz_mod)

# Check residuals:
res_model = simulateResiduals(fittedModel = gam_st_mod_1)
plot(res_model)

# Check effects:
jtools::effect_plot(gam_st_mod_1, pred = time, cat.geom = "line")
jtools::effect_plot(gam_st_mod_1, pred = den_water2, cat.geom = "line")
jtools::effect_plot(gam_st_mod_1, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(gam_st_mod_1, pred = capacity, cat.geom = "line")
jtools::effect_plot(gam_st_mod_1, pred = follow_echo, cat.geom = "line")

# Predict CPUE (grid)
varData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), h_sunrise = mean(model_df$h_sunrise), 
                      den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'),
                      ID = unique(model_df$ID))
varData = varData %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
predData = full_join(varData, extraDF, by = 'ID')
predData$cpue_pred = exp(predict(gam_st_mod_1, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData = predData %>% dplyr::select(-geometry)
PredGrid = left_join(MyGrid, predData, by = c('ID'))
PredGrid = PredGrid %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% na.omit

p1 = ggplot() +  
  geom_sf(data = PredGrid, aes(fill = cpue_pred, color = cpue_pred)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ yyqq, ncol = 8) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_predictions_gamst_1.jpg'), plot = p1, 
       width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(value = sum(cpue_pred*grid_area*portion_on_ocean))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'gamst_1')
# PredTime$value = scale(PredTime$value)[,1]

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value)) +
  geom_line() +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_gamst_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

