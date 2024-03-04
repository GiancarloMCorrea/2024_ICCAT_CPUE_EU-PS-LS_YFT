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
theme_set(theme_classic())
source('aux_functions.R')

# -------------------------------------------------------------------------

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images/gam'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/model_outputs'

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model (hurdle approach): 

# Define dataset:
model_df = joinDF %>% mutate(presence = catch > 0)

# Run model: presence/absence
gam_mod_2a = gam(presence ~ s(lon, lat) + year*quarter +
                   den_water2 + capacity + h_sunrise, data = model_df, family = binomial)
summary(gam_mod_2a)

# Plot and check
par(mfrow = c(2,3))
plot(gam_mod_2a)
gam.check(gam_mod_2a)

# Check residuals:
res_model = simulateResiduals(fittedModel = gam_mod_2a)
plot(res_model)

# Check effects:
jtools::effect_plot(gam_mod_2a, pred = den_water2, cat.geom = "line")
jtools::effect_plot(gam_mod_2a, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(gam_mod_2a, pred = capacity, cat.geom = "line")


# Run model
gam_mod_2b = gam(log(catch) ~ s(lon, lat) + year*quarter +
                  den_water2 + capacity + h_sunrise +
                  follow_echo, data = subset(model_df, presence == 1))

# Plot and check
par(mfrow = c(2,3))
plot(gam_mod_2b)
gam.check(gam_mod_2b)

# Check residuals:
res_model = simulateResiduals(fittedModel = gam_mod_2b)
plot(res_model)

# Check effects:
jtools::effect_plot(gam_mod_2b, pred = den_water2, cat.geom = "line")
jtools::effect_plot(gam_mod_2b, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(gam_mod_2b, pred = capacity, cat.geom = "line")
jtools::effect_plot(gam_mod_2b, pred = follow_echo, cat.geom = "line")


# Predict CPUE (grid)
varData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), h_sunrise = mean(model_df$h_sunrise), 
                      den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'),
                      ID = unique(model_df$ID))
predData = full_join(varData, extraDF, by = 'ID')
predData$pres_pred = predict(gam_mod_2a, newdata = predData, type = 'response')
predData$pos_pred = exp(predict(gam_mod_2b, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData$cpue_pred = predData$pres_pred*predData$pos_pred
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
ggsave(filename = file.path(plot_folder, 'grid_predictions_gam_2.jpg'), plot = p1,
       width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(value = sum(cpue_pred*grid_area*portion_on_ocean))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'gam_2')
# PredTime$value = scale(PredTime$value)[,1]

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value)) +
  geom_line() +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_gam_2.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
