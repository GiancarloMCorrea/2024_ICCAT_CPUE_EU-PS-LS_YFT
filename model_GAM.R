rm(list = ls())
source('data_preparation.R')

# -------------------------------------------------------------------------
# Load libraries:
require(jtools)
require(interactions)
require(DHARMa)
require(statmod)
require(cplm)
require(GlmSimulatoR)
require(glmmTMB)
require(bbmle)
require(lme4)
require(sjPlot)
require(exactextractr)
require(GLMMadaptive)
require(lme4)
require(pscl)
require(mgcv)
library("rnaturalearth")
library("rnaturalearthdata")

# -------------------------------------------------------------------------
# Model 1 (using only positive observations): 

# Define dataset:
model_df = joinDF %>% filter(catch > 0) %>%
  mutate(den_water2 = scale(den_water)[,1], ID2 = factor(ID))

# Run model
gam_mod_1 = gam(log(catch) ~ te(lon, lat, k=13) + year + quarter +
                     s(den_water2, k = 15) + s(capacity,k = 7) + s(h_sunrise, k = 11) +
                     follow_echo, data = model_df)

# Plot and check
par(mfrow = c(2,3))
plot(gam_mod_1)
gam.check(gam_mod_1)

# Check residuals:
res_model = simulateResiduals(fittedModel = gam_mod_1)
plot(res_model)

# Check effects:
jtools::effect_plot(gam_mod_1, pred = den_water2, cat.geom = "line")
jtools::effect_plot(gam_mod_1, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(gam_mod_1, pred = capacity, cat.geom = "line")
jtools::effect_plot(gam_mod_1, pred = follow_echo, cat.geom = "line")

# Predict CPUE (grid)
varData = expand.grid(year = as.factor(2010:2021), quarter = as.factor(1:4), h_sunrise = mean(model_df$h_sunrise), 
                       den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'),
                       ID = unique(model_df$ID))
predData = full_join(varData, extraDF, by = 'ID')
predData$cpue_pred = exp(predict(gam_mod_1, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData = predData %>% dplyr::select(-geometry)
PredGrid = left_join(MyGrid, predData, by = c('ID'))
PredGrid = PredGrid %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% na.omit

ggplot() +  
  geom_sf(data = PredGrid, aes(fill = cpue_pred, color = cpue_pred)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ yyqq, ncol = 8) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, model_folder, 'grid_predictions_gam_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% dplyr::summarise(weighted = weighted.mean(cpue_pred, portion_on_ocean),
                                                              unweighted = mean(cpue_pred))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'gam_1')
PredTime = tidyr::gather(PredTime, 'type_cpue', 'value', 3:4)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
  geom_line() +
  labs(color = 'CPUE type') +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, model_folder, 'time_predictions_gam_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model 2 (hurdle approach): 

# Define dataset:
model_df = joinDF %>% 
  mutate(presence = catch > 0, den_water2 = scale(den_water)[,1])
# Create cluster column
n_cluster = my_clust$nc
model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')


# Run model: presence/absence
gam_mod_2a = gam(presence ~ te(lon, lat, k=13) + year + quarter +
                   s(den_water2, k = 15) + s(capacity,k = 7) + s(h_sunrise, k = 11), data = model_df, family = binomial)
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
gam_mod_2b = gam(log(catch) ~ te(lon, lat, k=13) + year + quarter +
                  s(den_water2, k = 15) + s(capacity,k = 7) + s(h_sunrise, k = 11) +
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
varData = expand.grid(year = as.factor(2010:2021), quarter = as.factor(1:4), h_sunrise = mean(model_df$h_sunrise), 
                      den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'),
                      ID = unique(model_df$ID))
predData = full_join(varData, extraDF, by = 'ID')
predData$pres_pred = predict(gam_mod_2a, newdata = predData, type = 'response')
predData$pos_pred = exp(predict(gam_mod_2b, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData$cpue_pred = predData$pres_pred*predData$pos_pred
predData = predData %>% dplyr::select(-geometry)
PredGrid = left_join(MyGrid, predData, by = c('ID'))
PredGrid = PredGrid %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% na.omit

ggplot() +  
  geom_sf(data = PredGrid, aes(fill = cpue_pred, color = cpue_pred)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ yyqq, ncol = 8) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, model_folder, 'grid_predictions_gam_2.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% dplyr::summarise(weighted = weighted.mean(cpue_pred, portion_on_ocean),
                                                                     unweighted = mean(cpue_pred))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'gam_1')
PredTime = tidyr::gather(PredTime, 'type_cpue', 'value', 3:4)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
  geom_line() +
  labs(color = 'CPUE type') +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, model_folder, 'time_predictions_gam_2.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
