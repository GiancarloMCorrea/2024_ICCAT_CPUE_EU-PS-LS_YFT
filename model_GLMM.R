rm(list = ls())
# -------------------------------------------------------------------------
# Load libraries:
require(sjPlot)
require(jtools)
require(DHARMa)
require(ggplot2)
require(viridis)
require(lme4)
require(dplyr)
require(tidyr)
require(sf)
theme_set(theme_classic())

data_folder = 'data'
plot_folder = 'plots/glmm'
params = list(species = "SKJ", ORP = "IOTC")

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'catch_grid_df.RData'))
load(file.path(data_folder, 'clust_area_df.RData'))

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
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model 1 (only positives approach): 

# Define dataset:
model_df = joinDF %>% filter(catch > 0)
# Create cluster column
n_cluster = nlevels(catch_grid_df$cluster)
model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')

# Run model:
glmm_mod_1 = glmer(log(catch) ~ year + quarter + den_water2 + h_sunrise + capacity + follow_echo + (1|year:quarter:cluster),
                   family = gaussian, data = model_df)
# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_mod_1)
plot(res_model)

# Check effects:
jtools::effect_plot(glmm_mod_1, pred = den_water2, cat.geom = "line")
jtools::effect_plot(glmm_mod_1, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(glmm_mod_1, pred = capacity, cat.geom = "line")
jtools::effect_plot(glmm_mod_1, pred = follow_echo, cat.geom = "line")
plot_model(model = glmm_mod_1, type = 're')

# Predict CPUE (grid)
predData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), cluster = as.factor(levels(model_df$cluster)),
                       h_sunrise = mean(model_df$h_sunrise), 
                       den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
predData$cpue_pred = exp(predict(glmm_mod_1, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
PredGrid = left_join(clust_area_df, predData, by = c('cluster'))
PredGrid = PredGrid %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% na.omit

ggplot() +  
  geom_sf(data = PredGrid, aes(fill = cpue_pred, color = cpue_pred)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ yyqq, ncol = 8) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_predictions_glmm_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)


# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(weighted = weighted.mean(cpue_pred, cluster_area),
                                                              unweighted = mean(cpue_pred))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'glmm_1')
PredTime = tidyr::gather(PredTime, 'type_cpue', 'value', 3:4)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
  geom_line() +
  labs(color = 'CPUE type') +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_glmm_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model 2 (hurdle approach): 

# Define dataset:
model_df = joinDF %>%  mutate(presence = catch > 0)
# Create cluster column
n_cluster = nlevels(catch_grid_df$cluster)
model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')

# Run model: presence/absence
glmm_mod_2a = glm(presence ~ year*quarter + cluster + h_sunrise + capacity + den_water2, 
                    family = binomial, data = model_df)
summary(glmm_mod_2a)

# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_mod_2a)
plot(res_model)

# Check effects:
jtools::effect_plot(glmm_mod_2a, pred = den_water2, cat.geom = "line")
jtools::effect_plot(glmm_mod_2a, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(glmm_mod_2a, pred = capacity, cat.geom = "line")


# Run model: only positives:
glmm_mod_2b = glmer(log(catch) ~ year + quarter + den_water2 + h_sunrise + capacity + follow_echo + (1|year:quarter:cluster), 
                 data = subset(model_df, presence == 1))
summary(glmm_mod_2b)

# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_mod_2b)
plot(res_model)

# Predict CPUE (grid)
predData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), cluster = as.factor(levels(model_df$cluster)),
                       h_sunrise = mean(model_df$h_sunrise), 
                       den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
predData$pres_pred = predict(glmm_mod_2a, newdata = predData, type = 'response')
predData$pos_pred = exp(predict(glmm_mod_2b, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData$cpue_pred = predData$pres_pred*predData$pos_pred
PredGrid = left_join(clust_area_df, predData, by = c('cluster'))
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
ggsave(filename = file.path(plot_folder, 'grid_predictions_glmm_2.jpg'), plot = p1, 
       width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(value = sum(cpue_pred*cluster_area))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'glmm_2')
# PredTime$value = scale(PredTime$value)[,1]

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value)) +
  geom_line() +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_glmm_2.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

