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
require(sf)
theme_set(theme_classic())
source('aux_functions.R')

data_folder = 'data'
plot_folder = 'plots/glm'
model_folder = 'model_outputs'

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'catch_grid_df.RData'))
load(file.path(data_folder, 'clust_area_df.RData'))
load(file.path(data_folder, 'limites.RData'))

# -------------------------------------------------------------------------
# Map information for plotting:
xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
worldmap = map_data("world")
colnames(worldmap) = c("X", "Y", "PID", "POS", "region", "subregion")
yBreaks = seq(from = yLim[1], to = yLim[2], by = 25)
xBreaks = seq(from = xLim[1], to = xLim[2], by = 25)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model 1 (only positives approach): 
# 
# # Define dataset (use only positives or do +1?):
# model_df = joinDF %>% filter(catch > 0)
# # model_df = joinDF %>% mutate(catch = catch + 1)
# # Create cluster column
# n_cluster = nlevels(catch_grid_df$cluster)
# model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')
# st_geometry(model_df) = NULL
# 
# mod1 = glm(log(catch) ~ year + quarter + cluster, data = model_df)
# par(mfrow = c(2,2))
# plot(mod1)
# mod2 = lm(log(catch) ~ year:quarter + year:quarter:cluster, data = model_df)
# par(mfrow = c(2,2))
# plot(mod2)
# mod3 = lm(log(catch) ~ year*quarter + quarter*cluster + year*cluster, data = model_df)
# par(mfrow = c(2,2))
# plot(mod3)
# mod4 = lm(log(catch) ~ year*quarter + quarter*cluster, data = model_df)
# par(mfrow = c(2,2))
# plot(mod4)
# mod5 = lm(log(catch) ~ year*quarter + year*cluster, data = model_df)
# par(mfrow = c(2,2))
# plot(mod5)
# 
# # Compare:
# AIC(mod1); AIC(mod2); AIC(mod3); AIC(mod4); AIC(mod5)
# BIC(mod1); BIC(mod2); BIC(mod3); BIC(mod4); BIC(mod5)
# 
# # GLM:
# glm_mod_1 = glm(log(catch) ~ year*quarter + quarter*cluster + den_water2 + h_sunrise + capacity + follow_echo, 
#                 data = model_df)
# summary(glm_mod_1)
# 
# 
# # Check residuals:
# res_model = simulateResiduals(fittedModel = glm_mod_1)
# plot(res_model)
# par(mfrow = c(1,1))
# qqnorm(qresid(glm_mod_1))
# qqline(qresid(glm_mod_1))
# 
# # Check effects:
# jtools::effect_plot(glm_mod_1, pred = den_water2, cat.geom = "line")
# jtools::effect_plot(glm_mod_1, pred = h_sunrise, cat.geom = "line")
# jtools::effect_plot(glm_mod_1, pred = capacity, cat.geom = "line")
# jtools::effect_plot(glm_mod_1, pred = follow_echo, cat.geom = "line")
# 
# # Check effects interactions
# interactions::cat_plot(model = glm_mod_1, pred = year, modx = quarter)
# 
# # Check inclusion of some variable:
# # alt_model = glm(log(catch) ~ year*quarter + quarter*cluster + year*cluster + den_water2 + h_sunrise + capacity, 
# #                 data = model_df)
# # anova(glm_mod_1, alt_model, test="Chi")
# 
# # TODO: How to impute values when there is no observations for interactions Y*Q*A?
# 
# # Predict CPUE (grid)
# predData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), cluster = as.factor(levels(model_df$cluster)),
#                        h_sunrise = mean(model_df$h_sunrise), 
#                        den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
# predData$cpue_pred = exp(predict(glm_mod_1, newdata = predData, type = 'response')) # , se.fit = T
# PredGrid = left_join(clust_area_df, predData, by = c('cluster'))
# PredGrid = PredGrid %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% na.omit
# 
# ggplot() +  
#   geom_sf(data = PredGrid, aes(fill = cpue_pred, color = cpue_pred)) + 
#   scale_fill_viridis() + scale_color_viridis() +
#   geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
#   coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
#   xlab(NULL) + ylab(NULL) +
#   theme(legend.position = 'bottom') +
#   scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
#   facet_wrap(~ yyqq, ncol = 8) +
#   labs(fill = "Predicted mean CPUE") + guides(color = 'none')
# ggsave(filename = file.path(plot_folder, 'grid_predictions_glm_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
# 
# # Calculate index (weighted sum by area): or average?
# PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(weighted = weighted.mean(cpue_pred, cluster_area),
#                                                               unweighted = mean(cpue_pred))
# PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
#                                model = 'glm_1')
# PredTime = tidyr::gather(PredTime, 'type_cpue', 'value', 3:4)
# 
# # Plot time predictions:
# ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
#   geom_line() +
#   labs(color = 'CPUE type') +
#   theme(legend.position = c(0.85, 0.15)) +
#   ylab('Predicted CPUE (by quarter)') + xlab('Time')
# ggsave(filename = file.path(plot_folder, 'time_predictions_glm_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
# 


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model 2 (hurdle approach): 

# Define dataset:
model_df = joinDF %>% mutate(presence = catch > 0)
# Create cluster column
n_cluster = nlevels(catch_grid_df$cluster)
model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')
st_geometry(model_df) = NULL

# Define formulas:
formula_1 = 'presence ~ year*quarter + year*cluster'
formula_2 = 'log(catch) ~ year*quarter + year*cluster'

# Run model: presence/absence
glm_mod_2a = glm(formula_1, family = binomial, data = model_df)
summary(glm_mod_2a)

# Check residuals:
res_model = simulateResiduals(fittedModel = glm_mod_2a)
plot(res_model)
par(mfrow = c(1,1))
qqnorm(qresid(glm_mod_2a))
qqline(qresid(glm_mod_2a))

# Check effects:
# jtools::effect_plot(glm_mod_2a, pred = den_water2, cat.geom = "line")

# Run model: only positives:
glm_mod_2b = glm(formula_2, data = subset(model_df, presence == 1))
summary(glm_mod_2b)

# Check residuals:
res_model = simulateResiduals(fittedModel = glm_mod_2b)
plot(res_model)
par(mfrow = c(1,1))
qqnorm(qresid(glm_mod_2b))
qqline(qresid(glm_mod_2b))

# Predict CPUE (grid)
predData = expand.grid(year = as.factor(levels(model_df$year)), quarter = as.factor(levels(model_df$quarter)), cluster = as.factor(levels(model_df$cluster)))
                       # h_sunrise = mean(model_df$h_sunrise), 
                       # den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
predData$pres_pred = predict(glm_mod_2a, newdata = predData, type = 'response')
predData$pos_pred = exp(predict(glm_mod_2b, newdata = predData, type = 'response')) # , se.fit = T
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
ggsave(filename = file.path(plot_folder, 'grid_predictions_glm.jpg'), plot = p1, 
       width = 190, height = 150, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------

# Perform bootstrap to calculate ci:
st_geometry(PredGrid) = NULL
merged_models = list(glm_mod_2a, glm_mod_2b)
boot_preds = bootstrap_glm(models = merged_models, mod_df = model_df,
                           group_factors = c('year', 'quarter', 'cluster'), pred_df = PredGrid, n_boot = 5)
boot_preds = apply(X = boot_preds, MARGIN = c(2,3), FUN = prod)

# Calculate index per iteration:
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(label = 'GLM')
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
save_cpue = matrix(NA, ncol = ncol(boot_preds), nrow = nrow(PredTime))
for(k in 1:ncol(boot_preds)) {
  tmp_df = PredGrid %>% mutate(cpue_i = boot_preds[,k])
  tmp_index = tmp_df %>% group_by(year, quarter) %>% summarise(value = sum(cpue_i*cluster_area*portion_on_ocean), .groups = 'drop') # area in km2
  save_cpue[,k] = tmp_index$value
}

# Add columns to df:
PredTime$mean_cpue = apply(save_cpue, 1, mean)
PredTime$sd_cpue = apply(save_cpue, 1, sd)
PredTime$lower_cpue = apply(save_cpue, 1, quantile, 0.025)
PredTime$upper_cpue = apply(save_cpue, 1, quantile, 0.975)
saveRDS(object = PredTime, file = file.path(model_folder, 'GLM_results.rds'))

# Plot time predictions:
scale_factor = 1e-08 # not really important
ggplot(data = PredTime, aes(x = time, y = mean_cpue*scale_factor)) +
         geom_ribbon(aes(ymin = lower_cpue*scale_factor, ymax = upper_cpue*scale_factor), alpha = 0.5) +
         geom_line() +
         ylab('Standardized CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_glm.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
