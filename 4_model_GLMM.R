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
source('aux_functions.R')

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-FS'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images/glmm'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/model_outputs/glmm'
# create folders in case missing:
dir.create(plot_folder, showWarnings = FALSE)
dir.create(model_folder, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'clust_area_df.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model (hurdle approach): 

# Define dataset:
model_df = joinDF %>%  mutate(presence = catch > 0)
st_geometry(model_df) = NULL

# Build model:
glmm_model = glmmTMB::glmmTMB(formula = catch ~ year:quarter + (1|year:quarter:cluster), data = model_df, 
                              zi = ~ ., family = glmmTMB::lognormal)
summary(glmm_model)
save(glmm_model, file = file.path(model_folder, 'my_model.RData'))

# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_model)
jpeg(filename = file.path(plot_folder, 'residuals_check.jpg'), width = 190, height = 100, units = 'mm', res = 500)
par(mar = c(4, 4, 1, 0.5))
plot(res_model, title = NULL)
dev.off()

# Predict CPUE (per cluster)

# xa = arm::sim(glm_mod_2a)
# coef_sims = coef(xa)
# hist(coef_sims[,3])

predData = expand.grid(year = as.factor(sort(unique(as.character(model_df$year)))), 
                       quarter = as.factor(sort(unique(as.character(model_df$quarter)))), 
                       cluster = as.factor(sort(unique(as.character(model_df$cluster)))))
                       # h_sunrise = mean(model_df$h_sunrise), 
                       # den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
predData$cpue_pred = predict(glmm_model, newdata = predData, type = 'response', allow.new.levels = TRUE)
PredGrid = left_join(clust_area_df, predData, by = c('cluster'))
plot_data = PredGrid %>% group_by(year, cluster) %>% summarise(cpue_pred = mean(cpue_pred)) # aggregate by year and cluster, just to display in a map

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
ggsave(filename = file.path(plot_folder, 'grid_predictions_glmm.jpg'), plot = p1, 
       width = 190, height = 180, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area)
PredTime = PredGrid %>% group_by(year, quarter) %>% summarise(value = sum(cpue_pred*cluster_area*portion_on_ocean))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4, type = 'GLMM')
st_geometry(PredTime) = NULL
write.csv(PredTime, file = file.path(model_folder, 'CPUE_ts.csv'), row.names = FALSE)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value*1e-08)) +
  geom_line() +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_glmm_2.jpg'), width = 95, height = 70, units = 'mm', dpi = 500)

