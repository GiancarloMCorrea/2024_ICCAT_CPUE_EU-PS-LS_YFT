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

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images/glmm'
model_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/model_outputs'

# -------------------------------------------------------------------------
# Load some objects created before:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'catch_grid_df.RData'))
load(file.path(data_folder, 'clust_area_df.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Model (hurdle approach): 

# Define dataset:
model_df = joinDF %>%  mutate(presence = catch > 0)
# model_df = model_df %>%  mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4) %>% 
#   filter(time >= 2015 & time < 2022) # just to make things faster, remove later
# Create cluster column
n_cluster = nlevels(catch_grid_df$cluster)
model_df = left_join(model_df, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')
st_geometry(model_df) = NULL

mod1 = glmmTMB::glmmTMB(catch ~ year*quarter + (1|year:cluster), data = model_df, zi = ~ ., family = glmmTMB::lognormal)
summary(mod1)
res_model = simulateResiduals(fittedModel = mod1)
plot(res_model)

xa = arm::sim(glm_mod_2a)
coef_sims = coef(xa)
hist(coef_sims[,3])


# Define formulas:
formula_1 = 'presence ~ year*quarter + (1|year:cluster)'
formula_2 = 'log(catch) ~ year*quarter + (1|year:cluster)'

# Run model: presence/absence
glmm_mod_2a = glmer(formula_1, family = binomial, data = model_df)
summary(glmm_mod_2a)

# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_mod_2a)
plot(res_model)

# Check effects:
jtools::effect_plot(glmm_mod_2a, pred = den_water2, cat.geom = "line")
jtools::effect_plot(glmm_mod_2a, pred = h_sunrise, cat.geom = "line")
jtools::effect_plot(glmm_mod_2a, pred = capacity, cat.geom = "line")


# Run model: only positives:
glmm_mod_2b = glmer(formula_2, data = subset(model_df, presence == 1))
summary(glmm_mod_2b)

# Check residuals:
res_model = simulateResiduals(fittedModel = glmm_mod_2b)
plot(res_model)

# Predict CPUE (grid)
predData = expand.grid(year = as.factor(sort(unique(as.character(model_df$year)))), 
                       quarter = as.factor(sort(unique(as.character(model_df$quarter)))), 
                       cluster = as.factor(sort(unique(as.character(model_df$cluster)))))
                       # h_sunrise = mean(model_df$h_sunrise), 
                       # den_water2 = mean(model_df$den_water2), capacity = mean(model_df$capacity), follow_echo = as.factor('no follow'))
predData$pres_pred = predict(glmm_mod_2a, newdata = predData, type = 'response', allow.new.levels = TRUE)
predData$pos_pred = exp(predict(glmm_mod_2b, newdata = predData, type = 'response', allow.new.levels = TRUE)) # , se.fit = T
predData$cpue_pred = predData$pres_pred*predData$pos_pred


predData$cpue_pred = predict(mod1, newdata = predData, type = 'response', allow.new.levels=TRUE)
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
ggsave(filename = file.path(plot_folder, 'grid_predictions_glmm.jpg'), plot = p1, 
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

