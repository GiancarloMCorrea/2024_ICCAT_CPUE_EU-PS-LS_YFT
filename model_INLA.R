rm(list = ls())
# -------------------------------------------------------------------------
library(INLA)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
theme_set(theme_classic())

# my_border = readRDS('region_extent.rds')

# -------------------------------------------------------------------------
data_folder = 'data'
plot_folder = 'plots/inla'
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
model_df = joinDF %>% mutate(presence = catch > 0,
         year2 = as.numeric(as.character(year)), quarter2 = as.numeric(as.character(quarter)),
         time = year2 + (quarter2-1)/4)
model_df = model_df %>% filter(presence == TRUE) 

# Make some changes before running model:
st_geometry(model_df) = NULL
all_years = sort(unique(as.numeric(as.character(unique(model_df$year)))))
all_quarters = sort(unique(as.numeric(as.character(unique(model_df$quarter)))))
time_df = expand.grid(year = all_years, quarter = all_quarters)
time_df = time_df[order(time_df$year),]
time_df = time_df %>% mutate(t_step = 1:n(), yyqq = paste(year, quarter, sep = '-'))
n_times = nrow(time_df)

# Create mesh
coords = cbind(model_df$lon, model_df$lat)
prdomain = inla.nonconvex.hull(points = coords,
                               convex = -0.075, concave = -0.5,
                               resolution=c(30,30))
my_mesh = inla.mesh.2d(loc=coords, boundary = prdomain, 
                    max.edge=c(3, 6), # this change the outer mesh size
                    cutoff = 1, # this change the dense mesh size
                    offset=c(1, 1)) # cutoff

# Save mesh plot:
jpeg(filename = file.path(plot_folder, 'mesh.jpg'), width = 190, height = 120, units = 'mm', res = 500)
par(mfrow = c(1,2))
par(mar = c(0,0,0,0))
plot(my_mesh, main="")
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")
plot(my_mesh, main="", xlim = c(50, 51), ylim = c(-1,0))
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")
dev.off()

# Make A matrix:
spde = inla.spde2.pcmatern(mesh = my_mesh, alpha = 2, constr = T,
                           prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                           prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01
indexs = inla.spde.make.index('s', n.spde = spde$n.spde, n.group = n_times)
group = (model_df$time - min(model_df$time))*4 + 1 # multiply by 4 because we are working with quarters
A = inla.spde.make.A(mesh = my_mesh, loc = coords, group = group)

# Prepare data for prediction:
predPoints = extraDF
n_pred_points = nrow(predPoints)
st_geometry(predPoints) = NULL
predPoints = predPoints %>% slice(rep(1:n(), times = n_times)) %>% 
                dplyr::mutate(t_step = rep(1:n_times, each = n_pred_points))
predPoints = left_join(predPoints, time_df, by = 't_step')

# Spatial and temporal variables:
coop = as.matrix(predPoints[,c('lon', 'lat')])
groupp = as.matrix(predPoints[, 't_step'])
Apred = inla.spde.make.A(mesh = my_mesh, loc = coop, group = groupp)

# Stack data and predictions:
stk_e = inla.stack(
  tag = "est",
  data = list(y = model_df$catch),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(model_df)),
                            den_water2 = model_df$den_water2,
                            capacity = model_df$capacity,
                            h_sunrise = model_df$h_sunrise,
                            follow_echo = model_df$follow_echo), 
                 s = indexs)
)

stk_p = inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Apred),
  effects = list(data.frame(b0 = rep(1, nrow(predPoints)), 
                            den_water2 = mean(model_df$den_water2),
                            capacity = mean(model_df$capacity),
                            h_sunrise = mean(model_df$h_sunrise),
                            follow_echo = factor('no follow')), 
                 s = indexs)
)

stk_full = inla.stack(stk_e, stk_p)

# Run model:
#rprior = list(theta = list(prior = "pccor1", param = c(0, 0.9)))
my_formula = y ~ 0 + b0 + den_water2 + capacity + h_sunrise + follow_echo + 
                  f(s, model = spde, group = s.group, control.group = list(model = 'ar1'))
                    
inla_mod_1 = inla(my_formula, family='lognormal', control.compute=list(cpo = TRUE, dic=TRUE),
                   data = inla.stack.data(stk_full),
                   control.predictor = list(A=inla.stack.A(stk_full), compute=TRUE), verbose = TRUE)

save(inla_mod_1, file = 'inla_mod_1.RData')

# Prepare data for inla
# sdat = inla.stack(data = list(y = model_df$catch), 
#                   A = list(1,A),
#                   effects = list(data.frame(b0 = rep(1, nrow(model_df))), s = indexs),
#                   tag = 'est')

# Prediction matrix:
# A = inla.spde.make.A(mesh, loc=coords)
# Make SPDE:
# spde =  inla.spde2.matern(mesh, alpha=2)  
# iset = inla.spde.make.index(name = "spatial.field", spde$n.spde)
# 
# # Prepare data for inla
# stk.dat = inla.stack(data=list(y = model_df$catch), A=list(A,1), 
#                      effects = list(
#                        list(i = 1:spde$n.spde), data.frame(Intercept = rep(1, nrow(model_df)))
#                      ),
#                      tag='dat')
# 
# # Run model:
# my_formula = y ~ 0 + Intercept + f(i, model = spde) # Sin covariables
# inla_mod_1 = inla(my_formula, family='lognormal', control.compute=list(cpo = TRUE, dic=TRUE),
#                    data=inla.stack.data(stk.dat),
#                    control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE), verbose = TRUE)

# PREDDICIONES DE LOS RESULTADOS: -----------------------------------------

index = inla.stack.index(stack = stk_full, tag = "pred")$data
predPoints$cpue_pred = exp(inla_mod_1$summary.fitted.values[index, "mean"])
predPoints$cpue_sd = inla_mod_1$summary.fitted.values[index, "sd"]
PredGrid = left_join(MyGrid, predPoints, by = c('ID')) %>% na.omit
PredGrid = PredGrid %>% na.omit

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
ggsave(filename = file.path(plot_folder, 'grid_predictions_inla_1_mean.jpg'), plot = p1, 
       width = 190, height = 150, units = 'mm', dpi = 500)

p1 = ggplot() +  
  geom_sf(data = PredGrid, aes(fill = cpue_sd, color = cpue_sd)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ yyqq, ncol = 8) +
  labs(fill = "Predicted sd CPUE") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_predictions_inla_1_sd.jpg'), plot = p1,
       width = 190, height = 150, units = 'mm', dpi = 500)

# Calculate index (weighted sum by area): or average?
PredGrid = PredGrid %>% mutate(year = as.factor(year), quarter = as.factor(quarter))
PredTime = PredGrid %>% group_by(year, quarter) %>% dplyr::summarise(weighted = sum(cpue_pred*portion_on_ocean),
                                                                     unweighted = sum(cpue_pred))
PredTime = PredTime %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4,
                               model = 'inla_1')
PredTime = tidyr::gather(PredTime, 'type_cpue', 'value', 3:4)

# Plot time predictions:
ggplot(data = PredTime, aes(x = time, y = value, color = factor(type_cpue))) +
  geom_line() +
  labs(color = 'CPUE type') +
  theme(legend.position = c(0.85, 0.15)) +
  ylab('Predicted CPUE (by quarter)') + xlab('Time')
ggsave(filename = file.path(plot_folder, 'time_predictions_inla_1.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

