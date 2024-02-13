rm(list = ls())

# -------------------------------------------------------------------------
library(INLA)
library(splancs)
library(dplyr)
library(lattice)
library(sf)
library(ggplot2)
library(viridis)
theme_set(theme_classic())

my_border = readRDS('region_extent.rds')

data_folder = 'data'
plot_folder = 'plots'
params = list(species = "SKJ", ORP = "IOTC")

yBreaks = seq(from = -20, to = 20, by = 20)
xBreaks = seq(from = 40, to = 100, by = 20)
limites = read.table(file.path(data_folder, "limites.csv"), header=TRUE,sep=",", na.strings="NA", dec=".", strip.white=TRUE)
limites = subset(limites, ORP == params$ORP)
xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
grid_size = 1

# -------------------------------------------------------------------------
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
worldmap = map_data("world")
data.table::setnames(worldmap, c("X", "Y", "PID", "POS", "region", "subregion"))

# -------------------------------------------------------------------------
# Define dataset:
model_df = joinDF %>% 
  mutate(presence = catch > 0, den_water2 = scale(den_water)[,1], area_swept = 1,
         year2 = as.numeric(as.character(year)), quarter2 = as.numeric(as.character(quarter)),
         time = year2 + (quarter2-1)/4)
st_geometry(model_df) = NULL
model_df = model_df %>% filter(time >= 2020, presence == TRUE)
n_times = length(unique(model_df$time))

# Create mesh
coords = cbind(model_df$lon, model_df$lat)

# Mesh construction
prdomain = inla.nonconvex.hull(points = coords,
                               convex = -0.075, concave = -0.5,
                               resolution=c(30,30))

my_mesh = inla.mesh.2d(loc=coords, boundary = prdomain, 
                    max.edge=c(3, 6), # this change the outer mesh size
                    cutoff = 1, # this change the dense mesh size
                    offset=c(1, 1)) # cutoff
par(mfrow = c(1,1))
plot(my_mesh, main="")
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")
par(mfrow = c(1,1))
plot(my_mesh, main="", xlim = c(50, 51), ylim = c(-1,0))
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")

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
                dplyr::mutate(time = rep(1:n_times, each = n_pred_points))

coop = as.matrix(predPoints[,c('lon', 'lat')])
groupp = as.matrix(predPoints[, 'time'])
Apred = inla.spde.make.A(mesh = my_mesh, loc = coop, group = groupp)

# Stack data and predictions:
stk_e = inla.stack(
  tag = "est",
  data = list(y = model_df$catch),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(model_df))), s = indexs)
)

stk_p = inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, Apred),
  effects = list(data.frame(b0 = rep(1, nrow(predPoints))), s = indexs)
)

stk_full = inla.stack(stk_e, stk_p)

# Run model:
rprior = list(theta = list(prior = "pccor1", param = c(0, 0.9)))
my_formula = y ~ 0 + b0 + f(s, model = spde, group = s.group, 
                            control.group = list(model = 'ar1', hyper = rprior))
inla_mod_1 = inla(my_formula, family='lognormal', control.compute=list(cpo = TRUE, dic=TRUE),
                   data = inla.stack.data(stk_full),
                   control.predictor = list(A=inla.stack.A(stk_full), compute=TRUE), verbose = TRUE)


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
predPoints$pred_mean = exp(inla_mod_1$summary.fitted.values[index, "mean"])
predPoints$pred_sd = inla_mod_1$summary.fitted.values[index, "sd"]
PredGrid = left_join(MyGrid, predPoints, by = c('ID')) %>% na.omit

ggplot() +  
  geom_sf(data = PredGrid, aes(fill = pred_mean, color = pred_mean)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ time, ncol = 3) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none') 

ggplot() +  
  geom_sf(data = PredGrid, aes(fill = pred_sd, color = pred_sd)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  facet_wrap(~ time, ncol = 3) +
  labs(fill = "Predicted mean CPUE") + guides(color = 'none')

# -------------------------------------------------------------------------

# stepsize = 1 # 1 degree
# x_range = diff(range(my_border[,1]))
# y_range = diff(range(my_border[,2]))
# nxy = round(c(x_range, y_range) / stepsize)
# 
# projgrid <- inla.mesh.projector(mesh, xlim = range(my_border[,1]), 
#                                 ylim = range(my_border[,2]), dims = nxy)
# 
# xy.in <- inout(projgrid$lattice$loc, my_border)
# Aprd <- projgrid$proj$A[which(xy.in), ]
# 
# # Create extrapolation grid:
# outGrid = st_centroid(MyGrid) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
# ran_Lon = range(outGrid$lon)
# ran_Lat = range(outGrid$lat)
# 
# lon_centroids = seq(from = ran_Lon[1], to = ran_Lon[2], by = grid_size)
# n_lon = length(lon_centroids)
# lat_centroids = seq(from = ran_Lat[1], to = ran_Lat[2], by = grid_size)
# n_lat = length(lat_centroids)
# nxy = c(n_lon, n_lat)
# 
# projgrid = inla.mesh.projector(mesh, xlim = ran_Lon, ylim = ran_Lat, dims = nxy)
# outGrid = outGrid %>% dplyr::mutate(pred_grid = outGrid$ID %in% extraDF$ID)
# xy.in = outGrid$pred_grid
# Aprd = projgrid$proj$A[which(xy.in), ]
# 
# 
# 
# stk.mesh = inla.stack(data = list(y = NA),
#                       A = list(1, 1),  
#                       effects=list(i = 1:spde$n.spde,
#                                    data.frame(Intercept=rep(1, spde$n.spde))),
#                       tag = 'mesh') 
# 
# stk.b = inla.stack(stk.dat, stk.mesh)
# 
# pred_model = inla(my_formula, family = 'lognormal', 
#                     data = inla.stack.data(stk.b), 
#                     control.predictor = list(A = inla.stack.A(stk.b),
#                                              compute = TRUE, link = 1), 
#                     quantiles = NULL, 
#                     control.compute = list(config = TRUE)) # Needed to sample
# 
# sam_post = inla.posterior.sample(n = 1000, result = pred_model)
# id_prd_mesh = inla.stack.index(stk.b, 'mesh')$data
# pred_nodes = exp(sapply(sam_post, function(x) x$latent[id_prd_mesh]))
# 
# mean_pred = matrix(NA, nxy[1], nxy[2])
# sd_pred = matrix(NA, nxy[1], nxy[2])
# mean_pred[xy.in] = drop(Aprd %*% rowMeans(pred_nodes))
# sd_pred[xy.in] = drop(Aprd %*% apply(pred_nodes, 1, sd))
# 
# # Organize data to plot:
# rownames(mean_pred) = projgrid$x
# colnames(mean_pred) = projgrid$y
# pred_mean_df = as.data.frame(as.table(as.matrix(mean_pred)))
# pred_sd_df = as.data.frame(as.table(as.matrix(sd_pred)))
# pred_mean_df$sdx = pred_sd_df$Freq
# colnames(pred_mean_df) = c('lon', 'lat', 'meanx', 'sdx')
# pred_mean_df$lon = as.numeric(as.character(pred_mean_df$lon))
# pred_mean_df$lat = as.numeric(as.character(pred_mean_df$lat))
# pred_mean_df$ID = 1:nrow(pred_mean_df)
# PredGrid = left_join(MyGrid, pred_mean_df, by = c('ID')) %>% na.omit
# 
# ggplot() +  
#   geom_sf(data = PredGrid, aes(fill = meanx, color = meanx)) + 
#   scale_fill_viridis() + scale_color_viridis() +
#   geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
#   coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
#   xlab(NULL) + ylab(NULL) +
#   theme(legend.position = 'bottom') +
#   scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
#   labs(fill = "Predicted mean CPUE") + guides(color = 'none')
# 
# ggplot() +  
#   geom_sf(data = PredGrid, aes(fill = sdx, color = sdx)) + 
#   scale_fill_viridis() + scale_color_viridis() +
#   geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
#   coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
#   xlab(NULL) + ylab(NULL) +
#   theme(legend.position = 'bottom') +
#   scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
#   labs(fill = "Predicted mean CPUE") + guides(color = 'none')