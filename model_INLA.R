rm(list = ls())

# -------------------------------------------------------------------------
library(INLA)
library(splancs)
library(lattice)
my_border = readRDS('borders_extrapolation.rds')

data_folder = 'data'
plot_folder = 'plots'
params = list(species = "SKJ", ORP = "IOTC")

# -------------------------------------------------------------------------
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))


# -------------------------------------------------------------------------
# Define dataset:
model_df = joinDF %>% 
  mutate(presence = catch > 0, den_water2 = scale(den_water)[,1], area_swept = 1,
         year2 = as.numeric(as.character(year)))
st_geometry(model_df) = NULL
model_df = model_df %>% filter(year2 == 2019, presence == TRUE)


# Create mesh
coords = cbind(model_df$lon, model_df$lat)

# Mesh construction
prdomain = inla.nonconvex.hull(points = coords,
                               convex = -0.075, concave = -0.5,
                               resolution=c(30,30))

mesh = inla.mesh.2d(loc=coords, boundary = prdomain, 
                    #min.angle = c(20,20),
                    max.edge=c(1.5, 3), 
                    cutoff = 0.3,
                    offset=c(1, 1)) # cutoff define el lado del triangulo. 
par(mfrow = c(1,1))
plot(mesh, main="")
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")
par(mfrow = c(1,1))
plot(mesh, main="", xlim = c(50, 51), ylim = c(-1,0))
points(coords[,1], coords[,2], cex=0.2, lwd=2, col="red")

# Prediction matrix:
A = inla.spde.make.A(mesh, loc=coords)

# Make SPDE:
spde =  inla.spde2.matern(mesh, alpha=2)  
iset = inla.spde.make.index(name = "spatial.field", spde$n.spde)

# Prepare data for inla
stk.dat = inla.stack(data=list(y = model_df$catch), A=list(A,1), 
                     effects = list(
                       list(i = 1:spde$n.spde), data.frame(Intercept = rep(1, nrow(model_df)))
                     ),
                     tag='dat')

# Run model:
my_formula = y ~ 0 + Intercept + f(i, model = spde) # Sin covariables
inla_mod_1 = inla(my_formula, family='lognormal', control.compute=list(cpo = TRUE, dic=TRUE),
                   data=inla.stack.data(stk.dat),
                   control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE), verbose = TRUE)

# PREDDICIONES DE LOS RESULTADOS: -----------------------------------------
gproj <- inla.mesh.projector(mesh,  dims = c(300, 300))
g.mean <- inla.mesh.project(gproj, inla_mod_1$summary.random$i$mean)
g.sd <- inla.mesh.project(gproj, inla_mod_1$summary.random$i$sd)

grid.arrange(levelplot(g.mean, scales=list(draw=F), xlab='', ylab='', main='mean',col.regions = heat.colors(16)),
             levelplot(g.sd, scal=list(draw=F), xla='', yla='', main='sd' ,col.regions = heat.colors(16)), nrow=1)



reference.coordinates <- coordinates(reference.image)[in.country,]


# MAKE SPATIALPOLYGON OBJECT BORDERS PRED



stepsize = 1 # 1 degree
x_range = diff(range(user_region[, 1]))
y_range = diff(range(user_region[, 2]))
nxy = round(c(x_range, y_range) / stepsize)

projgrid = inla.mesh.projector(mesh, loc = as.matrix(user_region[,c(1,2)]), projection = 'longlat')

# Random effects:
xmean = inla.mesh.project(projgrid, inla_mod_1$summary.random$i$mean)
xsd = inla.mesh.project(projgrid, inla_mod_1$summary.random$i$sd)
xy.in = splancs::inout(projgrid$loc, user_region[,c(1,2)])
xmean[!xy.in] = NA
xsd[!xy.in] = NA

Aprediction = inla.spde.make.A(mesh = mesh, loc = as.matrix(user_region[,c(1,2)]))

Aprd = projgrid$proj$A[which(xy.in), ]
prdcoo = projgrid$lattice$loc[which(xy.in), ]

stk.mesh = inla.stack(data = list(y = NA),
                      A = list(Aprediction, 1),  
                      effects=list(c(list(Intercept=1)
                                     ,iset)),
                      tag='pred') 

stk.mesh = inla.stack(data = list(y = NA),
                      A = list(1, 1),  
                      effects=list(i = 1:spde$n.spde,
                                   data.frame(Intercept=rep(1, spde$n.spde))),
                      tag = 'mesh') 

stk.b <- inla.stack(stk.dat, stk.mesh)

new_model = inla(my_formula, family = 'lognormal', 
                    data = inla.stack.data(stk.b), 
                    control.predictor = list(A = inla.stack.A(stk.b),
                                             compute = TRUE, link = 1), 
                    quantiles = NULL, 
                    control.compute = list(config = TRUE)) # Needed to sample
sampl = inla.posterior.sample(n = 1000, result = new_model)
id.prd.mesh = inla.stack.index(stk.b, 'mesh')$data
pred.nodes = exp(sapply(sampl, function(x) x$latent[id.prd.mesh]))

dim(pred.nodes)
    
sd.prd.s = matrix(NA, nxy[1], nxy[2])
m.prd.s = matrix(NA, nxy[1], nxy[2])

m.prd.s[xy.in] = drop(Aprd %*% rowMeans(pred.nodes))
sd.prd.s[xy.in] = drop(Aprd %*% apply(pred.nodes, 1, sd))


image(m.prd.s)
predict = m.prd.s
Longitud = unique((projgrid$x))
Latitud = unique((projgrid$y))
maxi = round(max(predict*10, na.rm = T))
scale.color = designer.colors(maxi, c('white', 'lightblue','green', 'yellow', 'red','black'))


par(mfrow = c(1,1), mar=c(4,4,0.5,0.5))
image(Longitud, Latitud, predict, col=scale.color, main="",
      zlim = c(0,ceiling(max(predict, na.rm = T))), 
      axes=F, xlim = c(-83, -74), ylim = c(-16, -5), 
      xlab = "longitude", ylab = "latitude")
map("worldHires", add=T, fill=T, col=8)
axis(1, at = seq(-82, -74, by = 2),  labels = 
       seq(-82, -74, by = 2))
axis(2, at = seq(-16, -4, 2) ,  labels = 
       seq(-16, -4, 2), las = 2)
legend.krige(c(-82.5,-82.1),c(-15.5,-9.5), 
             predict, vertical=T, col=scale.color)
box()





predData = st_centroid(extraDF) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
st_geometry(predData) = NULL
xy.in <- inout(projgrid$lattice$loc, PRborder)




## Prediccion del campo aleatorio
gridsize = 60 # Tamano de grilla (en mn2)

j = gridsize
stepsize = j/60 
nxy = round(c(diff(range(user_region[,1])), diff(range(user_region[,2])))/stepsize)

# Projeccion
projgrid = inla.mesh.projector(mesh, xlim = range(user_region[,1]),
                               ylim = range(user_region[,2]), dims = nxy)

xmean = inla.mesh.project(projgrid, inla_mod_1$summary.random$i$mean) #posterior mean
xsd   = inla.mesh.project(projgrid, inla_mod_1$summary.random$i$sd) #posterior sd

#inside the borders 

xy.in = inout(projgrid$lattice$loc,
              cbind(user_region[,1], user_region[,2]))
xmean[!xy.in] = xsd[!xy.in] = NA  #doy valores de N.A fuera del borde

#plot prediction      
jet.colors <- colorRampPalette(c('#00007F','blue','#007FFF','cyan','#7FFF7F','yellow',
                                 '#FF7F00','red','#7F0000') )

png(filename = paste0(outNameSurvey, "_", j, "mn", "-mapPred_Var.png"), width = 500, height = 750, res = 140)
do.call('grid.arrange',
        lapply(list(xmean, xsd),
               levelplot, col.regions=jet.colors(200),
               xlab='', ylab='', scales=list(draw=FALSE)))
dev.off()

Longitud = unique((projgrid$x))
Latitud = unique((projgrid$y))
predict = matrix(data=xmean, ncol=length(Latitud), nrow=length(Longitud))
varianza = matrix(data=xsd, ncol=length(Latitud), nrow=length(Longitud))

write.csv(predict, paste0(outNameSurvey,"_", j, "mn", "-predict_raw.csv"), row.names = F)
write.csv(varianza, paste0(outNameSurvey,"_", j, "mn", "-variance_raw.csv"), row.names = F)

# PREPARAR LOS DATOS PARA HACER FIGURA PUBLICADO
predict[which(predict < 0)] = 0
ztemp = predict * 10
maxi = round(max(ztemp, na.rm = T))

write.csv(predict, paste0(outNameSurvey, "-predict.csv"), row.names = F)
write.csv(Longitud, paste0(outNameSurvey, "-Longitud.csv"), row.names = F)
write.csv(Latitud, paste0(outNameSurvey, "-Latitud.csv"), row.names = F)

z_temp = as.vector(predict)
coordenadas = merge(x = Longitud, y = Latitud, all = T)
outMatrix = data.frame(lon = coordenadas$x, lat = coordenadas$y,
                       var = z_temp)
write.csv(outMatrix, paste0(outNameSurvey,"_", j, "mn", "-predData.csv"), row.names = F)


png(filename = paste0(outNameSurvey,"_", j, "mn", "-mapPred_show.png"), 
    width = 500, height = 600, res = 120)
par(mfrow = c(1,1), mar=c(4,4,0.5,0.5))
image(Longitud, Latitud, predict, col=scale.color, main="",
      zlim = c(0,ceiling(max(predict, na.rm = T))), 
      axes=F, xlim = c(-83, -74), ylim = c(-16, -5), 
      xlab = "longitude", ylab = "latitude")
map("worldHires", add=T, fill=T, col=8)
axis(1, at = seq(-82, -74, by = 2),  labels = 
       seq(-82, -74, by = 2))
axis(2, at = seq(-16, -4, 2) ,  labels = 
       seq(-16, -4, 2), las = 2)
legend.krige(c(-82.5,-82.1),c(-15.5,-9.5), 
             predict, vertical=T, col=scale.color)
box()
dev.off()

