rm(list = ls())
source('data_preparation.R')


# -------------------------------------------------------------------------
library(VAST)


# -------------------------------------------------------------------------
# Define dataset:
model_df = joinDF %>% 
  mutate(presence = catch > 0, den_water2 = scale(den_water)[,1], area_swept = 1,
         year2 = as.numeric(as.character(year)))
st_geometry(model_df) = NULL
model_df = subset(model_df, year2 == 2019)

settings = make_settings(n_x = 400, Region = 'User',
                          purpose = "index2", bias.correct=FALSE,
                          knot_method = 'grid')
settings$FieldConfig[2,] = 0 ## turn off temporal components

user_region = readRDS('user_region.rds')
vast_mod_1 = fit_model(settings = settings,
                       Lat_i = model_df$lat, Lon_i = model_df$lon,
                       t_i = model_df$year2, b_i = model_df$catch,
                       a_i = model_df$area_swept,
                       input_grid = user_region, run_model = FALSE)
plot_results(vast_mod_1, plot_set=3)




plot(vast_mod_1$spatial_list$MeshList$anisotropic_mesh)
plot(vast_mod_1$spatial_list$loc_i[,2], vast_mod_1$spatial_list$loc_i[,1])


vast.output = vast_mod_1

Extrapolation_List = vast.output$extrapolation_list
Spatial_List = vast.output$spatial_list
TmbData = vast.output$data_list

# recreate Data_Geostat
Data_Geostat = data.table::data.table(obs = TmbData$b_i, ts= as.vector(TmbData$t_iz)+1, knot = Spatial_List$knot_i, a_i = TmbData$a_i)
colnames(Data_Geostat) = c("obs","ts","knot","a_i")
Data_Geostat = cbind(Data_Geostat,Convert_EN_to_LL_Fn.ndd(Spatial_List$loc_i[,"E_km"], Spatial_List$loc_i[,"N_km"], crs.en = attr(Spatial_List$loc_i, "projargs"), crs.ll = attr(Spatial_List$loc_i, "origargs")))
Data_Geostat = as.data.frame(Data_Geostat)
color.palette=c("#017cdf","#db00ab","#00b4bb","#ec6900","#b74f66","#738a38")
spatial.agg.df = as.data.frame(Extrapolation_List$Data_Extrap)
spatial.agg.df$knot_i = Spatial_List$PolygonList$NN_Extrap$nn.idx[,1]
ll_to_EN = Convert_LL_to_EastNorth_Fn.ndd( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'],crs.en = attr(Spatial_List$loc_i, "projargs"), crs.ll = attr(Spatial_List$loc_i, "origargs"))


spat.agg.df.dup = spatial.agg.df
spat.agg.df.dup$Cell_ID = 1:nrow(spat.agg.df.dup)
spat.agg.df.dup = spat.agg.df.dup[,c("Lon","Lat","Cell_ID")]
sp::coordinates(spat.agg.df.dup) = ~ Lon + Lat
sp::gridded(spat.agg.df.dup) = TRUE
spat.agg.ras = raster::raster(spat.agg.df.dup)
raster::crs(spat.agg.ras) = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
# create key to map back and forth between extrap.df rows and raster cells	
spat.agg.index.key = spat.agg.ras[1:raster::ncell(spat.agg.ras)]  # this vector gives the row index of Extrapolation_List$Data_Extrap that corresponds with each raster cell
# i.e. spat.agg.index.key[1] = 10542 which means that raster cell 1 corresponds to spat.agg.ras row 10542
# assign knot membership to raster
spat.agg.ras[1:raster::ncell(spat.agg.ras)] = spatial.agg.df$knot_i[spat.agg.index.key]
# get "knot" adjacent to each extrapolation cell 
# ~30 seconds using Microsoft Open R
knot.adj.extrapcell.list = lapply(1:raster::ncell(spat.agg.ras),function(x)unique(spat.agg.ras[raster::adjacent(spat.agg.ras,x,pairs=FALSE)]))
# define knot adjacency matrix
knot.adj.mat = matrix(0,ncol=length(unique(spatial.agg.df$knot_i)),nrow=length(unique(spatial.agg.df$knot_i)))
rownames(knot.adj.mat) = colnames(knot.adj.mat) = 1:length(unique(spatial.agg.df$knot_i))
# fill-in the knot adjacency matrix
for(i in 1:nrow(knot.adj.mat))
{
  knot.adj.mat[i,sort(unique(unlist(knot.adj.extrapcell.list[which(spat.agg.ras[1:raster::ncell(spat.agg.ras)]==i)])))] = 1
  knot.adj.mat[i,i] = 0 # no self loops along the diagonal
}
# use color algorithm to get colors for plotting
# return vector of length (number of knots) where each entry is an integer represnting a color
# implement Welsh-Powell algorithm for graph coloring

knot.color.index = rep(0,nrow(knot.adj.mat))
knot.degree.index = as.numeric(names(sort(rowSums(knot.adj.mat),decreasing=TRUE)))
current.color = 1

'%ni%' = Negate('%in%')
while(length(which(knot.color.index == 0))>0)
{
  colored.index = c()
  remove.index = c()
  for(i in 1:length(knot.degree.index))
  {
    if(all(which(knot.adj.mat[knot.degree.index[i],]==1) %ni% colored.index))
    {
      remove.index = c(remove.index,i)
      colored.index = c(colored.index,knot.degree.index[i])
      knot.color.index[knot.degree.index[i]] = current.color
    }
  }
  knot.degree.index = knot.degree.index[-remove.index]
  current.color = current.color + 1
}

col.vec = color.palette[knot.color.index]


obs.loc = Data_Geostat[,c('Lon','Lat')]
knot.loc = Convert_EN_to_LL_Fn.ndd(Spatial_List$loc_x[,1],Spatial_List$loc_x[,2], crs.en = attr(Spatial_List$loc_i, "projargs"), crs.ll = attr(Spatial_List$loc_i, "origargs"))
points.to.plot = obs.loc
points.to.plot.cols = col.vec[Spatial_List$knot_i]

# main plotting
# define plotting window
plt.xmin = min(c(raster::extent(region.shp.list)@xmin, min(points.to.plot$Lon), min(knot.loc[,"Lon"])))
plt.xmax = max(c(raster::extent(region.shp.list)@xmax, max(points.to.plot$Lon), max(knot.loc[,"Lon"])))
plt.ymin = min(c(raster::extent(region.shp.list)@ymin, min(points.to.plot$Lat), min(knot.loc[,"Lat"])))
plt.ymax = max(c(raster::extent(region.shp.list)@ymax, max(points.to.plot$Lat), max(knot.loc[,"Lat"])))


if(missing(save.dir))
{
  plot(points.to.plot,type="n",xlab="",ylab="",xlim=c(plt.xmin,plt.xmax),ylim=c(plt.ymin,plt.ymax),frame.plot=FALSE,axes=FALSE,asp=1)
  sp::plot(coast.shp,col="gray90",border="gray80",add=TRUE)
  
  # points(knot.loc,pch=16,cex=1.25,col=col.vec)
  # draw links
  for(i in 1:nrow(knot.adj.mat))
  {
    for(j in 1:ncol(knot.adj.mat))
    {
      if(knot.adj.mat[i,j]==1)
      {
        lines(c(knot.loc[i,1],knot.loc[j,1]),c(knot.loc[i,2],knot.loc[j,2]),col="black",lwd=2)
      }
    }
  }
  points(points.to.plot,pch=16,cex=0.5,col=scales::alpha(points.to.plot.cols,0.8))
  sp::plot(region.shp.list,col=NA,border="white",lwd=8,add=TRUE)
  sp::plot(region.shp.list,col=NA,border="black",lwd=4,add=TRUE)
  points(knot.loc,pch=16,cex=2,col="black")
  
} else {
  if (! dir.exists(save.dir))dir.create(save.dir,recursive=TRUE)
  png(filename=paste0(save.dir,"knots.png"),width = 16, height = 9, res=300 ,units = "in",  bg = "white", type = "windows")
  plot(points.to.plot,type="n",xlab="",ylab="",xlim=c(plt.xmin,plt.xmax),ylim=c(plt.ymin,plt.ymax),frame.plot=FALSE,axes=FALSE,asp=1)
  sp::plot(coast.shp,col="gray90",border="gray80",add=TRUE)
  
  # points(knot.loc,pch=16,cex=1.25,col=col.vec)
  # draw links
  for(i in 1:nrow(knot.adj.mat))
  {
    for(j in 1:ncol(knot.adj.mat))
    {
      if(knot.adj.mat[i,j]==1)
      {
        lines(c(knot.loc[i,1],knot.loc[j,1]),c(knot.loc[i,2],knot.loc[j,2]),col="black",lwd=2)
      }
    }
  }
  points(points.to.plot,pch=16,cex=0.5,col=scales::alpha(points.to.plot.cols,0.8))
  sp::plot(region.shp.list,col=NA,border="white",lwd=8,add=TRUE)
  sp::plot(region.shp.list,col=NA,border="black",lwd=4,add=TRUE)
  points(knot.loc,pch=16,cex=2,col="black")
  dev.off()	
