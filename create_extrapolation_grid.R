# Load libraries:

library(sp) 
library(sf)
library(maps) 
plot_folder = 'plots'

# Get borders:
maps::map(database = "world", xlim = c(limites$xlim1, limites$xlim2), ylim = c(limites$ylim1, limites$ylim2)) 
points(main_df$lon, main_df$lat, col = 'red')
LL = locator()
saveRDS(LL, 'extent_LL.rds')

LL = readRDS('extent_LL.rds')
region_extent = data.frame(long=LL$x, lat=LL$y)
str(region_extent)

region_extent = rbind(region_extent, region_extent[1,])
## https://www.maths.lancs.ac.uk/~rowlings/Teaching/Sheffield2013/cheatsheet.html
poly = Polygon(region_extent)
polys = Polygons(list(poly), ID='all')
sps = SpatialPolygons(list(polys))
## I think the F_AREA could be dropped here
sps = SpatialPolygonsDataFrame(sps, data.frame(Id=factor('all'), F_AREA=1, row.names='all'))
proj4string(sps) = CRS("+proj=longlat +datum=WGS84")
sps = spTransform(sps, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
### Get UTM zone for conversion to UTM projection
## retrieves spatial bounding box from spatial data [,1] is
## longitude
lon = sum(bbox(sps)[1,])/2
## convert decimal degrees to utm zone for average longitude, use
## for new CRS
utmzone <- floor((lon + 180)/6)+1
crs_LL <- CRS('+proj=longlat +ellps=WGS84 +no_defs')
sps@proj4string <- crs_LL

### --------------------------------------------------
### Create the VAST extroplation grid for method 1 and 2
## Convert the final in polygon to UTM
crs_UTM <- CRS(paste0("+proj=utm +zone=",utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
region_polygon <- spTransform(sps, crs_UTM)

### Construct the extroplation grid for VAST using sf package
## Size of grid **in meters** (since working in UTM). Controls
## the resolution of the grid.
cell_size <- 111.12*1000*1 # 5 is number of degrees
## This step is slow at high resolutions
region_grid <- st_make_grid(region_polygon, cellsize = cell_size, what = "centers")
## Convert region_grid to Spatial Points to SpatialPointsDataFrame
region_grid <- as(region_grid, "Spatial")
region_grid_sp <- as(region_grid, "SpatialPointsDataFrame")
## combine shapefile data (region_polygon) with Spatial Points
## (region_grid_spatial) & place in SpatialPointsDataFrame data
## (this provides you with your strata identifier (here called
## Id) in your data frame))
region_grid_sp@data <- over(region_grid, region_polygon)

## Convert back to lon/lat coordinates as that is what VAST uses
region_grid_sp = spTransform(region_grid_sp, crs_LL)
saveRDS(region_grid_sp, file = "borders_extrapolation.rds")
region_grid_LL <- as.data.frame(region_grid_sp)
region_df <- with(region_grid_LL,
                  data.frame(Lon=coords.x1,
                             Lat=coords.x2, Id,
                             Area_km2=( (cell_size/1000)^2),
                             row=1:nrow(region_grid_LL)))
## Filter out the grid that does not overlap (outside extent)
region <- subset(region_df, !is.na(Id))
## This is the final file needed.
str(region)
## > 'data.frame':	106654 obs. of  5 variables:
##  $ Lon     : num  -166 -166 -166 -166 -166 ...
##  $ Lat     : num  53.9 53.9 54 53.9 53.9 ...
##  $ Id      : Factor w/ 1 level "all": 1 1 1 1 1 1 1 1 1 1 ...
##  $ Area_km2: num  4 4 4 4 4 4 4 4 4 4 ...
##  $ row     : int  401 402 975 976 977 978 1549 1550 1551 1552 ...

### Save it to be read in and passed to VAST later.
saveRDS(region, file = "user_region.rds")
### End of creating user extrapolation region object
### --------------------------------------------------

### Quick plots of the process for method 1
png(file.path(plot_folder, 'user_region.png'), width=7, height=7, units='in', res=200)
par(mfrow=c(2,2))
with(region_extent, plot(long, lat, main='Extent in points in LL'))
plot(region_polygon, main='Polygon in UTM', axes=TRUE)
plot(region_grid, col=ifelse(is.na(region_df$Id), 'red', 'black'),
     axes=TRUE, main='Extrapolation area UTM')
with(region, plot(Lon, Lat, main='Extrapolation region in LL', pch='.'))
dev.off()
