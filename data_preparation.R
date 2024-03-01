rm(list = ls())

# -------------------------------------------------------------------------
# Read libraries:
library(ggplot2)
library(dplyr)  
library(sf)
library(tibble)
library(readr)
source('aux_functions.R')

# -------------------------------------------------------------------------

grid_size = 1 # in degrees
data_folder = 'data'

# -------------------------------------------------------------------------

# Set working directory
mydir = getwd()
setwd(mydir)

# Read data:
main_df = read_delim(file = file.path(data_folder, "cae_ps_ecd_1991_2022.csv"), delim = ';', show_col_types = FALSE)

# Filter some variables:
main_df = main_df %>% filter(c_opera %in% c(0,1)) # code of operation, should be 0 (null) or 1 (positive)
main_df = main_df %>% filter(codeassocg == 2) # obj (1), free school (2), or ND (3)

# Explore variables:
summary(main_df)
dim(main_df)

# Remove duplicated rows:
main_df$dup = duplicated(main_df) # remove duplicates 
main_df = subset(main_df, dup==FALSE)

# FIX TEMPORARILY LON LAT MIN:
main_df = main_df %>% mutate(latmin = if_else(latmin < 60, true = latmin, false = 59),
                             lonmin = if_else(lonmin < 60, true = lonmin, false = 59))

# Create variables:
main_df = main_df %>% mutate(yyqq = as.factor(paste(year_d_act, quarter, sep="-")),
                             year = as.factor(year_d_act),
                             quarter = as.factor(quarter),
                             lat = latdeg + latmin/60,
                             lon = londeg + lonmin/60)
main_df = main_df %>% mutate(lat = if_else(quadrant %in% c(2,3), true = lat*-1, false = lat),
                             lon = if_else(quadrant %in% c(3,4), true = lon*-1, false = lon))

# Filter relevant data:
main_df = main_df %>% filter(flag %in% c(1,4), ocean==1, engine==1, v_nb_calees>0,
                             lon > -45, lat < 31) # remove observations in the caribe and off Portugal

# Create the catch/set column
main_df = main_df %>% mutate(catch = v_poids_capt_yft/v_nb_calees) # Make sure you select the right column

# Calculate density
# main_df = main_df %>% mutate(den = fr_avg_density_rf + es_avg_density_rf,
#                              den_water = den/(water_area_1x1_km2*max(total_area_1x1_km2)),
#                              follow = as.factor(follow))
# main_df = main_df %>% filter(!is.na(den_water)) # should we remove these observations at this stage?

# follow: follow==1| follow==0
# main_df = dplyr::rename(main_df, h_sunrise = hrs_since_local_sunrise)
# main_df = dplyr::rename(main_df, capacity = turbobat_cap_m3)
# main_df$capacity = as.numeric(main_df$capacity)

# Select relevant variables:
myvars = c("year", "quarter", "lat", "lon", "catch", "vescode")
main_df = main_df[myvars]
# main_df = main_df %>% mutate(den_water2 = scale(den_water)[,1], # scale to avoid convergence problems
#                              follow_echo = factor(ifelse(follow == '1', buoy_echo, 'no follow'),
#                                                   levels = c("no follow","no echo","echo_1freq","echo_2freq"))) %>%
#                 dplyr::select(-c(follow, buoy_echo)) # remove unused variables
save(main_df, file = file.path(data_folder, 'main_df.RData'))

# -------------------------------------------------------------------------
# Spatial data frames:

limites = data.frame(xlim1 = floor(min(main_df$lon)/5)*5, ylim1 = floor(min(main_df$lat)/5)*5,
                     xlim2 = ceiling(max(main_df$lon)/5)*5, ylim2 = ceiling(max(main_df$lat)/5)*5)
save(limites, file = file.path(data_folder, 'limites.RData'))

# Point and grid:
MyPoints = main_df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
MyGrid = st_make_grid(MyPoints, cellsize = c(grid_size, grid_size), offset = c(limites$xlim1, limites$ylim1)) %>% 
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())
save(MyPoints, file = file.path(data_folder, 'MyPoints.RData'))

# Join both Grid and Points:
joinDF = st_join(MyGrid, left = TRUE, MyPoints) %>% na.omit
save(joinDF, file = file.path(data_folder, 'joinDF.RData'))

# # Extrapolation grid:
extraDF = joinDF %>% group_by(ID) %>% dplyr::summarise()
extraDF = st_centroid(extraDF) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
save(extraDF, file = file.path(data_folder, 'extraDF.RData'))

# Calculate portion on ocean for grids (important when calculating CPUE stand):
df_area_land = MyGrid %>% 
  group_by(ID) %>% 
  group_map(~ calculate_area_on_land(.x)) %>% 
  unlist() %>%
  tibble::enframe(name = 'ID', value = 'area_on_land')
MyGrid = left_join(MyGrid, df_area_land)
MyGrid$grid_area = as.numeric(st_area(MyGrid))
MyGrid = MyGrid %>% mutate(portion_on_ocean = 1 - (area_on_land/grid_area))
MyGrid = MyGrid %>% mutate(grid_area = grid_area*1e-06) # in km2
save(MyGrid, file = file.path(data_folder, 'MyGrid.RData'))
