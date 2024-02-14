rm(list = ls())

# -------------------------------------------------------------------------
# Read libraries:
library(ggplot2)
library(dplyr)  
library(sf)
library(tibble)
source('aux_functions.R')

# -------------------------------------------------------------------------

params = list(species = "SKJ", ORP = "IOTC")
grid_size = 1 # in degrees
data_folder = 'data'

# -------------------------------------------------------------------------

# Set working directory
mydir = getwd()
setwd(mydir)

# Read data:
main_df = readRDS(file = file.path(data_folder, "ecd_fr_es_indian_2010_2021_all_info_v2023_08_14.merge.RDS"))
main_df = main_df  %>%  dplyr::as_tibble()  %>%  dplyr::select(-the_geom) # remove the_geom variable

# Remove duplicated rows:
main_df$dup = duplicated(main_df) # remove duplicates 
main_df = subset(main_df, dup==FALSE)

# Create variables:
main_df = main_df %>% mutate(yyqq = as.factor(paste(annee_de_peche, trimestre, sep="-")),
                             year = as.factor(annee_de_peche),
                             quarter = as.factor(trimestre),
                             pays = as.factor(pays),
                             numbat = as.factor(numbat),
                             lat = latitude_deg + latitude_min/60,
                             lon = longitude_deg + longitude_min/60)
main_df = main_df %>% mutate(lat = if_else(quadrant %in% c(2,3), true = lat*-1, false = lat),
                             lon = if_else(quadrant %in% c(3,4), true = lon*-1, false = lon))

# Select only one stock:
if(params$species=="SKJ") {main_df$catch = main_df$capture_skj_corrigee}
if(params$species=="YFT") {main_df$catch = main_df$capture_yft_corrigee}
if(params$species=="catch") {main_df$catch = main_df$capture_catch_corrigee}

# Filter relevant data:
main_df = main_df %>% filter(pays %in% c(1,4), code_assoc_groupe==1, ocean==2, engin==1, nombre_de_calees>0)
# main_df = main_df %>% filter(!is.na(hrs_since_local_sunrise), !is.na(follow)) # should we remove these observations at this stage?
# Calculate density
main_df = main_df %>% mutate(den = fr_avg_density_rf + es_avg_density_rf,
                             den_water = den/(water_area_1x1_km2*max(total_area_1x1_km2)),
                             follow = as.factor(follow))
# main_df = main_df %>% filter(!is.na(den_water)) # should we remove these observations at this stage?

# follow: follow==1| follow==0
main_df = dplyr::rename(main_df, h_sunrise = hrs_since_local_sunrise)
main_df = dplyr::rename(main_df, capacity = turbobat_cap_m3)
main_df$capacity = as.numeric(main_df$capacity)

# Select relevant variables:
myvars = c("year", "quarter", "h_sunrise", "lat", "lon", "pays", "den_water",
           "numbat", "capacity", "follow", "catch", "buoy_echo")
main_df = main_df[myvars]
main_df = main_df %>% mutate(den_water2 = scale(den_water)[,1], # scale to avoid convergence problems
                             follow_echo = factor(ifelse(follow == '1', buoy_echo, 'no follow'),
                                                  levels = c("no follow","no echo","echo_1freq","echo_2freq"))) %>%
                dplyr::select(-c(follow, buoy_echo)) # remove unused variables
save(main_df, file = file.path(data_folder, 'main_df.RData'))

# -------------------------------------------------------------------------
# Spatial data frames:

limites = read.table(file.path(data_folder, "limites.csv"), header=TRUE,sep=",", na.strings="NA", dec=".", strip.white=TRUE)
limites = subset(limites, ORP == params$ORP)

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
save(MyGrid, file = file.path(data_folder, 'MyGrid.RData'))
