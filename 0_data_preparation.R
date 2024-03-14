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
# Folder where EU-PS data provided by IEO and IRD is stored:
in_data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine'
# Folder where processed data sets will be stored:
out_data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-LS'

# -------------------------------------------------------------------------

# Set working directory
mydir = getwd()
setwd(mydir)

# Read data:
all_data = readRDS(file = file.path(in_data_folder, "ecd_atlantic_fr_es_2010_2022_followed_density_vms_ts.RDS"))
# Also explore metadata in `all_data`
main_df = all_data$data
st_geometry(main_df) = NULL # remove sf feature since will be added later
main_df$vms_geom = NULL

# Filter some variables:
main_df = main_df %>% filter(c_opera %in% c(0,1)) # code of operation, should be 0 (null) or 1 (positive)
main_df = main_df %>% filter(code_assoc_groupe == 1) # obj (1), free school (2), or ND (3)
main_df = main_df %>% filter(annee_de_peche >= 2010) # from 2010 
main_df = main_df %>% filter(nombre_de_calees > 0) # only observation with positive effort
main_df = main_df %>% filter(pays %in% c(1,4), ocean==1) # only SPA and FRA, and Atlantic Ocean

# Remove duplicated rows:
main_df$dup = duplicated(main_df) # remove duplicates 
main_df = subset(main_df, dup==FALSE)

# Explore variables:
glimpse(main_df)
dim(main_df)

# Create variables:
main_df = main_df %>% mutate(year = as.factor(annee_de_peche),
                             quarter = as.factor(trimestre),
                             time = annee_de_peche + (trimestre-1)/4,
                             country = as.factor(pays),
                             numbat = as.factor(numbat),
                             hold_cap = turbobat_cap_m3,
                             vessel_op = annee_de_peche-turbobat_an_serv, # number of years of vessel activity
                             avg_density = scale(avg_density_1x1_month_water)[,1]) # to avoid convergence issues
                             # lat = latdeg + latmin/60,
                             # lon = londeg + lonmin/60)
# Continuous values for lon and lat:
# main_df = main_df %>% mutate(lat = if_else(quadrant %in% c(2,3), true = lat*-1, false = lat),
#                              lon = if_else(quadrant %in% c(3,4), true = lon*-1, false = lon))

# Create the catch/set column:
main_df = main_df %>% mutate(catch = capture_yft_corrigee/nombre_de_calees) # Make sure you select the right species column

# Select relevant variables:
myvars = c("year", "quarter", "time", "lon", "lat", "catch", "numbat", "country", 
           "hold_cap", "vessel_op",
           "followed", "fads_echo_capable",
           "num_buoys_20nm", "num_owned_20nm", "num_owned_echo_20nm", "num_company_20nm", "num_company_echo_20nm",
           "num_buoys_250km", "num_owned_250km", "num_owned_echo_250km", "num_company_250km", "num_company_echo_250km",
           "avg_density")
# Variables related to time since sunrise have too many NAs, so wont be used here
main_df = main_df[myvars]
dim(main_df)
main_df = main_df %>% na.omit # WARNING: make sure you do not remove too many observations.
dim(main_df)

# Now make follow variable:
main_df$follow = NA
main_df = main_df %>% mutate(follow = ifelse(!followed, yes = 0, follow))
main_df = main_df %>% mutate(follow = ifelse(followed & !fads_echo_capable, yes = 1, follow))
main_df = main_df %>% mutate(follow = ifelse(followed & fads_echo_capable, yes = 2, follow))
main_df$follow = factor(main_df$follow)
main_df = main_df %>% dplyr::select(-c(followed, fads_echo_capable)) # remove unused variables

# -------------------------------------------------------------------------
# Spatial data frames:

# Point and grid (aggregated):
MyPoints = main_df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
min_lon = min(MyPoints$lon) - 0.5
min_lat = min(MyPoints$lat) - 0.5
MyGrid = st_make_grid(MyPoints, cellsize = c(grid_size, grid_size), offset = c(min_lon, min_lat)) %>% 
  st_set_crs(4326) %>% st_sf() %>% dplyr::mutate(ID = 1:n())

# Join both Grid and Points:
joinDF = st_join(MyGrid, left = TRUE, MyPoints) %>% na.omit

# Identify grid (and points inside) with some criteria: 
my_tab = xtabs(~ ID + year, data = joinDF) # find frequency of sets per grid and year
my_freq = apply(my_tab, 1, function(x) sum(x > 0)) # find recurrent grids over the years
include_grid = names(my_freq)[which(as.vector(my_freq) >= 3)] # grids to be excluded because infrequent sample (less than 3 years)
include_grid = as.numeric(include_grid)

# Remove grids:
joinDF = joinDF %>% filter(ID %in% include_grid)
save(joinDF, file = file.path(out_data_folder, 'joinDF.RData'))

# Calculate portion on ocean for grids (important when calculating standardized CPUE):
MyGrid = MyGrid %>% filter(ID %in% include_grid)
IDvec = MyGrid$ID
MyGrid = MyGrid %>% mutate(area_on_land = NA)
for(i in seq_along(IDvec)) {
  tmp = MyGrid %>% filter(ID == IDvec[i])
  MyGrid$area_on_land[i] = calculate_area_on_land(tmp)
}
MyGrid$grid_area = as.numeric(st_area(MyGrid))*1e-06 # in km2
MyGrid = MyGrid %>% mutate(portion_on_ocean = 1 - (area_on_land/grid_area))
save(MyGrid, file = file.path(out_data_folder, 'MyGrid.RData'))

# Limits for plotting later:
MyGrid2 = st_centroid(MyGrid) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
limites = data.frame(xlim1 = floor(min(MyGrid2$lon)/5)*5, ylim1 = floor(min(MyGrid2$lat)/5)*5,
                     xlim2 = ceiling(max(MyGrid2$lon)/5)*5, ylim2 = ceiling(max(MyGrid2$lat)/5)*5)
save(limites, file = file.path(out_data_folder, 'limites.RData'))
