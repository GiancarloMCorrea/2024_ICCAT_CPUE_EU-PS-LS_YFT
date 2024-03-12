rm(list = ls())
library(distances)
library(fpc)
library(sf)
library(ggplot2)
library(dplyr)
theme_set(theme_classic())
source('aux_functions.R')

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-FS'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS_YFT/images'

# -------------------------------------------------------------------------
# Grid size already defined in data_preparation.R
load(file.path(data_folder, 'joinDF.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# Average catch per grid across years:
catch_grid_df = joinDF %>% group_by(ID) %>% dplyr::summarise(avg_catch = mean(catch)) # WARNING: logscale
tmp_grid = st_centroid(catch_grid_df) %>% dplyr::mutate(lon = sf::st_coordinates(.)[,1], lat = sf::st_coordinates(.)[,2])
sel_var_mat = tmp_grid %>% dplyr::select(avg_catch, lon, lat) %>% st_drop_geometry
sel_var_mat = scale(sel_var_mat) # standardize
my_distances = distances(sel_var_mat, weights = c(1,1,1)) # equal weights to each variable
my_clust = fpc::pamk(data = my_distances, krange = 2:6)
catch_grid_df = catch_grid_df %>% mutate(cluster = as.factor(as.vector(my_clust$pamobject$clustering)))

# calculate area of each cluster
clust_area_df = catch_grid_df %>% group_by(cluster) %>% dplyr::summarize(geometry = st_union(geometry))
df_area_land = clust_area_df %>% 
  group_by(cluster) %>% 
  group_map(~ calculate_area_on_land(.x)) %>% 
  unlist() %>%
  tibble::enframe(name = 'cluster', value = 'area_on_land')
df_area_land = df_area_land %>% mutate(cluster = factor(cluster))
clust_area_df = left_join(clust_area_df, df_area_land)
clust_area_df$cluster_area = as.numeric(st_area(clust_area_df))*1e-06 # in km2
clust_area_df = clust_area_df %>% mutate(portion_on_ocean = 1 - (area_on_land/cluster_area))
save(clust_area_df, file = file.path(data_folder, 'clust_area_df.RData'))

# Plot clusters:
p1 = ggplot() +  
  geom_sf(data = catch_grid_df, aes(fill = cluster, color = cluster)) +
  scale_fill_brewer(palette = 'Set1') + scale_color_brewer(palette = 'Set1') +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.9, 0.8)) +
  labs(fill = "Cluster") + guides(color = 'none') 
ggsave(filename = file.path(plot_folder, 'grid_cluster.jpg'), plot = p1, width = 170, height = 160, units = 'mm', dpi = 500)

# Update joinDF with cluster column:
st_geometry(catch_grid_df) = NULL
joinDF = left_join(joinDF, catch_grid_df %>% dplyr::select(ID, cluster), by = 'ID')
save(joinDF, file = file.path(data_folder, 'joinDF.RData'))
