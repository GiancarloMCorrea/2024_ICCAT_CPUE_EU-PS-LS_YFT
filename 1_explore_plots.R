rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(DescTools)
library(gridExtra)
library(viridis)
library(purrr)
source('aux_functions.R')
theme_set(theme_classic())

data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT-LS'
plot_folder = 'C:/Use/OneDrive - AZTI/My_working_papers/ICCAT/2024/CPUE_EU-PS-LS_YFT/images'

# -------------------------------------------------------------------------
# Read some data inputs created in data_preparation.R:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
# Read plot parameters:
source('plot_parameters.R')

# -------------------------------------------------------------------------
# Histograms

p1 = ggplot(data = joinDF, aes(x = catch)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("Catch (tonnes)")
p2 = ggplot(data = joinDF %>% filter(catch > 0), aes(x = log(catch))) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("log(catch) (only positive sets)")
hist_plot = grid.arrange(p1, p2, nrow = 1)
ggsave(filename = paste0('hist_catch', img_type), path = plot_folder, plot = hist_plot, 
       width = 170, height = 70, units = 'mm', dpi = img_res)

# -------------------------------------------------------------------------
# Calculate spatial indices:
ind_df = joinDF %>% group_by(year, quarter) %>% 
            dplyr::summarise(moran = calculate_moran(lon = lon, lat = lat, zval = catch),
                      gini = Gini(x = catch),
                      clark = calculate_clarkevans(lon = lon, lat = lat),
                      covarea = calculate_covarea(cur_data()),
                      c_lon = center_gravity(coord_vec = lon, zval = catch),
                      c_lat = center_gravity(coord_vec = lat, zval = catch))
st_geometry(ind_df) = NULL
plot_df = tidyr::gather(ind_df, 'variable', 'value', 3:ncol(ind_df))
plot_df = plot_df %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
plot_df$variable = factor(plot_df$variable, levels = c('clark', 'covarea', 'c_lon', 'c_lat', 'moran', 'gini'),
                               labels = c('Clark-Evans', 'Covered area', 'CG (lon)', 'CG (lat)', 'Moran index', 'Gini index'))

ggplot(data = plot_df, aes(time, value)) +
  geom_point() +
  geom_line() +
  ylab('Value') + xlab('Time') +
  facet_wrap(~ variable, scales = 'free_y', ncol = 3)
ggsave(filename = paste0('spat_ind', img_type), path = plot_folder, 
       width = 170, height = 120, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Time plots:

# Catch per set:
p1 = ggplot(data = joinDF, aes(year, catch)) +
  geom_boxplot(fill = 'skyblue', width = 0.5, outlier.size = 0.4) +
  ylab('Catch (tonnes) per set') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))
p2 = ggplot(data = subset(joinDF, catch > 0), aes(year, log(catch))) +
  geom_boxplot(fill = 'skyblue', width = 0.5, outlier.size = 0.4) +
  ylab('log(catch) per set (only positives)') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))


# Proportion of zeros and effort:
plot_df = joinDF %>% group_by(year) %>% 
              dplyr::summarise(n_obs=n(), prop_zero=length(which(catch == 0))/n())
plot_df = plot_df %>% mutate(time = as.numeric(as.character(year)))

p3 = ggplot(plot_df, aes(time, n_obs)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = pretty(plot_df$time)) +
  ylab('Number of sets') + xlab('Time') 
p4 = ggplot(plot_df, aes(time, prop_zero)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = pretty(plot_df$time)) +
  ylab('Proportion of null sets') + xlab('Time') 

merged_plot = grid.arrange(p1, p2, p3, p4)
ggsave(filename = paste0('time_catch', img_type), path = plot_folder,  plot = merged_plot,
       width = 170, height = 120, units = 'mm', dpi = img_res)

# All grids and points --------------------------------------------------------

MyPoints = joinDF
st_geometry(MyPoints) = NULL
MyPoints = MyPoints %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

p1 = ggplot() +  
  geom_sf(data = MyPoints, color = 'black', alpha = 0.5) + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) 
ggsave(filename = paste0('all_sets', img_type), path = plot_folder, plot = p1,
       width = 170, height = 160, units = 'mm', dpi = img_res)

# Grids extrapolation --------------------------------------------------------

p1 = ggplot() +  
  geom_sf(data = MyGrid, fill = 'white') + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) 
ggsave(filename = paste0('grid_extrapo', img_type), path = plot_folder, plot = p1,
       width = 170, height = 160, units = 'mm', dpi = img_res)


# Effort per grid (aggregated) -----------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID) %>% dplyr::summarise(n_obs=n())

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = n_obs, color = n_obs)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.9, 0.8)) +
  labs(fill = "Effort") + guides(color = 'none') 
ggsave(filename = paste0('grid_eff', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 160, units = 'mm', dpi = img_res)

# Effort per grid (per year) -----------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year) %>% dplyr::summarise(n_obs=n()) 
plot_df = plot_df %>% filter(!is.na(year))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = n_obs, color = n_obs)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.8, 0.08), legend.direction="horizontal") +
  labs(fill = "Effort") + guides(color = 'none') +
  facet_wrap(~ factor(year))
ggsave(filename = paste0('grid_eff_time', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 170, units = 'mm', dpi = img_res)


# Avg catch per grid (aggregated) -----------------------------------------

# Summarise data to plot: TODO filter only positive
plot_df = joinDF %>% group_by(ID) %>% dplyr::summarise(avg_catch=mean(catch))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = avg_catch, color = avg_catch)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.9, 0.8)) +
  labs(fill = "Avg catch (t)") + guides(color = 'none')
ggsave(filename = paste0('grid_catch', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 160, units = 'mm', dpi = img_res)


# Avg catch per grid (per year) -------------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year) %>% dplyr::summarise(avg_catch=mean(catch))
plot_df = plot_df %>% filter(!is.na(year))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = avg_catch, color = avg_catch)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.8, 0.08), legend.direction="horizontal") +
  labs(fill = "Avg catch (t)") + guides(color = 'none') +
  facet_wrap(~ factor(year))
ggsave(filename = paste0('grid_catch_time', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 170, units = 'mm', dpi = img_res)

# Prop of zeros (aggregated) -----------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID) %>% dplyr::summarise(prop_zero=length(which(catch == 0))/n())

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = prop_zero, color = prop_zero)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.9, 0.8)) +
  labs(fill = "Prop of zeros") + guides(color = 'none')
ggsave(filename = paste0('grid_zeroprop', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 160, units = 'mm', dpi = img_res)
 

# Prop zero per grid (per year) -------------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year) %>% dplyr::summarise(prop_zero=length(which(catch == 0))/n()) 
plot_df = plot_df %>% filter(!is.na(year))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = prop_zero, color = prop_zero)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.8, 0.08), legend.direction="horizontal") +
  labs(fill = "Prop of zeros") + guides(color = 'none') +
  facet_wrap(~ factor(year))
ggsave(filename = paste0('grid_zeroprop_time', img_type), path = plot_folder, plot = p1, 
       width = 170, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs effort:
plot_df = joinDF %>% group_by(year, quarter) %>% dplyr::summarise(n_obs=n())

p1 = ggplot(data = plot_df, aes(x = quarter, y = n_obs)) +
  geom_col(aes(fill = quarter)) +
  xlab('Quarter') + ylab('Effort') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year)
ggsave(filename = paste0('eff_quarter', img_type), path = plot_folder, plot = p1,
       width = 170, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs log(catch):

p1 = ggplot(data = subset(joinDF, catch > 0), aes(x = quarter, y = log(catch))) +
  geom_boxplot(aes(fill = quarter), outlier.size = 0.5) +
  xlab('Quarter') + ylab('log(catch) (only positive sets)') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_quarter', img_type), path = plot_folder, plot = p1,
       width = 170, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs propZero:
plot_df = joinDF %>% group_by(year, quarter) %>% dplyr::summarise(prop_zero=length(which(catch == 0))/n()) 

p1 = ggplot(data = plot_df, aes(x = quarter, y = prop_zero)) +
  geom_col(aes(fill = quarter)) +
  xlab('Quarter') + ylab('Proportion of null sets') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year)
ggsave(filename = paste0('propzero_quarter', img_type), path = plot_folder, plot = p1,
       width = 170, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Evaluate relationship between variables vs log(catch):

ggplot(data = subset(joinDF, catch > 0), aes(x = lon, y = log(catch))) +
  geom_point() +
  xlab('longitude') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_lon', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = lat, y = log(catch))) +
  geom_point() +
  xlab('latitude') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_lat', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = avg_density, y = log(catch))) +
  geom_point() +
  xlab('avg_density') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_avg_density', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = hold_cap, y = log(catch))) +
  geom_point() +
  xlab('vessel_cap') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_hold_cap', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = vessel_op, y = log(catch))) +
  geom_point() +
  xlab('vessel_op') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_vessel_op', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = num_buoys_20nm, y = log(catch))) +
  geom_point() +
  xlab('num_buoys20nm') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_num_buoys20nm', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)


ggplot(data = subset(joinDF, catch > 0), aes(x = num_buoys_250km, y = log(catch))) +
  geom_point() +
  xlab('num_buoys250km') + ylab('log(catch) (only positive sets)') +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_num_buoys250km', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)


ggplot(data = subset(joinDF, catch > 0), aes(x = country, y = log(catch))) +
  geom_boxplot() +
  xlab('country') + ylab('log(catch) (only positive sets)') +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_country', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)

ggplot(data = subset(joinDF, catch > 0), aes(x = follow, y = log(catch))) +
  geom_boxplot() +
  xlab('followed') + ylab('log(catch) (only positive sets)') +
  facet_wrap(~ year)
ggsave(filename = paste0('catch_followed', img_type), path = plot_folder, width = 170, height = 170, units = 'mm', dpi = img_res)


# -------------------------------------------------------------------------
# Evaluate relationship between time vs log(catch) per grid in a map:

plot_df = joinDF
n_all_quarters = length(unique(plot_df$time))

mods = plot_df %>% split(f = plot_df$ID) %>% 
          purrr::map(~ slope_grid(.x, prop_obs = 0.33))

new_data = data.frame(ID = as.numeric(names(mods)), slope = as.vector(unlist(mods)))
MyGrid2 = left_join(MyGrid, new_data, by = 'ID')
MyGrid2 = MyGrid2 %>% mutate(type_slope = if_else(slope < 0, true = -1, false = 1)) %>% na.omit

p1 = ggplot() +  
  geom_sf(data = MyGrid2, aes(fill = factor(type_slope), alpha = abs(slope)), color = 'transparent') + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = c(0.9, 0.8)) +
  guides(fill="none", color = 'none', alpha=guide_legend(title = 'log(catch) trend')) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  scale_fill_brewer(palette = 'Set1', direction = -1)
ggsave(filename = paste0('grid_trend', img_type), path = plot_folder,  plot = p1, 
       width = 170, height = 160, units = 'mm', dpi = img_res)
