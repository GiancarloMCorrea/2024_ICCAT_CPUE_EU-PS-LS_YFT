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

data_folder = 'data'
plot_folder = 'plots'

# -------------------------------------------------------------------------
# Read some data inputs created in data_preparation.R:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'main_df.RData'))
load(file.path(data_folder, 'MyPoints.RData'))
load(file.path(data_folder, 'MyGrid.RData'))
load(file.path(data_folder, 'limites.RData'))

# -------------------------------------------------------------------------
# Some data to plot later:
datos_sp = main_df

# Map information for plotting:
xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
yBreaks = seq(from = yLim[1], to = yLim[2], by = 25)
xBreaks = seq(from = xLim[1], to = xLim[2], by = 25)
worldmap = map_data("world")
colnames(worldmap) = c("X", "Y", "PID", "POS", "region", "subregion")

# -------------------------------------------------------------------------
# Histograms

p1 = ggplot(data = datos_sp, aes(x = catch)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("Catch (tonnes)")
p2 = ggplot(data = datos_sp %>% filter(catch > 0), aes(x = log(catch))) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("log(catch) (only positive sets)")
hist_plot = grid.arrange(p1, p2, nrow = 1)
ggsave(filename = file.path(plot_folder, 'hist_catch.jpg'), plot = hist_plot, 
       width = 190, height = 70, units = 'mm', dpi = 500)

# -------------------------------------------------------------------------
# Calculate spatial indices:
ind_df = datos_sp %>% group_by(year, quarter) %>% 
            dplyr::summarise(moran = calculate_moran(lon = lon, lat = lat, zval = catch),
                      gini = Gini(x = catch),
                      clark = calculate_clarkevans(lon = lon, lat = lat),
                      covarea = calculate_covarea(cur_data()),
                      c_lon = center_gravity(coord_vec = lon, zval = catch),
                      c_lat = center_gravity(coord_vec = lat, zval = catch))
plot_df = tidyr::gather(ind_df, 'variable', 'value', 3:ncol(ind_df))
plot_df = plot_df %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
plot_df$variable = factor(plot_df$variable, levels = c('clark', 'covarea', 'c_lon', 'c_lat', 'moran', 'gini'),
                               labels = c('Clark-Evans', 'Covered area', 'CG (lon)', 'CG (lat)', 'Moran index', 'Gini index'))

ggplot(data = plot_df, aes(time, value)) +
  geom_point() +
  geom_line() +
  ylab('Value') + xlab('Time') +
  facet_wrap(~ variable, scales = 'free_y', ncol = 3)
ggsave(filename = file.path(plot_folder, 'spat_ind.jpg'), width = 190, height = 140, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Time plots:

# Catch per set:
p1 = ggplot(data = datos_sp, aes(year, catch)) +
  geom_boxplot(fill = 'skyblue', width = 0.5, outlier.size = 0.4) +
  ylab('Catch (tonnes) per set') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))
p2 = ggplot(data = subset(datos_sp, catch > 0), aes(year, log(catch))) +
  geom_boxplot(fill = 'skyblue', width = 0.5, outlier.size = 0.4) +
  ylab('log(catch) per set (only positives)') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))


# Proportion of zeros and effort:
plot_df = datos_sp %>% group_by(year) %>% 
              dplyr::summarise(n_obs=n(), prop_zero=length(which(catch == 0))/n())
plot_df = plot_df %>% mutate(time = as.numeric(as.character(year)))

p3 = ggplot(plot_df, aes(time, n_obs)) +
  geom_point() +
  geom_line() +
  ylab('Number of sets') + xlab('Time') 
p4 = ggplot(plot_df, aes(time, prop_zero)) +
  geom_point() +
  geom_line() +
  ylab('Proportion of null sets') + xlab('Time') 

merged_plot = grid.arrange(p1, p2, p3, p4)
ggsave(filename = file.path(plot_folder, 'time_catch.jpg'), plot = merged_plot,
       width = 190, height = 160, units = 'mm', dpi = 500)

# All grids and points --------------------------------------------------------

ggplot() +  
  geom_sf(data = MyGrid, fill = 'white') + 
  geom_sf(data = MyPoints, color = 'black', alpha = 0.5) + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) 
ggsave(filename = file.path(plot_folder, 'grid_sets.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)

# Grids extrapolation --------------------------------------------------------

ggplot() +  
  geom_sf(data = joinDF, fill = 'white') + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) 
ggsave(filename = file.path(plot_folder, 'grid_extrapo.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)


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
ggsave(filename = file.path(plot_folder, 'grid_eff.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)

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
  facet_wrap(~ factor(year), ncol = 7)
ggsave(filename = file.path(plot_folder, 'grid_eff_time.jpg'), plot = p1, width = 190, height = 170, units = 'mm', dpi = 500)


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
ggsave(filename = file.path(plot_folder, 'grid_catch.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)


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
  facet_wrap(~ factor(year), ncol = 7)
ggsave(filename = file.path(plot_folder, 'grid_catch_time.jpg'), plot = p1, width = 190, height = 170, units = 'mm', dpi = 500)

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
ggsave(filename = file.path(plot_folder, 'grid_zeroprop.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)
 

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
  facet_wrap(~ factor(year), ncol = 7)
ggsave(filename = file.path(plot_folder, 'grid_zeroprop_time.jpg'), plot = p1, width = 190, height = 170, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs effort:
plot_df = joinDF %>% group_by(year, quarter) %>% dplyr::summarise(n_obs=n())

ggplot(data = plot_df, aes(x = quarter, y = n_obs)) +
  geom_col(aes(fill = quarter)) +
  xlab('Quarter') + ylab('Effort') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year, ncol = 7)
ggsave(filename = file.path(plot_folder, 'eff_quarter.jpg'), width = 190, height = 160, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs log(catch):

ggplot(data = subset(datos_sp, catch > 0), aes(x = quarter, y = log(catch))) +
  geom_boxplot(aes(fill = quarter), outlier.size = 0.5) +
  xlab('Quarter') + ylab('log(catch) (only positive sets)') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year, ncol = 7)
ggsave(filename = file.path(plot_folder, 'catch_quarter.jpg'), width = 190, height = 160, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs propZero:
plot_df = joinDF %>% group_by(year, quarter) %>% dplyr::summarise(prop_zero=length(which(catch == 0))/n()) 

ggplot(data = plot_df, aes(x = quarter, y = prop_zero)) +
  geom_col(aes(fill = quarter)) +
  xlab('Quarter') + ylab('Proportion of null sets') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year, ncol = 7)
ggsave(filename = file.path(plot_folder, 'propzero_quarter.jpg'), width = 190, height = 160, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between variables vs log(catch) per yyqq:

plot_df = main_df 

ggplot(data = subset(plot_df, catch > 0), aes(x = lon, y = log(catch))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = file.path(plot_folder, 'catch_lon.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

ggplot(data = subset(plot_df, catch > 0), aes(x = lat, y = log(catch))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ year)
ggsave(filename = file.path(plot_folder, 'catch_lat.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

# ggplot(data = plot_df, aes(x = h_sunrise, y = log(catch + 1))) +
#   geom_point() +
#   geom_smooth(method = loess, se = FALSE) +
#   facet_wrap(~ yyqq)
# ggsave(filename = file.path(plot_folder, 'catch_hsunr.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)
# 
# ggplot(data = plot_df, aes(x = den_water, y = log(catch + 1))) +
#   geom_point() +
#   geom_smooth(method = loess, se = FALSE) +
#   facet_wrap(~ yyqq)
# ggsave(filename = file.path(plot_folder, 'catch_denwater.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)
# 
# ggplot(data = plot_df, aes(x = capacity, y = log(catch + 1))) +
#   geom_point() +
#   geom_smooth(method = loess, se = FALSE) +
#   facet_wrap(~ yyqq)
# ggsave(filename = file.path(plot_folder, 'catch_capacity.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)


# plot_df = main_df %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% filter(!is.na(yyqq) & !is.na(follow_echo))
# 
# ggplot(data = plot_df, aes(x = follow_echo, y = log(catch + 1))) +
#   geom_boxplot(aes(fill = follow_echo)) +
#   scale_fill_brewer(palette = 'Set2') +
#   theme(legend.position = 'none') +
#   xlab(NULL) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
#   facet_wrap(~ yyqq)
# ggsave(filename = file.path(plot_folder, 'catch_follow.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between time vs log(catch) per grid:

# plot_df = joinDF %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% 
#   filter(!is.na(yyqq)) %>%
#   mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
# 
# ggplot(data = plot_df, aes(x = time, y = log(catch + 1))) +
#   geom_point() +
#   geom_smooth(method = loess, se = FALSE) +
#   xlab('Time') +
#   facet_wrap(~ factor(ID))
# ggsave(filename = file.path(plot_folder, 'catch_time_grid.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between time vs log(catch) per grid in a map:

plot_df = joinDF %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
n_all_quarters = length(unique(plot_df$time))

mods = plot_df %>% split(f = plot_df$ID) %>% 
          purrr::map(~ slope_grid(.x, prop_obs = 0.2))

new_data = data.frame(ID = as.numeric(names(mods)), slope = as.vector(unlist(mods)))
MyGrid2 = left_join(MyGrid, new_data, by = 'ID')
MyGrid2 = MyGrid2 %>% mutate(type_slope = if_else(slope < 0, true = -1, false = 1)) %>% na.omit

ggplot() +  
  geom_sf(data = MyGrid2, aes(fill = factor(type_slope), alpha = abs(slope)), color = 'transparent') + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = c(0.9, 0.8)) +
  guides(fill="none", color = 'none', alpha=guide_legend(title = 'log(catch) trend')) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  scale_fill_brewer(palette = 'Set1', direction = -1)
ggsave(filename = file.path(plot_folder, 'grid_trend.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
