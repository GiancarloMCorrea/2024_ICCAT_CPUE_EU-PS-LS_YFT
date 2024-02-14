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
params = list(species = "SKJ", ORP = "IOTC")

# -------------------------------------------------------------------------
# Read some data inputs created in data_preparation.R:
load(file.path(data_folder, 'joinDF.RData'))
load(file.path(data_folder, 'extraDF.RData'))
load(file.path(data_folder, 'main_df.RData'))
load(file.path(data_folder, 'MyPoints.RData'))
load(file.path(data_folder, 'MyGrid.RData'))

# -------------------------------------------------------------------------
# Some data to plot later:
datos_sp = main_df %>% dplyr::select(lat, lon, year, quarter, catch) %>%
              mutate(yyqq = paste(year, quarter, sep = '-'))# to avoid replacing all objects below

# Map information for plotting:
limites = read.table(file.path(data_folder, "limites.csv"), header=TRUE,sep=",", na.strings="NA", dec=".", strip.white=TRUE)
limites = subset(limites, ORP == params$ORP)
xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
yBreaks = seq(from = -20, to = 20, by = 20)
xBreaks = seq(from = 40, to = 100, by = 20)
worldmap = map_data("world")
colnames(worldmap) = c("X", "Y", "PID", "POS", "region", "subregion")

# -------------------------------------------------------------------------
# Histograms

p1 = ggplot(data = datos_sp, aes(x = catch)) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("Catch (tonnes)")
p2 = ggplot(data = datos_sp %>% filter(catch > 0), aes(x = log(catch))) +
  geom_histogram(aes(y=after_stat(density)), position="identity", alpha=0.5)+
  ylab("Density") + xlab("log(Catch)")
hist_plot = grid.arrange(p1, p2, nrow = 1)
ggsave(filename = file.path(plot_folder, 'hist_catch.jpg'), plot = hist_plot, 
       width = 190, height = 70, units = 'mm', dpi = 500)

# -------------------------------------------------------------------------
# Calculate spatial indices:
ind_df = datos_sp %>% group_by(year, quarter) %>% 
            dplyr::summarise(moran = calculate_moran(lon = lon, lat = lat, zval = catch),
                      gini = Gini(x = catch),
                      #clark = calculate_clarkevans(lon = lon, lat = lat),
                      covarea = calculate_covarea(cur_data()))
plot_df = tidyr::gather(ind_df, 'variable', 'value', 3:ncol(ind_df))
plot_df = plot_df %>% mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
plot_df$variable = factor(plot_df$variable, levels = c('moran', 'gini', 'covarea'),
                               labels = c('Moran index', 'Gini index', 'Covered area'))

ggplot(data = plot_df, aes(time, value)) +
  geom_point() +
  geom_line() +
  ylab('Value') + xlab('Time') +
  facet_wrap(~ variable, scales = 'free_y', ncol = 3)
ggsave(filename = file.path(plot_folder, 'spat_ind.jpg'), width = 190, height = 70, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Time plots:

# Catch per set:
p1 = ggplot(data = datos_sp, aes(yyqq, catch)) +
  geom_boxplot(fill = 'skyblue', width = 0.5, outlier.size = 0.5) +
  ylab('Catch (tonnes) per set') + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8))

# Proportion of zeros and effort:
plot_df = datos_sp %>% group_by(year, quarter) %>% 
              dplyr::summarise(n_obs=n(), prop_zero=length(which(catch == 0))/n())
plot_df = plot_df %>% mutate(year = as.numeric(as.character(year)),
                             quarter = as.numeric(as.character(quarter)),
                             time = year + (quarter-1)/4)

p2 = ggplot(plot_df, aes(time, n_obs)) +
  geom_point() +
  geom_line() +
  ylab('Number of sets') + xlab('Time') 
p3 = ggplot(plot_df, aes(time, prop_zero)) +
  geom_point() +
  geom_line() +
  ylab('Proportion of null sets') + xlab('Time') 

merged_plot = grid.arrange(p1, p2, p3)
ggsave(filename = file.path(plot_folder, 'time_catch.jpg'), plot = merged_plot,
       width = 120, height = 170, units = 'mm', dpi = 500)

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
  theme(legend.position = c(0.9, 0.2)) +
  labs(fill = "Effort") + guides(color = 'none') 
ggsave(filename = file.path(plot_folder, 'grid_eff.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)

# Effort per grid (per year-quarter) -----------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year, quarter) %>% dplyr::summarise(n_obs=n()) %>% mutate(yyqq = paste(year, quarter, sep = '-'))
plot_df = plot_df %>% filter(!is.na(yyqq))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = n_obs, color = n_obs)) + 
  scale_fill_viridis() + scale_color_viridis() +
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = 'bottom') +
  labs(fill = "Effort") + guides(color = 'none') +
  facet_wrap(~ factor(yyqq), ncol = 8)
ggsave(filename = file.path(plot_folder, 'grid_eff_time.jpg'), plot = p1, width = 190, height = 170, units = 'mm', dpi = 500)


# Avg catch per grid (aggregated) -----------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID) %>% dplyr::summarise(avg_catch=mean(catch))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = avg_catch, color = avg_catch)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = c(0.9, 0.2)) +
  labs(fill = "Avg catch (t)") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_catch.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)


# Avg catch per grid (per year-quarter) -------------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year, quarter) %>% dplyr::summarise(avg_catch=mean(catch)) %>% mutate(yyqq = paste(year, quarter, sep = '-'))
plot_df = plot_df %>% filter(!is.na(yyqq))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = avg_catch, color = avg_catch)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = 'bottom') +
  labs(fill = "Avg catch (t)") + guides(color = 'none') +
  facet_wrap(~ factor(yyqq), ncol = 8)
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
  theme(legend.position = c(0.9, 0.2)) +
  labs(fill = "Prop of zeros") + guides(color = 'none')
ggsave(filename = file.path(plot_folder, 'grid_zeroprop.jpg'), plot = p1, width = 190, height = 150, units = 'mm', dpi = 500)
 

# Prop zero per grid (per year-quarter) -------------------------------------------

# Summarise data to plot:
plot_df = joinDF %>% group_by(ID, year, quarter) %>% dplyr::summarise(prop_zero=length(which(catch == 0))/n()) %>% mutate(yyqq = paste(year, quarter, sep = '-'))
plot_df = plot_df %>% filter(!is.na(yyqq))

p1 = ggplot() +  
  geom_sf(data = plot_df, aes(fill = prop_zero, color = prop_zero)) + 
  scale_fill_viridis() + scale_color_viridis() + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  theme(legend.position = 'bottom') +
  labs(fill = "Prop of zeros") + guides(color = 'none') +
  facet_wrap(~ factor(yyqq), ncol = 8)
ggsave(filename = file.path(plot_folder, 'grid_zeroprop_time.jpg'), plot = p1, width = 190, height = 170, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between quarter vs log(catch):

plot_df = datos_sp %>% filter(!is.na(yyqq))

ggplot(data = plot_df, aes(x = quarter, y = log(catch + 1))) +
  geom_boxplot(aes(fill = quarter)) +
  xlab('Quarter') +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  facet_wrap(~ year)
ggsave(filename = file.path(plot_folder, 'catch_season.jpg'), width = 190, height = 160, units = 'mm', dpi = 500)


# -------------------------------------------------------------------------
# Evaluate relationship between variables vs log(catch) per yyqq:

plot_df = main_df %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% filter(!is.na(yyqq))

ggplot(data = plot_df, aes(x = lon, y = log(catch + 1))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_lon.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

ggplot(data = plot_df, aes(x = lat, y = log(catch + 1))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_lat.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

ggplot(data = plot_df, aes(x = h_sunrise, y = log(catch + 1))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_hsunr.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

ggplot(data = plot_df, aes(x = den_water, y = log(catch + 1))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_denwater.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)

ggplot(data = plot_df, aes(x = capacity, y = log(catch + 1))) +
  geom_point() +
  geom_smooth(method = loess, se = FALSE) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_capacity.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)


plot_df = main_df %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% filter(!is.na(yyqq) & !is.na(follow_echo))

ggplot(data = plot_df, aes(x = follow_echo, y = log(catch + 1))) +
  geom_boxplot(aes(fill = follow_echo)) +
  scale_fill_brewer(palette = 'Set2') +
  theme(legend.position = 'none') +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8)) +
  facet_wrap(~ yyqq)
ggsave(filename = file.path(plot_folder, 'catch_follow.jpg'), width = 190, height = 200, units = 'mm', dpi = 500)


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

plot_df = joinDF %>% mutate(yyqq = paste(year, quarter, sep = '-')) %>% 
              filter(!is.na(yyqq)) %>%
              mutate(time = as.numeric(as.character(year)) + (as.numeric(as.character(quarter))-1)/4)
n_all_quarters = length(unique(plot_df$yyqq))

mods = plot_df %>% split(f = plot_df$ID) %>% 
          purrr::map(~ slope_grid(.x))

new_data = data.frame(ID = as.numeric(names(mods)), slope = as.vector(unlist(mods)))
MyGrid2 = left_join(MyGrid, new_data, by = 'ID')
MyGrid2 = MyGrid2 %>% mutate(type_slope = if_else(slope < 0, true = -1, false = 1)) %>% na.omit

ggplot() +  
  geom_sf(data = MyGrid2, aes(fill = factor(type_slope), alpha = abs(slope)), color = 'transparent') + 
  geom_polygon(data = worldmap, aes(X, Y, group=PID), fill = "gray60", color=NA) +
  coord_sf(expand = FALSE, xlim = xLim, ylim = yLim) +
  xlab(NULL) + ylab(NULL) +
  theme(legend.position = c(0.9, 0.2)) +
  guides(fill="none", color = 'none', alpha=guide_legend(title = 'log(catch) trend')) +
  scale_x_continuous(breaks = xBreaks) + scale_y_continuous(breaks = yBreaks) +
  scale_fill_brewer(palette = 'Set1', direction = -1)
ggsave(filename = file.path(plot_folder, 'grid_trend.jpg'), width = 190, height = 150, units = 'mm', dpi = 500)
