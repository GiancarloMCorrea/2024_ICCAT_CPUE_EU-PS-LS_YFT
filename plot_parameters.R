
# -------------------------------------------------------------------------

# Plot parameters:
# This will apply for all models.
data_folder = 'C:/Use/OneDrive - AZTI/Data/ICCAT/2024/EU_Purse-seine/YFT'
load(file.path(data_folder, 'limites.RData'))

xLim = c(limites$xlim1, limites$xlim2)
yLim = c(limites$ylim1, limites$ylim2)
worldmap = map_data("world")
colnames(worldmap) = c("X", "Y", "PID", "POS", "region", "subregion")
yBreaks = seq(from = yLim[1], to = yLim[2], by = 25)
xBreaks = seq(from = xLim[1], to = xLim[2], by = 25)

