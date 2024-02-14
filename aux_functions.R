
calculate_moran = function(lon, lat, zval, output_type = 'observed') {

  library(ape)
  library(sp)
  these_dist = sp::spDists(x = cbind(lon, lat))
  these_dist_inv = 1/these_dist
  these_dist_inv[which(is.infinite(these_dist_inv))] = 0
  moran_res = ape::Moran.I(zval, these_dist_inv)
  out = moran_res[[output_type]] # observed, p.value
  return(out)
  
}

calculate_clarkevans = function(lon, lat, output_type = 'statistic') {
  
  library(sf)
  points = data.frame(lon = lon, lat = lat)
  MyPoints_temp = points %>% st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  MyPoints_trnsf = sf::st_transform(MyPoints_temp, 32619)
  MyPoints_ppp = as.ppp(MyPoints_trnsf)
  clark_test = clarkevans.test(MyPoints_ppp, correction = "donnelly")
  out = clark_test[[output_type]] # statistic, p.value
  return(out)
  
}

calculate_covarea = function(dat) {
  
  library(sf)
  tmp_df = dat %>% 
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
    st_set_crs("+proj=utm +zone=42N +datum=WGS84 +units=km") 
  cov_area = st_buffer(tmp_df, 1) %>% st_union() %>% st_area() # in km2
  return(cov_area)
  
}

calculate_area_on_land = function(dat) {
  
  library(rgeos)
  library(rnaturalearth)
  world_map = rnaturalearth::ne_countries(scale = 'small', returnclass = c("sf"))
  wm = as(world_map, "Spatial")
  cs = gUnaryUnion(wm, id=as.character(wm$continent))
  cs_sf = st_as_sf(cs)
  inter_grid = st_intersection(cs_sf, dat)
  if(nrow(inter_grid) > 0) area_on_land = sum(as.numeric(st_area(inter_grid)))
  else area_on_land = 0
  
  return(area_on_land)
  
}

hurdle_fn = function(data, i) {
  dat_boot <- data[i, ]
  m1 <- glm(non_zero ~ 1, data = dat_boot,
            family = binomial(link = logit))
  m2 <- glm(y ~ 1, data = subset(dat_boot, non_zero == 1),
            family = Gamma(link = log))
  bin_coef <- plogis(coef(m1)[[1]])
  gamma_coef <- exp(coef(m2)[[1]])
  exp(log(bin_coef) + log(gamma_coef))
}

slope_grid = function(df) {
  n_quarters = length(unique(df$yyqq))
  if(n_quarters > n_all_quarters*0.4) { # cover at least 40% of yyqq
    mod1 = lm(log(catch + 1) ~ time, data = df)
    slopeMod = coef(mod1)[2] 
  } else {
    slopeMod = NA
  }
  return(slopeMod)
}