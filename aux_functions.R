

# -------------------------------------------------------------------------


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
  library(spatstat)
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
  if(nrow(inter_grid) > 0) area_on_land = sum(as.numeric(st_area(inter_grid)))*1e-06 # in km2
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

slope_grid = function(df, prop_obs = 0.4) {
  n_quarters = length(unique(df$time))
  if(n_quarters > n_all_quarters*prop_obs) { # cover at least 40% of yyqq
    mod1 = lm(log(catch + 1) ~ time, data = df)
    slopeMod = coef(mod1)[2] 
  } else {
    slopeMod = NA
  }
  return(slopeMod)
}


# Function to get CI for GLM:
bootstrap_glm = function(models, group_factors, mod_df, pred_df, prop_size = 0.8, n_boot = 100) {
  
  # Prepare prediction matrix:
  n_models = length(models)
  pred_matrix = array(NA, dim = c(n_models, nrow(pred_df), n_boot))
  
  require(dplyr)
  require(stringr)
  for(j in 1:n_boot) {
    sample_df = mod_df %>% group_by_at(group_factors) %>% sample_frac(size = prop_size)
    
    boot_models = list()
    for(i in seq_along(models)) {
      this_family = models[[i]]$family$family
      this_formula = models[[i]]$formula
      this_link = models[[i]]$family$link
      resp_var_name = rownames(attr(models[[i]]$terms, 'factors'))[1]
      resp_var_name = str_remove_all(resp_var_name, paste(c('log', '\\(', '\\)'), collapse = "|"))
      check_log = substr(x = this_formula, start = 1, stop = 3) # way to check if response variable has been log()
      # Run model and make predictions:
      if(check_log == 'log' | this_link == 'log') { 
        boot_models[[i]] = glm(formula = this_formula, data = sample_df %>% filter(!!as.symbol(resp_var_name) > 0), family = this_family) 
        # This way does not work with different links. Check this later
        pred_matrix[i,,j] = exp(predict(boot_models[[i]], newdata = pred_df, type = 'response'))
      } else { 
        boot_models[[i]] = glm(formula = this_formula, data = sample_df, family = this_family)
        pred_matrix[i,,j] = predict(boot_models[[i]], newdata = pred_df, type = 'response')
      }
    } # model loop
    cat(paste0("END iteration: ", j, "\n"))
  } # boot loop
  
  return(pred_matrix)
  
}


center_gravity = function(coord_vec, zval = 1) {
  
  out = weighted.mean(x = coord_vec, w = zval)
  return(out)

}

# Function that will compute a prediction interval for an arbitrary GAM
# Author: David Kaplan
pi_sims_gam = function(model,newdata,n=1e4,alpha=0.05,
                       probs=c(alpha/2,1-alpha/2),return.sims=FALSE) {
  family = sub("[(].*[)]$","",model$family$family)
  invlink = model$family$linkinv
  
  # Model family must be known to this function so as to be able
  # to identify how to draw values 
  # from the distribution family, simulated means and other parameters contained
  # in the model (e.g., the scale parameter). 
  
  # Get simulated coefficients
  beta <- coef(model) # coefficients
  V <- vcov(model) # covariance of coefficients
  Cv <- Matrix::chol(V) # Cholevsky factorisation of covariance
  nus <- rnorm(n * length(beta)) # Unit normally distributed values
  beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = n) # Simulated coefficients
  
  # Get predictions in linear space
  lpmat <- predict(model, newdata = newdata, type = "lpmatrix")
  linpred_sims <- lpmat %*% beta_sims
  
  # Simulated responses after applying inverse link function
  val_sims <- invlink(linpred_sims)
  
  # Simulated predictions drawn from distribution
  y_sims = switch(family,
                  binomial=rbinom(n = prod(dim(val_sims)), size = 1,
                                  prob = val_sims),
                  gaussian=rnorm(n = prod(dim(val_sims)), mean = val_sims, 
                                 sd = sqrt(summary(model)$scale)),
                  poisson=rpois(n = prod(dim(val_sims)), lambda = val_sims),
                  "Negative Binomial"=,
                  nb=rnbinom(n = prod(dim(val_sims)), mu = val_sims,
                             size=model$family$getTheta(TRUE)),
                  "Tweedie"=,
                  tw=mgcv::rTweedie(mu = val_sims,
                                    p = model$family$getTheta(TRUE),
                                    phi = summary(model)$scale),
                  stop("Uknown model family: ",family)) |> 
    matrix(ncol=ncol(val_sims),nrow=nrow(val_sims))
  
  qq = apply(y_sims,1,quantile,probs=probs) |> t() |> as.data.frame()
  
  if (return.sims) {
    res = list(quantiles=qq,sims=y_sims)
  } else {
    res = qq
  }
  
  return(res)
}
