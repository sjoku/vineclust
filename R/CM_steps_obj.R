#' internal function
#' @noRd
CM_step_mixture_probs <- function(z_values) apply(z_values, 2, mean)

#' internal function
#' @noRd
CM_step_margin_params <- function(pars, data, marginal_par, p, marginal_families,
                                  z_values, vine_model){
  total_features <- dim(data)[2]
  pars <- matrix(pars, 4, 1)
  if(marginal_families[p]=='Normal' || marginal_families[p]=='Lognormal' || marginal_families[p]=='Logistic'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- max(1e-10, exp(pars[2]))
  }
  else if(marginal_families[p]=='Skew Normal'  || marginal_families[p]=='Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- max(1e-10, pars[2])
    marginal_par[3,p] <- pars[3]
  }
  else if(marginal_families[p]=='Skew Student-t'){
    marginal_par[1,p] <- pars[1]
    marginal_par[2,p] <- max(1e-10, pars[2])
    marginal_par[3,p] <- pars[3]
    marginal_par[4,p] <- pars[4]
  }
  else{
    marginal_par[1,p] <- max(1e-10, exp(pars[1]))
    marginal_par[2,p] <- max(1e-10, exp(pars[2]))
  }
  u_data <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                marginal_par[,x], 'cdf'))
                                                                
  # Ensure u_data is strictly within (0, 1) to avoid infinite density issues with rvinecopulib
  u_data <- pmax(pmin(u_data, 1 - 1e-10), 1e-10)
  
  rvine_densities <- rvinecopulib::dvinecop(u_data, vine_model)
  
  margin_densities <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_families[x],
                                                                         marginal_par[,x], 'pdf'))
  margin_density <- rep(1, dim(margin_densities)[1])
  for(t in 1:total_features){
    margin_density <- margin_density * margin_densities[,t]
  }
  density <- rvine_densities*margin_density
  density[which(density == 0 | is.na(density))] <- 1e-100
  -sum(z_values*log(density))
}

#' internal function
#' @noRd
CM_steps <- function(data, vine_model, z_value, marginal_fam, marginal_par, maxit, bicop_mapped){
  total_features <- dim(data)[2]
  data_min <- apply(data, 2, min)
  data_max <- apply(data, 2, max)
  data_sd <- pmax(apply(data, 2, sd), sqrt(.Machine$double.eps))
  
  #CM-step 2
  for(p in 1:total_features){
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    else if(marginal_fam[p]!='Skew Normal'  && marginal_fam[p]!='Student-t' && marginal_fam[p]!='Skew Student-t'){
      marginal_par[1,p] <- log(marginal_par[1,p])
      marginal_par[2,p] <- log(marginal_par[2,p])
    }
    pars <- marginal_par[,p]
    if(marginal_fam[p]=='Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 2.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_model=vine_model,
                           method = "L-BFGS-B", control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Normal'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 0.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_model=vine_model,
                           method = "L-BFGS-B", control = list(maxit=maxit))
    }
    else if(marginal_fam[p]=='Skew Student-t'){
      opt_margins <- optim(par=pars, CM_step_margin_params, lower = c(data_min[p], 0.01 * data_sd[p], 2.0001, 0.0001),
                           upper = c(data_max[p], 100 * data_sd[p], 100, 100), data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_model=vine_model,
                           method = "L-BFGS-B", control = list(maxit=maxit))
    }
    else{
      opt_margins <- optim(par=pars, CM_step_margin_params, data=data, marginal_par=marginal_par, p=p,
                           marginal_families=marginal_fam, z_values=z_value, vine_model=vine_model,
                           method = "BFGS", control = list(maxit=maxit))

    }
    optimized_par <- opt_margins$par
    if(marginal_fam[p]=='Normal' || marginal_fam[p]=='Lognormal' || marginal_fam[p]=='Logistic'){
      marginal_par[2,p]  <- max(1e-10, exp(optimized_par[2]))
      marginal_par[1,p]  <- optimized_par[1]
    }
    else if(marginal_fam[p]=='Skew Normal'  || marginal_fam[p]=='Student-t'){
      marginal_par[1,p]  <- optimized_par[1]
      marginal_par[2,p]  <- max(1e-10, optimized_par[2])
      marginal_par[3,p]  <- optimized_par[3]
    }
    else if(marginal_fam[p]=='Skew Student-t'){
      marginal_par[1,p]  <- optimized_par[1]
      marginal_par[2,p]  <- max(1e-10, optimized_par[2])
      marginal_par[3,p]  <- optimized_par[3]
      marginal_par[4,p]  <- optimized_par[4]
    }
    else{
      marginal_par[2,p]  <- max(1e-10, exp(optimized_par[2]))
      marginal_par[1,p]  <- max(1e-10, exp(optimized_par[1]))
    }
  }
  
  #CM-step 3
  udata <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fam[x],
                                                               marginal_par[,x], 'cdf'))
  
  # Re-estimate copula with fixed structure and given weights
  new_vine_model <- rvinecopulib::vinecop(udata, family_set = bicop_mapped,
                                          structure = vine_model$structure, 
                                          weights = z_value, keep_data = FALSE, cores = 1)
                                          
  result <- list("marginal_par"=marginal_par, "vine_model"=new_vine_model, "u_data"=udata)
  result
}
