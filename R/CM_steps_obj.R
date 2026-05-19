#' internal function
#' @noRd
CM_step_mixture_probs <- function(z_values) apply(z_values, 2, mean)



#' internal function
#' @noRd
CM_steps <- function(data, vine_model, z_value, marginal_fam, marginal_par, maxit, bicop_mapped, iteration, burn_in_iters, trunc_lvl, tau_threshold, cores = 1, global_min, global_max, global_sd){
  total_features <- dim(data)[2]
  
  if (sum(z_value) < 1e-8) {
      warning("Component weight near zero. Returning previous state to prevent crash.")
      udata <- eval_all_margins_cpp(as.matrix(data), marginal_fam, marginal_par, "cdf")
      udata[udata < 1e-10] <- 1e-10
      udata[udata > 1 - 1e-10] <- 1 - 1e-10
      return(list("marginal_par"=marginal_par, "vine_model"=vine_model, "u_data"=udata))
  }
  
  #CM-step 2
  for(p in 1:total_features){
    if (marginal_fam[p] == 'Normal') {
      w_sum <- sum(z_value)
      mu <- sum(z_value * data[,p]) / w_sum
      var_w <- sum(z_value * (data[,p] - mu)^2) / w_sum
      marginal_par[1:2, p] <- c(mu, max(sqrt(var_w), 1e-6))
      next
    } else if (marginal_fam[p] == 'Lognormal') {
      w_sum <- sum(z_value)
      log_d <- log(pmax(data[,p], 1e-10))
      mu <- sum(z_value * log_d) / w_sum
      var_w <- sum(z_value * (log_d - mu)^2) / w_sum
      marginal_par[1:2, p] <- c(mu, max(sqrt(var_w), 1e-6))
      next
    }
    
    if(marginal_fam[p] %in% c('Normal', 'Lognormal', 'Logistic', 'Cauchy')){
      pars <- marginal_par[1:2,p]
    } else if(marginal_fam[p] %in% c('Gamma', 'Loglogistic')){
      pars <- marginal_par[1:2,p]
    } else if(marginal_fam[p] %in% c('Skew Normal', 'Student-t')){
      pars <- marginal_par[1:3,p]
    } else if(marginal_fam[p] == 'Skew Student-t'){
      pars <- marginal_par[1:4,p]
    }
    
    req_pars <- length(pars)
    lower_b <- rep(-Inf, req_pars)
    upper_b <- rep(Inf, req_pars)
    
    if(marginal_fam[p] %in% c('Normal', 'Lognormal', 'Logistic', 'Cauchy')) {
      lower_b[2] <- max(0.01 * global_sd[p], 1e-6) # Scale/SD must be strictly positive
      upper_b[2] <- max(100 * global_sd[p], 10)
    } else if(marginal_fam[p] %in% c('Gamma', 'Loglogistic')) {
      lower_b[1:2] <- 1e-6                       # Shape and Rate/Scale strictly positive
    } else if(marginal_fam[p] == 'Student-t') {
      lower_b <- c(global_min[p], 0.01 * global_sd[p], 2.0001)
      upper_b <- c(global_max[p], 100 * global_sd[p], 100)
    } else if(marginal_fam[p] == 'Skew Normal') {
      lower_b <- c(global_min[p], 0.01 * global_sd[p], 0.0001)
      upper_b <- c(global_max[p], 100 * global_sd[p], 100)
    } else if(marginal_fam[p] == 'Skew Student-t') {
      lower_b <- c(global_min[p], 0.01 * global_sd[p], 2.0001, 0.0001)
      upper_b <- c(global_max[p], 100 * global_sd[p], 100, 100)
    }
    
    pars[pars < lower_b] <- lower_b[pars < lower_b]
    pars[pars > upper_b] <- upper_b[pars > upper_b]
    
    opt_margins <- optim(par=pars, IFM_margin_obj_cpp, lower = lower_b,
                         upper = upper_b, data_p=data[,p], z_values=z_value, family=marginal_fam[p],
                         method = "L-BFGS-B", control = list(maxit=maxit))
    
    optimized_par <- opt_margins$par
    marginal_par[1:req_pars, p] <- optimized_par
  }
  
  udata <- eval_all_margins_cpp(as.matrix(data), marginal_fam, marginal_par, "cdf")
  udata[udata < 1e-10] <- 1e-10
  udata[udata > 1 - 1e-10] <- 1 - 1e-10
  
  # Re-estimate copula with fixed structure and given weights
  if (iteration <= burn_in_iters) {
     new_vine_model <- tryCatch({
       rvinecopulib::vinecop(udata, family_set = bicop_mapped,
                             weights = z_value, keep_data = FALSE, cores = cores,
                             trunc_lvl = trunc_lvl, threshold = tau_threshold)
     }, error = function(e) {
       if (grepl("negative weight", e$message)) {
         warning("rvinecopulib MST negative weight error detected. Retrying with tree_criterion='rho'")
         tryCatch({
            rvinecopulib::vinecop(udata, family_set = bicop_mapped,
                                  weights = z_value, keep_data = FALSE, cores = cores,
                                  trunc_lvl = trunc_lvl, threshold = tau_threshold,
                                  tree_criterion = "rho")
         }, error = function(e2) {
            warning("rvinecopulib rho tree failed. Returning previous vine_model.")
            vine_model
         })
       } else {
         warning(paste("rvinecopulib crash:", e$message, "| Returning previous vine_model."))
         vine_model
       }
     })
  } else {
     new_vine_model <- tryCatch({
       rvinecopulib::vinecop(udata, family_set = bicop_mapped,
                             structure = rvinecopulib::as_rvine_structure(vine_model), 
                             weights = z_value, keep_data = FALSE, cores = cores,
                             trunc_lvl = trunc_lvl, threshold = tau_threshold)
     }, error = function(e) {
       if (grepl("negative weight", e$message)) {
         warning("rvinecopulib MST negative weight error detected. Retrying with tree_criterion='rho'")
         tryCatch({
            rvinecopulib::vinecop(udata, family_set = bicop_mapped,
                                  structure = rvinecopulib::as_rvine_structure(vine_model), 
                                  weights = z_value, keep_data = FALSE, cores = cores,
                                  trunc_lvl = trunc_lvl, threshold = tau_threshold,
                                  tree_criterion = "rho")
         }, error = function(e2) {
            warning("rvinecopulib rho tree failed. Returning previous vine_model.")
            vine_model
         })
       } else {
         warning(paste("rvinecopulib crash:", e$message, "| Returning previous vine_model."))
         vine_model
       }
     })
  }
                                          
  result <- list("marginal_par"=marginal_par, "vine_model"=new_vine_model, "u_data"=udata)
  result
}
