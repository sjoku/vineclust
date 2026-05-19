#' internal function
#' @noRd
final_selection <- function(data, total_cluster, final_cvine, final_vinestr, final_trunclevel, mix_probs, p_probs,
                            iteration, init_method, final_mar, final_bicop, trunc_lvl_param, tau_threshold){
  if(is.na(final_cvine)) final_cvine <- 0
  if(is.na(final_trunclevel)) final_trunclevel <- ncol(data) - 1
  
  final_bicop_mapped <- map_family(final_bicop)
  
  data <- cbind(data, apply(p_probs,1,function(x) which(x==max(x))))
  total_obs <- dim(data)[1]
  total_features <- dim(data)[2]-1
  data_cluster <- list()
  
  u_data <- array(0, dim=c(total_obs, total_features, total_cluster))
  vine_models <- list()
  marginal_fams <- matrix(0,total_features, total_cluster)
  marginal_params <- array(0, dim=c(4, total_features, total_cluster))
  
  rvine_densities <- matrix(0, total_obs, total_cluster)
  total_margin_dens <- matrix(0, total_obs, total_cluster)
  margin_densities <- array(0, dim=c(total_obs, total_features, total_cluster))
  lik_points <- matrix(0,dim(data)[1],total_cluster)
  
  total_cop_pars <- 0
  
  p <- progressr::progressor(steps = total_cluster)
  
  for(j in 1:total_cluster){
    data_cluster[[j]] <- data[data[,(total_features+1)] == j,1:total_features, drop=FALSE]
    if(nrow(data_cluster[[j]]) < 5){
      data_cluster[[j]] <- data[,1:total_features, drop=FALSE]
    }
    for(i in 1:total_features){
      min_value <- min(data_cluster[[j]][,i])
      model_margin <- fit_margin(data_cluster[[j]][,i], min_value, final_mar)
      marginal_fams[i,j] <- model_margin$fam
      n_pars <- length(model_margin$par_mar)
      marginal_params[1:n_pars,i,j] <- model_margin$par_mar
    }
    u_data[,,j] <- eval_all_margins_cpp(as.matrix(data[,1:total_features]), marginal_fams[,j], marginal_params[,,j], "cdf")
                                                                       
    trunc_lvl <- NA
    if (!is.na(final_trunclevel)) trunc_lvl <- final_trunclevel
                                                                       
    u_data_cluster <- matrix(u_data[data[,(total_features+1)] == j,,j], ncol=total_features)
    if(nrow(u_data_cluster) < 5){
      u_data_cluster <- u_data[,,j]
    }
    
    if(is.matrix(final_vinestr) || inherits(final_vinestr, "rvine_structure")){
      struct <- rvinecopulib::as_rvine_structure(final_vinestr)
      fit_rvine <- rvinecopulib::vinecop(u_data_cluster, family_set = final_bicop_mapped,
                                         structure = struct, trunc_lvl = trunc_lvl,
                                         keep_data = FALSE, cores = 1)
    }else{
      fit_rvine <- rvinecopulib::vinecop(u_data_cluster, family_set = final_bicop_mapped,
                                         trunc_lvl = trunc_lvl, keep_data = FALSE, cores = 1)
    }
    
    vine_models[[j]] <- fit_rvine
    total_cop_pars <- total_cop_pars + fit_rvine$npars
    
    p(sprintf("Final section %d", j))
  }
  
  data <- data[,1:total_features]
  
  u_data_safe <- u_data
  u_data_safe[u_data_safe < 1e-10] <- 1e-10
  u_data_safe[u_data_safe > 1 - 1e-10] <- 1 - 1e-10
  rvine_densities <- sapply(1:total_cluster, function(j) rvinecopulib::dvinecop(u_data_safe[,,j], vine_models[[j]]))
  
  for(j in 1:total_cluster){
    margin_densities[,,j] <- eval_all_margins_cpp(as.matrix(data), marginal_fams[,j], marginal_params[,,j], "pdf")
  }
  
  margin_densities[margin_densities < 1e-300] <- 1e-300
  rvine_densities[rvine_densities < 1e-300] <- 1e-300
  
  log_lik_points <- matrix(0, nrow=total_obs, ncol=total_cluster)
  for(j in 1:total_cluster) {
    log_m_dens <- rowSums(log(matrix(margin_densities[,,j], nrow=total_obs)))
    log_c_dens <- log(rvine_densities[,j])
    log_lik_points[,j] <- log(mix_probs[j]) + log_m_dens + log_c_dens
  }
  
  max_log_lik <- apply(log_lik_points, 1, max)
  exp_diff <- exp(log_lik_points - max_log_lik)
  sum_exp <- rowSums(exp_diff)
  z_values <- exp_diff / sum_exp
  z_values[is.na(z_values)] <- 1 / total_cluster
  z_values[z_values < 0] <- 0
  z_values[z_values > 1] <- 1
  
  lik_per_obs <- max_log_lik + log(sum_exp)
  loglik <- sum(lik_per_obs)
  
  total_mar_pars <- 0
  for(j in 1:total_cluster){
    for(i in 1:total_features){
      if(marginal_fams[i,j]=='Skew Normal'  || marginal_fams[i,j]=='Student-t'){total_mar_pars <- total_mar_pars + 3}
      else if(marginal_fams[i,j]=='Skew Student-t'){total_mar_pars <- total_mar_pars + 4}
      else{total_mar_pars <- total_mar_pars + 2}
    }
  }
  
  total_mix_pars <- total_cluster-1
  total_pars <- total_mar_pars + total_cop_pars + total_mix_pars
  bic_cop <- (-2)*loglik + log(total_obs)*total_pars
  class <- apply(z_values,1,function(x) which(x==max(x)))
  const <- 0
  if(length(unique(class)) == total_cluster){
    for(i in 1:total_obs){
      cl <- class[i]
      const <- const + log(z_values[i, cl])
    }
    icl <- bic_cop-2*const
  }
  else{icl <- bic_cop}
  
  output <- list("loglik"=loglik, "bic"=bic_cop, "icl"=icl, "init_clustering"=init_method, "iteration"=iteration,
                 "total_pars"=total_pars, "mixture_prob"=mix_probs, "margin"=marginal_fams,
                 "marginal_param"=marginal_params, "vine_models"=vine_models, "z_values"=z_values)
  output
}
