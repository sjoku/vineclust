#' internal function
#' @noRd
final_selection <- function(data, total_cluster, final_cvine, final_vinestr, final_trunclevel, mix_probs, p_probs,
                            iteration, init_method, final_mar, final_bicop){
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
    data_cluster[[j]] <- data[data[,(total_features+1)] == j,1:total_features]
    for(i in 1:total_features){
      min_value <- min(data_cluster[[j]][,i])
      model_margin <- fit_margin(data_cluster[[j]][,i], min_value, final_mar)
      marginal_fams[i,j] <- model_margin$fam
      marginal_params[1,i,j] <- model_margin$par_mar[1]
      marginal_params[2,i,j] <- model_margin$par_mar[2]
      marginal_params[3,i,j] <- model_margin$par_mar[3]
      marginal_params[4,i,j] <- model_margin$par_mar[4]
    }
    u_data[,,j] <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                       marginal_params[,x,j], 'cdf'))
                                                                       
    trunc_lvl <- NA
    if (!is.na(final_trunclevel)) trunc_lvl <- final_trunclevel
                                                                       
    u_data_cluster <- u_data[data[,(total_features+1)] == j,,j]
    
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
  
  u_data_safe <- pmax(pmin(u_data, 1 - 1e-10), 1e-10)
  rvine_densities <- sapply(1:total_cluster, function(j) rvinecopulib::dvinecop(u_data_safe[,,j], vine_models[[j]]))
  
  for(j in 1:total_cluster){
    margin_densities[,,j]<-sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                               marginal_params[,x,j], 'pdf'))
  }
  total_margin_dens <- sapply(1:total_cluster, function(j) apply(margin_densities[,,j], 1, prod))
  lik_points <- sapply(1:total_cluster, function(j) mix_probs[j]*total_margin_dens[,j]*rvine_densities[,j])
  
  lik_per_obs <- apply(lik_points, 1, sum)
  z_values <- t(apply(lik_points, 1, function(row) row / sum(row)))
  z_values[is.na(z_values)] <- 1 / total_cluster
  loglik <- sum(log(lik_per_obs))
  
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
