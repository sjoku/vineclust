#' internal function
#' @noRd
final_selection <- function(data, total_cluster, mix_probs, z_values, iteration, init_method,
                            marginal_fams, marginal_params, vine_models, loglik) {
  
  total_obs <- dim(data)[1]
  total_features <- dim(data)[2]
  
  total_cop_pars <- 0
  for(j in 1:total_cluster){
    total_cop_pars <- total_cop_pars + vine_models[[j]]$npars
  }
  
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
  class <- max.col(z_values, ties.method = "first")
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
