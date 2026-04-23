#' internal function
#' @noRd
initial_clustering <- function(data, total_cluster, is_cvine, init_vinestr, init_trunclevel, init_mar,
                               init_bicop, clustering_method){
  if(is.na(is_cvine)) is_cvine <- 0
  
  init_bicop_mapped <- map_family(init_bicop)
  
  total_features <- dim(data)[2]
  total_obs <- dim(data)[1]
  u_data_to_cluster <- list()
  data_to_cluster <- list()
  u_data <- array(0, dim=c(total_obs, total_features, total_cluster))
  
  vine_models <- list()
  
  marginal_fams <- matrix(0,total_features, total_cluster)
  marginal_params <- array(0, dim=c(4, total_features, total_cluster))
  if(clustering_method == 'gmm'){
    gmm_fit <- mclust::Mclust(data, total_cluster, verbose=FALSE)
  }
  if(clustering_method == 'kmeans'){
    kmeans_fit <- kmeans(scale(data), total_cluster)
  }
  if(clustering_method == 'hcVVV'){
    hcVVV_fit <- mclust::hcVVV(data=scale(data), alpha = 1)
    hcVVV_cl <- mclust::hclass(hcVVV_fit, total_cluster)
  }
  mix_probs <- vector()
  
  # if progressr is used higher up, this will signal progress
  p <- progressr::progressor(steps = total_cluster)
  
  for(j in 1:total_cluster){
    if(clustering_method == 'gmm'){
      data_to_cluster[[j]] <- data[gmm_fit$classification == j,]
    }
    if(clustering_method == 'kmeans'){
      data_to_cluster[[j]] <- data[kmeans_fit$cluster == j,]
    }
    if(clustering_method == 'hcVVV'){
      data_to_cluster[[j]] <- data[hcVVV_cl == j,]
    }
    for(i in 1:total_features){
      min_value <- min(data_to_cluster[[j]][,i])
      model_margin <- fit_margin(data_to_cluster[[j]][,i], min_value, init_mar)
      marginal_fams[i,j] <- model_margin$fam
      marginal_params[1,i,j] <- model_margin$par_mar[1]
      marginal_params[2,i,j] <- model_margin$par_mar[2]
      if(model_margin$fam=='Skew Normal' || model_margin$fam=='Student-t'){
        marginal_params[3,i,j] <- model_margin$par_mar[3]
        }
      if(model_margin$fam=='Skew Student-t'){
        marginal_params[3,i,j] <- model_margin$par_mar[3]
        marginal_params[4,i,j] <- model_margin$par_mar[4]
      }
    }
    u_data[,,j] <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                       marginal_params[,x,j], 'cdf'))
    
    trunc_lvl <- NA
    if (!is.na(init_trunclevel)) trunc_lvl <- init_trunclevel
    
    cluster_u_data <- NULL
    if(clustering_method == 'gmm'){
      cluster_u_data <- u_data[gmm_fit$classification == j,,j]
    } else if(clustering_method == 'kmeans'){
      cluster_u_data <- u_data[kmeans_fit$cluster == j,,j]
    } else if(clustering_method == 'hcVVV'){
      cluster_u_data <- u_data[hcVVV_cl == j,,j]
    }
    
    if(is.matrix(init_vinestr) || inherits(init_vinestr, "rvine_structure")){
      struct <- rvinecopulib::as_rvine_structure(init_vinestr)
      fit_rvm <- rvinecopulib::vinecop(cluster_u_data, family_set = init_bicop_mapped, 
                                       structure = struct, trunc_lvl = trunc_lvl, 
                                       keep_data = FALSE, cores = 1)
    } else {
      fit_rvm <- rvinecopulib::vinecop(cluster_u_data, family_set = init_bicop_mapped, 
                                       trunc_lvl = trunc_lvl, keep_data = FALSE, cores = 1)
    }
    
    vine_models[[j]] <- fit_rvm
    mix_probs[j] <- length(data_to_cluster[[j]][,1])/total_obs
    p(sprintf("Initial clustering (%s) component %d", clustering_method, j))
  }
  
  result <- list("u_data"=u_data, "mix_probs"=mix_probs, "marginal_fams" = marginal_fams,
                 "marginal_params"=marginal_params, "vine_models" = vine_models)
  result
}
