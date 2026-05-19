#' Model-Based Clustering with Vine Copulas
#'
#' It fits vine copula based mixture model distributions to the continuous data
#' for a given number of components as described in Sahin and Czado (2021)
#' and use its results for clustering.
#'
#' @param data 	A matrix or data frame of observations. Categorical/discrete variables not (yet) allowed.
#' If a matrix or data frame, rows correspond to observations (i) and columns correspond to variables (p).
#' @param total_comp An integer specifying the numbers of mixture components (clusters)
#' @param is_cvine An integer specifying if the type of components' vine tree structure is C-vine
#' before/after the ECM phase of clustering.
#' * 0 = R-vine (default)
#' * 1 = C-vine
#' @param vinestr A matrix specifying vine tree structures before/after the ECM phase of clustering.
#' The default is automatic selection.
#' \link[rvinecopulib]{rvine_structure} checks for a valid R-vine structure.
#' @param trunclevel An integer showing the level of truncation for vine tree structures before the ECM phase of clustering.
#' The default is 1.
#' @param mar A vector of character strings indicating the parametric univariate marginal distributions
#' to be fitted before/after the ECM phase of clustering.
#' The default is c('cauchy','gamma','llogis','lnorm','logis','norm','snorm','std', 'sstd').
#' Other distributions not (yet) allowed.
#' @param bicop A vector of integers or strings denoting the parametric bivariate copula families to be fitted
#' before/after the ECM phase of clustering.
#' The default is c(1,2,3,4,5,6,7,8,10,13,14,16,17,18,20,23,24,26,27,28,30,33,34,36,37,38,40).
#' \link[rvinecopulib]{bicop_dist} describes the available families with their specifications.
#' @param methods A vector of character strings indicating initial clustering method(s) to have a partition
#' for model selection before the ECM phase of clustering. Current options:
#' * 'kmeans' (default)
#' * c('kmeans', 'gmm', 'hcVVV')
#' @param threshold A numeric, stopping the ECM phase of clustering. The default is 1e-4.
#' @param maxit An integer, specifying the maximum number of iterations in the CM-step 2 optimization. The default is 10.
#' @param cores An integer, showing the number of cores to use for parallel computing.
#' @param verbose A boolean indicating whether to log detailed debugging steps. Defaults to `FALSE`.
#' @param burn_in_iters An integer specifying the number of initial iterations to perform full vine tree structure estimation before freezing it. Defaults to 5.
#' @param trunc_lvl An integer showing the level of truncation for vine tree structures during the CM-steps. Defaults to 2.
#' @param tau_threshold A numeric threshold for Kendall's tau below which pair-copulas are set to independence. Defaults to 0.1.
#' @param batch_size An integer specifying the random sub-sample size for stochastic mini-batch EM. If NULL, uses the full dataset. Defaults to NULL.
#' @param ema_alpha A numeric between 0 and 1 specifying the smoothing factor for Exponential Moving Average updates across batches. Defaults to 0.8.
#'
#' @return An object of class vcmm result. It contains the elements
#' \describe{
#' \item{cluster}{the vector with the classification of observations}
#' \item{output}{a list containing the fitted VCMM.}
#' } Use `print.vcmm_res()` to obtain log-likelihood, BIC, ICL, number of estimated parameters, initial clustering method used
#'  and total number of ECM iterations for the fitted VCMM. `summary.vcmm_res()` shows the fitted vine tree structures and
#'  univariate marginal distributions, bivariate copula families with the estimated parameters, as well as
#'  mixture proportions of each component.
#'
#' @references
#' Sahin and Czado (2021), Vine copula mixture models and clustering for non-Gaussian data, Econometrics and Statistics.
#' doi: 10.1016/j.ecosta.2021.08.011
#'
#' @seealso [dvcmm()], [rvcmm()]
#'
#' @examples
#' \dontrun{
#' # Example: fit parametric 4 dimensional vine copula based mixture model with 2 components
#' # data from UCI Machine Learning Repository 
#' data_wisc <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = FALSE)
#' 
#' # Fit the model
#' fit <- vcmm(data=data_wisc[,c(15,27,29,30)], total_comp=2, verbose=FALSE)
#' 
#' # Display model statistics
#' print(fit)
#' 
#' # Extract vine tree structure information
#' summary(fit)
#' 
#' # Evaluate the density of the fitted model at a given point
#' RVMs_fitted <- fit$output$vine_models
#' dens <- dvcmm(c(2.747, 0.1467, 0.13, 0.05334), fit$output$margin, 
#'               fit$output$marginal_param, RVMs_fitted, fit$output$mixture_prob)
#' }
#'
#' @export
#'
#' @import rvinecopulib
#' @import progressr
#' @import mclust
#' @import univariateML
#' @importFrom fGarch psnorm dsnorm pstd dstd psstd dsstd
#' @importFrom stats dgamma dlnorm dlogis dnorm dcauchy kmeans optim pgamma plnorm plogis pnorm pcauchy sd

vcmm <- function(data, total_comp, is_cvine=NA, vinestr=NA, trunclevel=1, mar=NA, bicop=NA,
                 methods=c('kmeans'),  threshold=0.0001, maxit=10, cores=1, verbose=FALSE,
                 burn_in_iters=5, trunc_lvl=2, tau_threshold=0.1, batch_size=NULL, ema_alpha="decay", max_iter=500){
  initial_df_check(data)
  data <- as.matrix(data)
  initial_args_check(data, total_comp, is_cvine, vinestr, trunclevel, mar, bicop,
                     methods, threshold, maxit, cores)
  final_cvine <- is_cvine
  final_vinestr <- vinestr
  final_trunclevel <- NA
  final_mar <- mar
  final_bicop <- bicop
  final_bicop_mapped <- map_family(bicop)
  winner_bic <- 1000000

  
  global_min <- apply(data, 2, min)
  global_max <- apply(data, 2, max)
  global_sd <- apply(data, 2, sd)
  global_sd[global_sd < sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  
  progressr::with_progress({
    for(method in methods){
      initial_out <- initial_clustering(data, total_comp, is_cvine, vinestr, trunc_lvl, mar, bicop, method, tau_threshold, trunc_lvl_param=trunc_lvl, cores=cores)
      marginal_params <- initial_out$marginal_params
      marginal_fams <- initial_out$marginal_fams
      u_data <- initial_out$u_data
      vine_models <- initial_out$vine_models
      mix_probs <- initial_out$mix_probs
      total_obs <- dim(data)[1]
      total_features <- dim(data)[2]
      iteration <- 1
      loglik_res <- vector()
      smoothed_loglik_res <- vector()
      cond <- TRUE
      
      p_ecm <- progressr::progressor(steps = 1)
      
      prev_marginal_params <- NULL
      prev_vine_models <- NULL
      prev_mix_probs <- NULL

      use_batching <- !is.null(batch_size) && batch_size < total_obs
      actual_burn_in <- if(use_batching) max(burn_in_iters, ceiling(total_obs / batch_size) * 2) else burn_in_iters

      while(cond==TRUE){
        batch_indices <- if(use_batching) sample(1:total_obs, batch_size, replace = FALSE) else 1:total_obs
        data_batch <- data[batch_indices, , drop=FALSE]
        total_obs_batch <- length(batch_indices)
        
        # Calculate u_data on the fly for the batch to save memory/time
        u_data_batch <- array(0, dim=c(total_obs_batch, total_features, total_comp))
        for(j in 1:total_comp){
          u_data_batch[,,j] <- eval_all_margins_cpp(data_batch, marginal_fams[,j], marginal_params[,,j], "cdf")
        }
        
        rvine_densities <- matrix(0, total_obs_batch, total_comp)
        
        # Avoid 1/0 boundary evaluation problems for rvinecopulib
        u_data_safe <- u_data_batch
        u_data_safe[u_data_safe < 1e-10] <- 1e-10
        u_data_safe[u_data_safe > 1 - 1e-10] <- 1 - 1e-10
        
        rvine_densities <- sapply(1:total_comp, function(j) rvinecopulib::dvinecop(u_data_safe[,,j], vine_models[[j]], cores=cores))
        
        margin_densities <- array(0, dim=c(total_obs_batch, total_features, total_comp))
        for(j in 1:total_comp){
          margin_densities[,,j] <- eval_all_margins_cpp(data_batch, marginal_fams[,j], marginal_params[,,j], "pdf")
        }
        
        margin_densities[margin_densities < 1e-300] <- 1e-300
        rvine_densities[rvine_densities < 1e-300] <- 1e-300
        
        log_lik_points <- matrix(0, nrow=total_obs_batch, ncol=total_comp)
        for(j in 1:total_comp) {
          log_m_dens <- rowSums(log(matrix(margin_densities[,,j], nrow=total_obs_batch)))
          log_c_dens <- log(rvine_densities[,j])
          log_lik_points[,j] <- log(mix_probs[j]) + log_m_dens + log_c_dens
        }
        
        max_log_lik <- apply(log_lik_points, 1, max)
        exp_diff <- exp(log_lik_points - max_log_lik)
        sum_exp <- rowSums(exp_diff)
        
        lik_per_obs <- max_log_lik + log(sum_exp)
        loglik <- sum(lik_per_obs)
        if (use_batching) loglik <- loglik * (total_obs / total_obs_batch)
        loglik_res[iteration] <- loglik
        
        # Calculate dynamic alpha if decaying
        if (use_batching) {
            current_alpha <- if (is.numeric(ema_alpha)) ema_alpha else (iteration + 5)^(-0.6)
        } else {
            current_alpha <- 1.0 # Full batch EM overrides EMA and uses exactly 100% of the new iteration
        }
        
        if (iteration == 1) {
          smoothed_loglik_res[iteration] <- loglik
        } else {
          smoothed_loglik_res[iteration] <- (1 - current_alpha) * smoothed_loglik_res[iteration-1] + current_alpha * loglik
        }
        
        if (verbose) {
           if (use_batching) {
               message(sprintf("Method %s | ECM iteration %d | smoothed batch logLik %.4f", method, iteration, smoothed_loglik_res[iteration]))
           } else {
               message(sprintf("Method %s | ECM iteration %d | logLik %.4f", method, iteration, loglik))
           }
        }
        p_ecm(sprintf("Method %s: ECM Iteration %d (logLik %.4f)", method, iteration, loglik), amount = 0)
        
        patience <- if (use_batching) 10 else 1
        min_iters_required <- if (use_batching) max(actual_burn_in + patience, ceiling(total_obs / batch_size)) else 2
        
        if (iteration > min_iters_required) {
          loglik_to_check <- if(use_batching) smoothed_loglik_res else loglik_res
          if ((abs(loglik_to_check[iteration] - loglik_to_check[iteration - patience]) / abs(loglik_to_check[iteration - patience])) <= threshold){
            cond <- FALSE
            p_ecm("Converged!", amount = 1)
            break
          }
        }
        
        if (iteration >= max_iter) {
           cond <- FALSE
           p_ecm(sprintf("Max Iterations Reached (%d)", iteration), amount = 1)
           warning(sprintf("ECM algorithm reached maximum permitted %d iterations without convergence", max_iter))
           break
        }
        
        #E-step
        z_values <- exp_diff / sum_exp
        z_values[is.na(z_values)] <- 1 / total_comp
        z_values[z_values < 0] <- 0
        z_values[z_values > 1] <- 1
        
        #CM-steps:
        #CM-step 1
        mix_probs_new <- CM_step_mixture_probs(z_values)
        
        current_maxit <- if(iteration < actual_burn_in) 2 else maxit
        
        #CM-step 2 and 3
        CMS <- lapply(1:total_comp, function(x) try(
          CM_steps(data_batch, vine_models[[x]], z_values[,x], marginal_fams[,x], marginal_params[,,x],
                   current_maxit, final_bicop_mapped, iteration, actual_burn_in, trunc_lvl, tau_threshold, cores, global_min, global_max, global_sd),
          silent = TRUE
        ))
        
        failed_components <- which(vapply(CMS, inherits, logical(1), "try-error"))
        if(length(failed_components) > 0){
          failure_messages <- vapply(CMS[failed_components], function(err) conditionMessage(attr(err, "condition")), character(1))
          stop("CM-step failed for component(s) ",
               paste(failed_components, collapse = ", "),
               ": ",
               paste(unique(failure_messages), collapse = " | "))
        }
        
        if (use_batching && iteration > 1) {
           mix_probs <- (1 - current_alpha) * prev_mix_probs + current_alpha * mix_probs_new
           for(j in 1:total_comp) {
              marginal_params[,,j] <- (1 - current_alpha) * prev_marginal_params[,,j] + current_alpha * CMS[[j]]$marginal_par
              vine_models[[j]] <- CMS[[j]]$vine_model
           }
        } else {
           mix_probs <- mix_probs_new
           for(j in 1:total_comp){
              marginal_params[,,j] <- CMS[[j]]$marginal_par
              vine_models[[j]] <- CMS[[j]]$vine_model
           }
        }
        
        prev_mix_probs <- mix_probs
        prev_marginal_params <- marginal_params
        prev_vine_models <- vine_models
        iteration <- iteration + 1
      }
      
      iteration <- iteration - 1
      
      # For final selection, calculate full u_data and z_values
      u_data <- array(0, dim=c(total_obs, total_features, total_comp))
      for(j in 1:total_comp){
         u_data[,,j] <- eval_all_margins_cpp(data, marginal_fams[,j], marginal_params[,,j], "cdf")
      }
      u_data_safe <- u_data
      u_data_safe[u_data_safe < 1e-10] <- 1e-10
      u_data_safe[u_data_safe > 1 - 1e-10] <- 1 - 1e-10
      rvine_densities <- sapply(1:total_comp, function(j) rvinecopulib::dvinecop(u_data_safe[,,j], vine_models[[j]], cores=cores))
      margin_densities <- array(0, dim=c(total_obs, total_features, total_comp))
      for(j in 1:total_comp){
         margin_densities[,,j] <- eval_all_margins_cpp(data, marginal_fams[,j], marginal_params[,,j], "pdf")
      }
      margin_densities[margin_densities < 1e-300] <- 1e-300
      rvine_densities[rvine_densities < 1e-300] <- 1e-300
      
      log_lik_points <- matrix(0, nrow=total_obs, ncol=total_comp)
      for(j in 1:total_comp) {
        log_m_dens <- rowSums(log(matrix(margin_densities[,,j], nrow=total_obs)))
        log_c_dens <- log(rvine_densities[,j])
        log_lik_points[,j] <- log(mix_probs[j]) + log_m_dens + log_c_dens
      }
      
      max_log_lik <- apply(log_lik_points, 1, max)
      exp_diff <- exp(log_lik_points - max_log_lik)
      sum_exp <- rowSums(exp_diff)
      z_values <- exp_diff / sum_exp
      z_values[is.na(z_values)] <- 1 / total_comp
      z_values[z_values < 0] <- 0
      z_values[z_values > 1] <- 1
      
      loglik <- sum(max_log_lik + log(sum_exp))
      
      final_out <- final_selection(data, total_comp, mix_probs, z_values, iteration, method, marginal_fams, marginal_params, vine_models, loglik)
      
      vcmm_bic <- final_out$bic
      if(vcmm_bic < winner_bic){
        winner_bic <- vcmm_bic
        out <- final_out
        vcmm_class <- apply(out$z_values,1,function(x) which(x==max(x)))
      }
    }
  })
  
  winner_bic <- round(winner_bic, 0)
  out_list <- list("output"=out, "cluster"=vcmm_class)
  class(out_list) <- "vcmm_res"
  out_list
}

#' @export
print.vcmm_res <- function(x, ...) {
  fit_info(x)
  invisible(x)
}

#' @export
summary.vcmm_res <- function(object, ...) {
  list(
    margins = object$output$margin,
    marginal_pars = object$output$marginal_param,
    vine_models = object$output$vine_models,
    mixture_probs = object$output$mixture_prob
  )
}

#' @export
predict.vcmm_res <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$cluster)
  }
  
  newdata <- as.matrix(newdata)
  total_comp <- length(object$output$mixture_prob)
  total_obs <- nrow(newdata)
  total_features <- ncol(newdata)
  
  rvine_densities <- matrix(0, total_obs, total_comp)
  margin_densities <- array(0, dim=c(total_obs, total_features, total_comp))
  u_data <- array(0, dim=c(total_obs, total_features, total_comp))
  
  for(j in 1:total_comp){
    u_data[,,j] <- eval_all_margins_cpp(newdata, object$output$margin[,j], object$output$marginal_param[,,j], "cdf")
    margin_densities[,,j] <- eval_all_margins_cpp(newdata, object$output$margin[,j], object$output$marginal_param[,,j], "pdf")
    u_data_safe <- u_data[,,j]
    u_data_safe[u_data_safe < 1e-10] <- 1e-10
    u_data_safe[u_data_safe > 1 - 1e-10] <- 1 - 1e-10
    rvine_densities[,j] <- rvinecopulib::dvinecop(u_data_safe, object$output$vine_models[[j]])
  }
  
  margin_densities[margin_densities < 1e-300] <- 1e-300
  rvine_densities[rvine_densities < 1e-300] <- 1e-300
  
  log_lik_points <- matrix(0, nrow=total_obs, ncol=total_comp)
  for(j in 1:total_comp) {
    log_m_dens <- rowSums(log(matrix(margin_densities[,,j], nrow=total_obs)))
    log_c_dens <- log(rvine_densities[,j])
    log_lik_points[,j] <- log(object$output$mixture_prob[j]) + log_m_dens + log_c_dens
  }
  
  max_log_lik <- apply(log_lik_points, 1, max)
  class <- apply(log_lik_points, 1, function(x) which.max(x))
  class
}
