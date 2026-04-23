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
                 methods=c('kmeans'),  threshold=0.0001, maxit=10, cores=1, verbose=FALSE){
  initial_df_check(data)
  initial_args_check(data, total_comp, is_cvine, vinestr, trunclevel, mar, bicop,
                     methods, threshold, maxit, cores)
  final_cvine <- is_cvine
  final_vinestr <- vinestr
  final_trunclevel <- NA
  final_mar <- mar
  final_bicop <- bicop
  final_bicop_mapped <- map_family(bicop)
  winner_bic <- 1000000
  use_future <- cores > 1 && total_comp > 1
  if(use_future){
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    if(.Platform$OS.type == "windows"){
      future::plan(future::multisession, workers = min(cores, total_comp))
    }
    else{
      future::plan(future::multicore, workers = min(cores, total_comp))
    }
  }
  
  progressr::with_progress({
    for(method in methods){
      initial_out <- initial_clustering(data, total_comp, is_cvine, vinestr, trunclevel, mar, bicop, method)
      marginal_params <- initial_out$marginal_params
      marginal_fams <- initial_out$marginal_fams
      u_data <- initial_out$u_data
      vine_models <- initial_out$vine_models
      mix_probs <- initial_out$mix_probs
      total_obs <- dim(data)[1]
      total_features <- dim(data)[2]
      iteration <- 1
      loglik_res <- vector()
      cond <- TRUE
      
      p_ecm <- progressr::progressor(steps = 1)
      
      while(cond==TRUE){
        rvine_densities <- matrix(0, total_obs, total_comp)
        
        # Avoid 1/0 boundary evaluation problems for rvinecopulib
        u_data_safe <- pmax(pmin(u_data, 1 - 1e-10), 1e-10)
        
        rvine_densities <- sapply(1:total_comp, function(j) rvinecopulib::dvinecop(u_data_safe[,,j], vine_models[[j]]))
        
        margin_densities <- array(0, dim=c(total_obs, total_features, total_comp))
        for(j in 1:total_comp){
          margin_densities[,,j] <- sapply(1:total_features, function(x) pdf_cdf_quant_margin(data[,x],marginal_fams[x,j],
                                                                                           marginal_params[,x,j], 'pdf'))
        }
        
        # Vectorized prod
        total_margin_dens <- sapply(1:total_comp, function(j) apply(margin_densities[,,j], 1, prod))
        
        lik_points <- sapply(1:total_comp, function(j) mix_probs[j]*total_margin_dens[,j]*rvine_densities[,j])
        lik_per_obs <- apply(lik_points, 1, sum)
        lik_per_obs[which(lik_per_obs == 0 | is.na(lik_per_obs))] <- 1e-100
        loglik <- sum(log(lik_per_obs))
        loglik_res[iteration] <- loglik
        
        if (verbose) message(sprintf("Method %s | ECM iteration %d | logLik %.4f", method, iteration, loglik))
        p_ecm(sprintf("Method %s: ECM Iteration %d (logLik %.4f)", method, iteration, loglik), amount = 0)
        
        if(iteration > 2){
          if ((abs(loglik_res[iteration]-loglik_res[iteration-1])/abs(loglik_res[iteration-1])) <= threshold){
            cond <- FALSE
            p_ecm("Converged!", amount = 1)
            break
          }
        }
        
        if (iteration >= 500) {
           cond <- FALSE
           p_ecm(sprintf("Max Iterations Reached (%d)", iteration), amount = 1)
           warning("ECM algorithm reached maximum permitted 500 iterations without convergence")
           break
        }
        
        #E-step
        z_values <- t(apply(lik_points, 1, function(row) row / sum(row)))
        z_values[is.na(z_values)] <- 1 / total_comp
        
        #CM-steps:
        #CM-step 1
        mix_probs <- CM_step_mixture_probs(z_values)
        
        #CM-step 2 and 3
        if(use_future){
          CMS <- future.apply::future_lapply(1:total_comp, function(x) try(
            CM_steps(data, vine_models[[x]], z_values[,x], marginal_fams[,x], marginal_params[,,x],
                     maxit, final_bicop_mapped),
            silent = TRUE
          ), future.scheduling = 1)
        }
        else{
          CMS <- lapply(1:total_comp, function(x) try(
            CM_steps(data, vine_models[[x]], z_values[,x], marginal_fams[,x], marginal_params[,,x],
                     maxit, final_bicop_mapped),
            silent = TRUE
          ))
        }
        
        failed_components <- which(vapply(CMS, inherits, logical(1), "try-error"))
        if(length(failed_components) > 0){
          failure_messages <- vapply(CMS[failed_components], function(err) conditionMessage(attr(err, "condition")), character(1))
          stop("CM-step failed for component(s) ",
               paste(failed_components, collapse = ", "),
               ": ",
               paste(unique(failure_messages), collapse = " | "))
        }
        
        for(j in 1:total_comp){
          marginal_params[,,j] <- CMS[[j]]$marginal_par
          vine_models[[j]] <- CMS[[j]]$vine_model
          u_data[,,j] <- CMS[[j]]$u_data
        }
        iteration <- iteration + 1
      }
      
      iteration <- iteration - 1
      mix_probs <- CM_step_mixture_probs(z_values)
      final_out <- final_selection(data, total_comp, final_cvine, final_vinestr, final_trunclevel, mix_probs, z_values,
                                   iteration, method, final_mar, final_bicop)
      
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
    u_data[,,j] <- sapply(1:total_features, function(i) pdf_cdf_quant_margin(newdata[,i], object$output$margin[i,j],
                                                                             object$output$marginal_param[,i,j], 'cdf'))
    margin_densities[,,j] <- sapply(1:total_features, function(i) pdf_cdf_quant_margin(newdata[,i], object$output$margin[i,j],
                                                                                       object$output$marginal_param[,i,j], 'pdf'))
    u_data_safe <- pmax(pmin(u_data[,,j], 1 - 1e-10), 1e-10)
    rvine_densities[,j] <- rvinecopulib::dvinecop(u_data_safe, object$output$vine_models[[j]])
  }
  
  lik_points <- sapply(1:total_comp, function(j) object$output$mixture_prob[j] * apply(margin_densities[,,j, drop=FALSE], 1, prod) * rvine_densities[,j])
  
  if (total_obs == 1) {
    lik_points <- matrix(lik_points, nrow=1)
  }
  
  z_values <- t(apply(lik_points, 1, function(row) row / sum(row)))
  z_values[is.na(z_values)] <- 1 / total_comp
  
  class <- apply(z_values, 1, function(x) which.max(x))
  class
}
