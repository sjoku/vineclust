# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("final model selection works properly", {
  z_vals <- matrix(c(rep(0.9, 100), rep(0.1, 100), rep(0.1, 100), rep(0.9, 100)), 200, 2)
  mfams <- matrix(c("Normal", "Normal", "Normal", "Normal"), 2, 2)
  mpars <- array(0, dim=c(4, 2, 2))
  vmodels <- list(list(npars=1), list(npars=1))
  fit <- final_selection(data=x_data, total_cluster=2, mix_probs=c(0.5, 0.5), z_values=z_vals,
                         iteration=1, init_method="gmm", marginal_fams=mfams, marginal_params=mpars, vine_models=vmodels, loglik=-100)
  expect_identical(
    names(fit),
    c(
      "loglik", "bic", "icl", "init_clustering", "iteration", "total_pars",
      "mixture_prob", "margin","marginal_param", "vine_models", "z_values"
    )
  )
  expect_type(fit$loglik, "double")
  expect_type(fit$bic, "double")
  expect_type(fit$icl, "double")
  expect_type(fit$init_clustering, "character")
  expect_type(fit$iteration, "double")
  expect_type(fit$total_pars, "double")
  expect_type(fit$mixture_prob, "double")
  expect_type(fit$margin, "character")
  expect_type(fit$marginal_param, "double")
  expect_type(fit$vine_models, "list")
  expect_type(fit$z_values, "double")
})
