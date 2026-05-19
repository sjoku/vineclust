# simulate data
set.seed(111)
cl1 <- sapply(1:2, function(i) rnorm(100, -10, 5))
cl2 <- sapply(1:2, function(i) rnorm(100, 10, 5))
x_data <- data.frame(rbind(cl1, cl2))

test_that("CM-step 1 works", {
  fit <- CM_step_mixture_probs(matrix(c(rep(0.9, 100), rep(0.1, 100), rep(0.1, 100), rep(0.9, 100)), 200, 2))
  expect_type(fit, "double")
  expect_lte(max(fit), 1)
  expect_gte(min(fit), 0)
  expect_equal(sum(fit), 1)
})


test_that("CM-steps 2 and 3 work", {
  vine_model <- rvinecopulib::vinecop_dist(
    pair_copulas = list(list(rvinecopulib::bicop_dist("gaussian", 0, 0.7))),
    structure = rvinecopulib::dvine_structure(1:2)
  )
  global_min <- apply(x_data, 2, min)
  global_max <- apply(x_data, 2, max)
  global_sd <- apply(x_data, 2, sd)
  
  fit <- CM_steps(data=x_data, vine_model=vine_model,
                  z_value=c(rep(0.9, 120), rep(0.1, 80)),
                  marginal_fam=c('Normal', 'Skew Normal'),
                  marginal_par=matrix(c(-9, 4.5, 0, 0, -3, 2, 4, 0), 4, 2), 
                  maxit=10, bicop_mapped="all", iteration=1, burn_in_iters=5, trunc_lvl=2, tau_threshold=0.1, 
                  cores=1, global_min=global_min, global_max=global_max, global_sd=global_sd)
  expect_identical(
    names(fit),
    c(
      "marginal_par", "vine_model", "u_data"
    )
  )
  expect_type(fit$marginal_par, "double")
  expect_s3_class(fit$vine_model, "vinecop_dist")
  expect_type(fit$u_data, "double")
})
