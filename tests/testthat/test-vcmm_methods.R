RVMs <- list()
dims <- 3
RVMs[[1]] <- rvinecopulib::vinecop_dist(
  pair_copulas = list(
    list(rvinecopulib::bicop_dist(family = "clayton", rotation = 0, parameters = 0.5),
         rvinecopulib::bicop_dist(family = "gumbel", rotation = 0, parameters = 2.5)),
    list(rvinecopulib::bicop_dist(family = "gumbel", rotation = 180, parameters = 5))
  ),
  structure = rvinecopulib::dvine_structure(1:3)
)
RVMs[[2]] <- rvinecopulib::vinecop_dist(
  pair_copulas = list(
    list(rvinecopulib::bicop_dist(family = "joe", rotation = 0, parameters = 2),
         rvinecopulib::bicop_dist(family = "frank", rotation = 0, parameters = 14)),
    list(rvinecopulib::bicop_dist(family = "clayton", rotation = 180, parameters = 1))
  ),
  structure = rvinecopulib::dvine_structure(1:3)
)
margin <- matrix(c('Normal', 'Gamma', 'Lognormal', 'Lognormal', 'Normal', 'Student-t'), 3, 2)
margin_pars <- array(0, dim=c(4, 3, 2))
margin_pars[,1,1] <- c(1, 2, 0, 0)
margin_pars[,1,2] <- c(1.5, 0.4, 0, 0)
margin_pars[,2,1] <- c(1, 0.2, 0, 0)
margin_pars[,2,2] <- c(18, 5, 0, 0)
margin_pars[,3,1] <- c(0.8, 0.8, 0, 0)
margin_pars[,3,2] <- c(4, 2, 5, 0)

# simulate data
set.seed(111)
x_data <- rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs)

test_that("d/r vcmm work", {
  expect_false(all(rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs) ==
                     rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs)))
  set.seed(111)
  expect_true(all(x_data == rvcmm(dims=3, obs=c(100,200),
                                  margin=margin, margin_pars=margin_pars, RVMs=RVMs)))
  expect_gte(min(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3,2/3))), 0)
})

test_that("catches incorrect arguments", {
  expect_error(dvcmm(x_data, margin[,1], margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars[,,1], RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs[[1]], c(1/3, 2/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3, 1/3, 1/3)))
  expect_error(dvcmm(x_data, margin, margin_pars, RVMs, c(1/3, 1)))
  expect_error(dvcmm(x_data,  matrix(c('xyz', 'xyz', 'xyz', 'xyz', 'xyz', 'xyz'), 3, 2),
                     margin_pars, RVMs, c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,NA),2,3), margin, margin_pars, RVMs,
                     c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,'x'),2,3), margin, margin_pars, RVMs,
                     c(1/3, 2/3)))
  expect_error(dvcmm(matrix(c(1,2,3,1,2,3),2,3), margin, margin_pars, RVMs, 1))
  expect_error(rvcmm(dims='x', obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=999, obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200,300), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=1, obs=c(100,200), margin, margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin[,1], margin_pars, RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin, margin_pars[,,1], RVMs))
  expect_error(rvcmm(dims=3, obs=c(100,200), margin, margin_pars, RVMs[[1]]))
  expect_error(rvcmm(dims=3, obs=c(100,200),  matrix(c('xyz', 'xyz', 'xyz', 'xyz', 'xyz', 'xyz'), 3, 2),
                     margin_pars, RVMs))
})
