if (!requireNamespace("rvinecopulib", quietly = TRUE)) install.packages("rvinecopulib", repos = "http://cran.us.r-project.org")
if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr", repos = "http://cran.us.r-project.org")
if (!requireNamespace("mclust", quietly = TRUE)) install.packages("mclust", repos = "http://cran.us.r-project.org")
if (!requireNamespace("univariateML", quietly = TRUE)) install.packages("univariateML", repos = "http://cran.us.r-project.org")
if (!requireNamespace("fGarch", quietly = TRUE)) install.packages("fGarch", repos = "http://cran.us.r-project.org")

files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
for (f in files) {
  tryCatch(source(f), error = function(e) message("Error sourcing ", f, ": ", e$message))
}

cat("Files sourced successfully. Running migration test on Wisconsin Breast Cancer Data...\n")
data_wisc <- read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data", header = FALSE)

cat("Fitting VCMM...\n")
set.seed(42)

tryCatch({
  fit <- vcmm(data = data_wisc[1:100, c(15,27,29,30)], total_comp=2, maxit=2, verbose=TRUE)
  cat("Success! Results:\n")
  print(fit)
}, error = function(e) {
  cat("Migration failed during execution. Error:\n", e$message, "\n")
})
