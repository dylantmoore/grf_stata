## Test 09: estimate.variance
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

rf <- regression_forest(X, Y, num.trees = 500, seed = 42, ci.group.size = 2)
result <- predict(rf, estimate.variance = TRUE)
preds <- result$predictions
vars  <- result$variance.estimates

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
df$r_var  <- vars
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test09_data.csv", row.names = FALSE)
cat("Test09 R done. Variance range:", range(vars), "\n")
