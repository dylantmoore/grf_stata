## Test 16: Combined: cluster + weights + estimatevariance
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)
cluster_id <- rep(1:10, length.out = n)
wts <- runif(n) + 0.1

rf <- regression_forest(X, Y, num.trees = 500, seed = 42,
                        clusters = cluster_id, sample.weights = wts,
                        ci.group.size = 2)
result <- predict(rf, estimate.variance = TRUE)
preds <- result$predictions
vars  <- result$variance.estimates

df <- as.data.frame(X)
df$y <- Y
df$cluster_id <- cluster_id
df$wts <- wts
df$r_pred <- preds
df$r_var  <- vars
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test16_data.csv", row.names = FALSE)
cat("Test16 R done. Predictions range:", range(preds), "\n")
cat("Test16 R done. Variance range:", range(vars), "\n")
