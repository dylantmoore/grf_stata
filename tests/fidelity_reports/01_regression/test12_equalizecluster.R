## Test 12: equalize.cluster.weights
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)
cluster_id <- rep(1:10, length.out = n)

rf <- regression_forest(X, Y, num.trees = 500, seed = 42,
                        clusters = cluster_id, equalize.cluster.weights = TRUE)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$cluster_id <- cluster_id
df$r_pred <- preds
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test12_data.csv", row.names = FALSE)
cat("Test12 R done. Predictions range:", range(preds), "\n")
