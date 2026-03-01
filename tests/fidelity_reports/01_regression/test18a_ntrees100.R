## Test 18a: ntrees=100
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

rf <- regression_forest(X, Y, num.trees = 100, seed = 42)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test18a_data.csv", row.names = FALSE)
cat("Test18a R done. Predictions range:", range(preds), "\n")
