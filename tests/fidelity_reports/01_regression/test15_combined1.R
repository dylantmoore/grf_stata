## Test 15: Combined: nohonesty + mtry=3 + minnodesize=15
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

rf <- regression_forest(X, Y, num.trees = 500, seed = 42,
                        honesty = FALSE, mtry = 3, min.node.size = 15)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test15_data.csv", row.names = FALSE)
cat("Test15 R done. Predictions range:", range(preds), "\n")
