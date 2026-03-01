## Test 08: boost.trees.tune=50 (more trees per tuning step)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3] + rnorm(n, sd=0.5)

brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42, boost.trees.tune=50)
preds <- predict(brf)$predictions
cat("Test08 R done. boost.trees.tune=50\n")
cat("  Actual boost steps:", brf$boost.steps, "\n")
cat("  Predictions range:", range(preds), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
df$true_mu <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test08_data.csv", row.names=FALSE)
cat("  R-vs-true corr:", cor(preds, df$true_mu), "\n")
