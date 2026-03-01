## Test 05: boost.max.steps=10 (allow more auto-tune steps)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3] + rnorm(n, sd=0.5)

brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42, boost.max.steps=10)
preds <- predict(brf)$predictions
cat("Test05 R done. boost.max.steps=10\n")
cat("  Predictions range:", range(preds), "\n")
cat("  Actual boost steps:", brf$boost.steps, "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
df$true_mu <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test05_data.csv", row.names=FALSE)
cat("  R-vs-true corr:", cor(preds, df$true_mu), "\n")
