## Test 14: Compare boosted vs plain regression_forest (OOB MSE)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3] + rnorm(n, sd=0.5)
true_mu <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3]

# Plain regression forest
rf <- regression_forest(X, Y, num.trees=500, seed=42)
rf_preds <- predict(rf)$predictions
rf_mse <- mean((rf_preds - true_mu)^2)

# Boosted regression forest (auto-tune)
brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42)
brf_preds <- predict(brf)$predictions
brf_mse <- mean((brf_preds - true_mu)^2)

cat("Test14: Boosted vs plain regression forest\n")
cat("  Regression forest OOB MSE (vs true mu):", rf_mse, "\n")
cat("  Boosted forest OOB MSE (vs true mu):", brf_mse, "\n")
cat("  Boost steps used:", brf$boost.steps, "\n")
cat("  Boosted improvement:", (rf_mse - brf_mse)/rf_mse * 100, "%\n")
cat("  rf-vs-true corr:", cor(rf_preds, true_mu), "\n")
cat("  brf-vs-true corr:", cor(brf_preds, true_mu), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred_brf <- brf_preds
df$r_pred_rf <- rf_preds
df$true_mu <- true_mu
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test14_data.csv", row.names=FALSE)
