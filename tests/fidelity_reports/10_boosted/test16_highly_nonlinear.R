## Test 16: Highly nonlinear data Y = sin(3*X1)*cos(2*X2) - boosting should help
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(3*X[,1])*cos(2*X[,2]) + rnorm(n, sd=0.3)
true_mu <- sin(3*X[,1])*cos(2*X[,2])

# Plain regression forest
rf <- regression_forest(X, Y, num.trees=500, seed=42)
rf_preds <- predict(rf)$predictions
rf_mse <- mean((rf_preds - true_mu)^2)

# Boosted
brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42)
brf_preds <- predict(brf)$predictions
brf_mse <- mean((brf_preds - true_mu)^2)

cat("Test16: Highly nonlinear data\n")
cat("  RF MSE:", rf_mse, "  BRF MSE:", brf_mse, "\n")
cat("  Boost steps:", brf$boost.steps, "\n")
cat("  Improvement:", (rf_mse - brf_mse)/rf_mse * 100, "%\n")
cat("  rf-vs-true corr:", cor(rf_preds, true_mu), "\n")
cat("  brf-vs-true corr:", cor(brf_preds, true_mu), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- brf_preds
df$true_mu <- true_mu
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test16_data.csv", row.names=FALSE)
