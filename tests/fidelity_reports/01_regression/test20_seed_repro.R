## Test 20: Seed reproducibility (same seed → same predictions)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

# Run twice with same seed
rf1 <- regression_forest(X, Y, num.trees = 500, seed = 123)
rf2 <- regression_forest(X, Y, num.trees = 500, seed = 123)
preds1 <- predict(rf1)$predictions
preds2 <- predict(rf2)$predictions

cat("R internal reproducibility (same seed, corr):", cor(preds1, preds2), "\n")
cat("R internal exact match:", all(preds1 == preds2), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds1
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_data.csv", row.names = FALSE)

# Save seed check result separately
seed_check <- data.frame(r_corr = cor(preds1, preds2), exact_match = as.integer(all(preds1 == preds2)))
write.csv(seed_check, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_r_seed_check.csv", row.names = FALSE)
cat("Test20 R done. Predictions range:", range(preds1), "\n")
