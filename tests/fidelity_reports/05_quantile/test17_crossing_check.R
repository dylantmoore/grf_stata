## Test 17: Quantile crossing check — do q10 <= q50 <= q90 hold?
## Uses default data; this test verifies ordering in both R and Stata outputs
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

qf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 500, seed = 42)
preds <- predict(qf, quantiles = c(0.1, 0.5, 0.9))$predictions

# Check ordering in R
n_cross_10_50 <- sum(preds[,1] > preds[,2])
n_cross_50_90 <- sum(preds[,2] > preds[,3])
cat("Test17 (crossing check):\n")
cat("  q10 > q50 violations:", n_cross_10_50, "\n")
cat("  q50 > q90 violations:", n_cross_50_90, "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_q10 <- preds[, 1]
df$r_q50 <- preds[, 2]
df$r_q90 <- preds[, 3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test17_data.csv", row.names = FALSE)
cat("Test17 R done. r_cross_10_50=", n_cross_10_50, " r_cross_50_90=", n_cross_50_90, "\n")
