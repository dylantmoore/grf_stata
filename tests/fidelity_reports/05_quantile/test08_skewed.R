## Test 08: Skewed data — Y = exp(X1) + rnorm(n), right-skewed outcome
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- exp(X[,1]) + rnorm(n, sd = 0.3)

qf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9), num.trees = 500, seed = 42)
preds <- predict(qf, quantiles = c(0.1, 0.5, 0.9))$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_q10 <- preds[, 1]
df$r_q50 <- preds[, 2]
df$r_q90 <- preds[, 3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test08_data.csv", row.names = FALSE)
cat("Test08 R done (skewed). q50 range:", range(preds[,2]), "\n")
