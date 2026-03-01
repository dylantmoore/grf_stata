## Test 02: Single quantile (0.5 - median only)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

qf <- quantile_forest(X, Y, quantiles = c(0.5), num.trees = 500, seed = 42)
preds <- predict(qf, quantiles = c(0.5))$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_q50 <- preds[, 1]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test02_data.csv", row.names = FALSE)
cat("Test02 R done. q50 range:", range(preds[,1]), "\n")
