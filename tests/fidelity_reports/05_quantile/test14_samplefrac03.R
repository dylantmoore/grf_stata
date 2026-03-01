## Test 14: samplefrac=0.3 — smaller subsample
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

qf <- quantile_forest(X, Y, quantiles = c(0.1, 0.5, 0.9),
                      num.trees = 500, seed = 42, sample.fraction = 0.3)
preds <- predict(qf, quantiles = c(0.1, 0.5, 0.9))$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_q10 <- preds[, 1]
df$r_q50 <- preds[, 2]
df$r_q90 <- preds[, 3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test14_data.csv", row.names = FALSE)
cat("Test14 R done (samplefrac=0.3). q50 range:", range(preds[,2]), "\n")
