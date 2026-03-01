## Test 04: Extreme quantiles (0.01, 0.99)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

qs <- c(0.01, 0.99)
qf <- quantile_forest(X, Y, quantiles = qs, num.trees = 500, seed = 42)
preds <- predict(qf, quantiles = qs)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_q1  <- preds[, 1]
df$r_q99 <- preds[, 2]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test04_data.csv", row.names = FALSE)
cat("Test04 R done. q1 range:", range(preds[,1]), " q99 range:", range(preds[,2]), "\n")
