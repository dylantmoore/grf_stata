## Test 06: regression.splitting + multiple quantiles (0.1, 0.25, 0.5, 0.75, 0.9)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

qs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
qf <- quantile_forest(X, Y, quantiles = qs,
                      num.trees = 500, seed = 42, regression.splitting = TRUE)
preds <- predict(qf, quantiles = qs)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_q10 <- preds[, 1]
df$r_q25 <- preds[, 2]
df$r_q50 <- preds[, 3]
df$r_q75 <- preds[, 4]
df$r_q90 <- preds[, 5]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test06_data.csv", row.names = FALSE)
cat("Test06 R done (regsplit + 5 quantiles). q50 range:", range(preds[,3]), "\n")
