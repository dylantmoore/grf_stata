## Test 11: Combined options: cluster + weights + nohonesty
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
prob <- pnorm(X[,1] + 0.5 * X[,2])
Y <- rbinom(n, 1, prob)
Y_factor <- as.factor(Y)

clusters <- rep(1:50, each = 10)
set.seed(123)
weights <- runif(n, 0.5, 2.0)

rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42,
                         clusters = clusters,
                         sample.weights = weights,
                         honesty = FALSE)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$cluster_id <- clusters
df$wt <- weights
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test11_data.csv", row.names = FALSE)
cat("Test11 Combined R done.\n")
cat("  Sum check max |sum-1|:", max(abs(rowSums(preds) - 1)), "\n")
