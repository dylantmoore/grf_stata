## Test 12: Probabilities sum to 1 (3-class)
## Verify that class probabilities sum to 1 within tolerance
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)

eta1 <- X[,1] + X[,2]
eta2 <- -X[,1] + X[,3]
p1 <- exp(eta1) / (1 + exp(eta1) + exp(eta2))
p2 <- exp(eta2) / (1 + exp(eta1) + exp(eta2))
p0 <- 1 - p1 - p2
Y3 <- sapply(1:n, function(i) sample(0:2, 1, prob = c(p0[i], p1[i], p2[i])))
Y_factor <- as.factor(Y3)

rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions

cat("R sum check:\n")
cat("  max |sum - 1| =", max(abs(rowSums(preds) - 1)), "\n")
cat("  All within 0.01?", all(abs(rowSums(preds) - 1) < 0.01), "\n")

df <- as.data.frame(X)
df$y <- Y3
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]
df$r_pred_c2 <- preds[, 3]
df$r_sum <- rowSums(preds)

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test12_data.csv", row.names = FALSE)
cat("Test12 ProbSum R done.\n")
