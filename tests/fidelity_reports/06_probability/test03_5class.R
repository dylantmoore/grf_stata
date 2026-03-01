## Test 03: 5-class multinomial
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)

# 5-class via softmax on 4 linear combinations
eta <- cbind(
  rep(0, n),
  X[,1] + 0.5 * X[,2],
  -X[,1] + 0.5 * X[,3],
  0.5 * X[,2] - X[,4],
  X[,1] - X[,3] + 0.3 * X[,5]
)
probs <- exp(eta) / rowSums(exp(eta))
Y5 <- apply(probs, 1, function(p) sample(0:4, 1, prob = p))
Y_factor <- as.factor(Y5)

rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions   # n x 5

df <- as.data.frame(X)
df$y <- Y5
for (c in 0:4) {
  df[[paste0("r_pred_c", c)]] <- preds[, c + 1]
}

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test03_data.csv", row.names = FALSE)
cat("Test03 5-class R done.\n")
cat("  Y distribution:", table(Y5), "\n")
cat("  Sum check max |sum-1|:", max(abs(rowSums(preds) - 1)), "\n")
