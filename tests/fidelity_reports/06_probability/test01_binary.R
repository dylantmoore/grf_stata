## Test 01: Binary classification (2 classes)
## Y = I(X1 + X2 > 0), logistic-type DGP
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
prob <- pnorm(X[,1] + 0.5 * X[,2])
Y <- rbinom(n, 1, prob)
Y_factor <- as.factor(Y)

rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions   # n x 2 matrix: col 1 = P(Y=0), col 2 = P(Y=1)

df <- as.data.frame(X)
df$y <- Y
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test01_data.csv", row.names = FALSE)
cat("Test01 Binary R done.\n")
cat("  Pred c0 range:", range(preds[,1]), "\n")
cat("  Pred c1 range:", range(preds[,2]), "\n")
cat("  Sum check max |sum-1|:", max(abs(rowSums(preds) - 1)), "\n")
cat("  Y distribution:", table(Y), "\n")
