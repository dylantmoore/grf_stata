## Test 08: nohonesty (honesty = FALSE)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
prob <- pnorm(X[,1] + 0.5 * X[,2])
Y <- rbinom(n, 1, prob)
Y_factor <- as.factor(Y)

rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42,
                         honesty = FALSE)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test08_data.csv", row.names = FALSE)
cat("Test08 NoHonesty R done.\n")
cat("  Sum check max |sum-1|:", max(abs(rowSums(preds) - 1)), "\n")
