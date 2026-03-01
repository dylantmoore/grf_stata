## Test 04: nclasses specified vs auto-detected (binary)
## Uses same binary DGP as test01, compares explicit nclasses(2) vs auto
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
prob <- pnorm(X[,1] + 0.5 * X[,2])
Y <- rbinom(n, 1, prob)
Y_factor <- as.factor(Y)

# R always auto-detects: use same predictions
rf <- probability_forest(X, Y_factor, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test04_data.csv", row.names = FALSE)
cat("Test04 nclasses specified R done.\n")
