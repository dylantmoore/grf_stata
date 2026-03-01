## Test 02: 3-class multinomial
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
preds <- predict(rf)$predictions   # n x 3

df <- as.data.frame(X)
df$y <- Y3
df$r_pred_c0 <- preds[, 1]
df$r_pred_c1 <- preds[, 2]
df$r_pred_c2 <- preds[, 3]

write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/06_probability/test02_data.csv", row.names = FALSE)
cat("Test02 3-class R done.\n")
cat("  Y distribution:", table(Y3), "\n")
cat("  Sum check max |sum-1|:", max(abs(rowSums(preds) - 1)), "\n")
