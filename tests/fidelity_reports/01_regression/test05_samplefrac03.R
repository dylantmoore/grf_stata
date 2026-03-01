## Test 05: sample.fraction=0.3
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

rf <- regression_forest(X, Y, num.trees = 500, seed = 42, sample.fraction = 0.3)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test05_data.csv", row.names = FALSE)
cat("Test05 R done. Predictions range:", range(preds), "\n")
