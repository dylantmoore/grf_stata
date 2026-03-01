## Test 10: With clusters
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3] + rnorm(n, sd=0.5)
# 50 clusters of 10 each
clusters <- rep(1:50, each=10)

brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42, clusters=clusters, boost.steps=2)
preds <- predict(brf)$predictions
cat("Test10 R done. clusters=50 groups, boost.steps=2\n")
cat("  Predictions range:", range(preds), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
df$cluster_id <- clusters
df$true_mu <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test10_data.csv", row.names=FALSE)
cat("  R-vs-true corr:", cor(preds, df$true_mu), "\n")
