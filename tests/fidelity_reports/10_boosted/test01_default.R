## Test 01: Default options (auto-tune boost steps)
## DGP: Y = sin(2*X1) + X2^2 + 0.5*X3 + N(0,0.5)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3] + rnorm(n, sd=0.5)

# booststeps=NULL => auto-tune
brf <- boosted_regression_forest(X, Y, num.trees=500, seed=42)
preds <- predict(brf)$predictions
boost_steps_used <- length(brf$forests)
cat("Test01 R done.\n")
cat("  Boost steps used:", boost_steps_used, "\n")
cat("  Predictions range:", range(preds), "\n")
cat("  Predictions mean:", mean(preds), "\n")
cat("  Predictions sd:", sd(preds), "\n")

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
df$true_mu <- sin(X[,1]*2) + X[,2]^2 + 0.5*X[,3]
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test01_data.csv", row.names=FALSE)
cat("  R-vs-true corr:", cor(preds, df$true_mu), "\n")
cat("  R boost_steps:", boost_steps_used, "\n")
