## Test 13: nomia (no missing indicator approach)
## Note: grf 2.5.0 does NOT have enable.missing.indicator parameter.
## MIA is always on in R when NAs are present. In Stata, nomia disables it.
## This test compares: R (MIA always on) vs Stata nomia (MIA off) on data WITHOUT NAs.
## When data has no NAs, MIA on vs off should give identical results.
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

# No missing data: MIA on vs off should be equivalent
rf <- regression_forest(X, Y, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test13_data.csv", row.names = FALSE)
cat("Test13 R done (no NAs, MIA irrelevant). Predictions range:", range(preds), "\n")
cat("NOTE: grf 2.5.0 has no enable.missing.indicator param. MIA always on in R.\n")
