## Test 19: Missing data (NAs in X), MIA on (default)
library(grf)
set.seed(42)
n <- 500; p <- 5
X <- matrix(runif(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- X[,1] + 2 * X[,2] + rnorm(n)

# Introduce ~10% missing in x1 and x3
miss_idx1 <- sample(n, 50)
miss_idx3 <- sample(n, 50)
X[miss_idx1, 1] <- NA
X[miss_idx3, 3] <- NA

rf <- regression_forest(X, Y, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions

df <- as.data.frame(X)
df$y <- Y
df$r_pred <- preds

# Write with explicit NA string as "." for Stata, but Stata handles missing numerics
# Use write.table with na="" for empty cells that Stata reads as .
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/01_regression/test19_data.csv",
          row.names = FALSE, na = "")
cat("Test19 R done. Predictions range:", range(preds), "\n")
cat("NAs in x1:", sum(is.na(df$x1)), "\n")
cat("NAs in x3:", sum(is.na(df$x3)), "\n")
