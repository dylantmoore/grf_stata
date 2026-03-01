## Comprehensive R fidelity tests for survival_forest and grf_expected_survival
## Generates data CSV files for comparison with Stata
library(grf)

OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/07_survival"

## ---- DGP ----
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
T_true <- rexp(n, rate = exp(0.5 * X[, 1]))   # Cox-type hazard
C       <- rexp(n, rate = 0.3)                 # censoring time
Y       <- pmin(T_true, C)
D       <- as.integer(T_true <= C)

cat("DGP: n=", n, "events=", sum(D), "censored=", sum(1-D),
    "event_rate=", round(mean(D), 3), "\n")

## Helper: integrate survival curve -> expected survival
expected_survival_r <- function(sf) {
  surv_pred   <- predict(sf)$predictions        # n x K
  ft          <- sf$failure.times               # K failure times
  K           <- length(ft)
  # Trapezoidal: E[T|X] = S(0)*t1 + sum_{j=2}^K 0.5*(S(t_{j-1})+S(t_j))*(t_j-t_{j-1})
  # i.e., [0,t1] interval: 0.5*(1 + S(t1))*t1
  esurv <- 0.5 * (1 + surv_pred[, 1]) * ft[1]
  if (K > 1) {
    for (j in 2:K) {
      dt <- ft[j] - ft[j - 1]
      esurv <- esurv + 0.5 * (surv_pred[, j - 1] + surv_pred[, j]) * dt
    }
  }
  esurv
}

## Helper: write test data frame
write_test <- function(df, name) {
  fpath <- file.path(OUTDIR, paste0(name, "_data.csv"))
  write.csv(df, fpath, row.names = FALSE)
  cat("Written:", fpath, "\n")
}

## ============================================================
## TEST 01: Default – noutput=20, Kaplan-Meier
## ============================================================
cat("\n=== TEST 01: Default (noutput=20, KM) ===\n")
sf01 <- survival_forest(X, Y, D, num.trees = 500, seed = 42)
pred01 <- predict(sf01)$predictions  # n x K_actual
ft01 <- sf01$failure.times
K01 <- min(20, length(ft01))
cat("failure times:", length(ft01), "| using first", K01, "\n")

df01 <- as.data.frame(X)
df01$time <- Y; df01$status <- D
for (j in 1:K01) df01[[paste0("r_s", j)]] <- pred01[, j]
write_test(df01, "test01_default")

## ============================================================
## TEST 02: noutput=50
## ============================================================
cat("\n=== TEST 02: noutput=50 ===\n")
sf02 <- survival_forest(X, Y, D, num.trees = 500, seed = 42)
pred02 <- predict(sf02)$predictions
ft02 <- sf02$failure.times
K02 <- min(50, length(ft02))
cat("failure times:", length(ft02), "| using first", K02, "\n")

df02 <- as.data.frame(X)
df02$time <- Y; df02$status <- D
for (j in 1:K02) df02[[paste0("r_s", j)]] <- pred02[, j]
write_test(df02, "test02_noutput50")

## ============================================================
## TEST 03: noutput=100
## ============================================================
cat("\n=== TEST 03: noutput=100 ===\n")
sf03 <- survival_forest(X, Y, D, num.trees = 500, seed = 42)
pred03 <- predict(sf03)$predictions
ft03 <- sf03$failure.times
K03 <- min(100, length(ft03))
cat("failure times:", length(ft03), "| using first", K03, "\n")

df03 <- as.data.frame(X)
df03$time <- Y; df03$status <- D
for (j in 1:K03) df03[[paste0("r_s", j)]] <- pred03[, j]
write_test(df03, "test03_noutput100")

## ============================================================
## TEST 04: predtype=0 (Nelson-Aalen)
## ============================================================
cat("\n=== TEST 04: predtype=0 Nelson-Aalen ===\n")
sf04 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        prediction.type = "Nelson-Aalen")
pred04 <- predict(sf04)$predictions
ft04 <- sf04$failure.times
K04 <- min(20, length(ft04))
cat("failure times:", length(ft04), "| using first", K04, "\n")

df04 <- as.data.frame(X)
df04$time <- Y; df04$status <- D
for (j in 1:K04) df04[[paste0("r_s", j)]] <- pred04[, j]
write_test(df04, "test04_predtype0")

## ============================================================
## TEST 05: predtype=1 (Kaplan-Meier, explicit)
## ============================================================
cat("\n=== TEST 05: predtype=1 KM explicit ===\n")
sf05 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        prediction.type = "Kaplan-Meier")
pred05 <- predict(sf05)$predictions
ft05 <- sf05$failure.times
K05 <- min(20, length(ft05))
cat("failure times:", length(ft05), "| using first", K05, "\n")

df05 <- as.data.frame(X)
df05$time <- Y; df05$status <- D
for (j in 1:K05) df05[[paste0("r_s", j)]] <- pred05[, j]
write_test(df05, "test05_predtype1")

## ============================================================
## TEST 06: nofastlogrank (fast.logrank=FALSE is the R default; test it explicitly)
## In Stata: default is fastlogrank (fast.logrank=TRUE), nofastlogrank turns it off
## In R grf 2.5.0: fast.logrank defaults to FALSE
## ============================================================
cat("\n=== TEST 06: nofastlogrank ===\n")
sf06 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        fast.logrank = FALSE)
pred06 <- predict(sf06)$predictions
ft06 <- sf06$failure.times
K06 <- min(20, length(ft06))
cat("failure times:", length(ft06), "| using first", K06, "\n")

df06 <- as.data.frame(X)
df06$time <- Y; df06$status <- D
for (j in 1:K06) df06[[paste0("r_s", j)]] <- pred06[, j]
write_test(df06, "test06_nofastlogrank")

## ============================================================
## TEST 07: cluster()
## ============================================================
cat("\n=== TEST 07: cluster() ===\n")
set.seed(42)
cluster_ids <- rep(1:50, each = 10)  # 50 clusters of 10
sf07 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        clusters = cluster_ids)
pred07 <- predict(sf07)$predictions
ft07 <- sf07$failure.times
K07 <- min(20, length(ft07))
cat("failure times:", length(ft07), "| using first", K07, "\n")

df07 <- as.data.frame(X)
df07$time <- Y; df07$status <- D; df07$cluster_id <- cluster_ids
for (j in 1:K07) df07[[paste0("r_s", j)]] <- pred07[, j]
write_test(df07, "test07_cluster")

## ============================================================
## TEST 08: weights()
## ============================================================
cat("\n=== TEST 08: weights() ===\n")
set.seed(42)
wts <- runif(n, 0.5, 2.0)
sf08 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        sample.weights = wts)
pred08 <- predict(sf08)$predictions
ft08 <- sf08$failure.times
K08 <- min(20, length(ft08))
cat("failure times:", length(ft08), "| using first", K08, "\n")

df08 <- as.data.frame(X)
df08$time <- Y; df08$status <- D; df08$wt <- wts
for (j in 1:K08) df08[[paste0("r_s", j)]] <- pred08[, j]
write_test(df08, "test08_weights")

## ============================================================
## TEST 09: nohonesty
## ============================================================
cat("\n=== TEST 09: nohonesty ===\n")
sf09 <- survival_forest(X, Y, D, num.trees = 500, seed = 42, honesty = FALSE)
pred09 <- predict(sf09)$predictions
ft09 <- sf09$failure.times
K09 <- min(20, length(ft09))
cat("failure times:", length(ft09), "| using first", K09, "\n")

df09 <- as.data.frame(X)
df09$time <- Y; df09$status <- D
for (j in 1:K09) df09[[paste0("r_s", j)]] <- pred09[, j]
write_test(df09, "test09_nohonesty")

## ============================================================
## TEST 10: mtry=2
## ============================================================
cat("\n=== TEST 10: mtry=2 ===\n")
sf10 <- survival_forest(X, Y, D, num.trees = 500, seed = 42, mtry = 2)
pred10 <- predict(sf10)$predictions
ft10 <- sf10$failure.times
K10 <- min(20, length(ft10))
cat("failure times:", length(ft10), "| using first", K10, "\n")

df10 <- as.data.frame(X)
df10$time <- Y; df10$status <- D
for (j in 1:K10) df10[[paste0("r_s", j)]] <- pred10[, j]
write_test(df10, "test10_mtry2")

## ============================================================
## TEST 11: minnodesize=20
## ============================================================
cat("\n=== TEST 11: minnodesize=20 ===\n")
sf11 <- survival_forest(X, Y, D, num.trees = 500, seed = 42, min.node.size = 20)
pred11 <- predict(sf11)$predictions
ft11 <- sf11$failure.times
K11 <- min(20, length(ft11))
cat("failure times:", length(ft11), "| using first", K11, "\n")

df11 <- as.data.frame(X)
df11$time <- Y; df11$status <- D
for (j in 1:K11) df11[[paste0("r_s", j)]] <- pred11[, j]
write_test(df11, "test11_minnodesize20")

## ============================================================
## TEST 12: Heavy censoring (80%)
## ============================================================
cat("\n=== TEST 12: Heavy censoring (~80%) ===\n")
set.seed(42)
C_heavy <- rexp(n, rate = 0.05)  # very light censoring rate -> heavy censoring
Y12 <- pmin(T_true, C_heavy)
D12 <- as.integer(T_true <= C_heavy)
cat("Events:", sum(D12), "/ censored:", sum(1-D12),
    "event_rate:", round(mean(D12), 3), "\n")

sf12 <- survival_forest(X, Y12, D12, num.trees = 500, seed = 42)
pred12 <- predict(sf12)$predictions
ft12 <- sf12$failure.times
K12 <- min(20, length(ft12))
cat("failure times:", length(ft12), "| using first", K12, "\n")

df12 <- as.data.frame(X)
df12$time <- Y12; df12$status <- D12
for (j in 1:K12) df12[[paste0("r_s", j)]] <- pred12[, j]
write_test(df12, "test12_heavy_censoring")

## ============================================================
## TEST 13: Light censoring (~10%)
## ============================================================
cat("\n=== TEST 13: Light censoring (~10%) ===\n")
set.seed(42)
C_light <- rexp(n, rate = 2.0)  # heavy censoring rate -> light censoring (most observe event)
Y13 <- pmin(T_true, C_light)
D13 <- as.integer(T_true <= C_light)
cat("Events:", sum(D13), "/ censored:", sum(1-D13),
    "event_rate:", round(mean(D13), 3), "\n")

sf13 <- survival_forest(X, Y13, D13, num.trees = 500, seed = 42)
pred13 <- predict(sf13)$predictions
ft13 <- sf13$failure.times
K13 <- min(20, length(ft13))
cat("failure times:", length(ft13), "| using first", K13, "\n")

df13 <- as.data.frame(X)
df13$time <- Y13; df13$status <- D13
for (j in 1:K13) df13[[paste0("r_s", j)]] <- pred13[, j]
write_test(df13, "test13_light_censoring")

## ============================================================
## TEST 14: Expected survival (compare R manual integration vs Stata)
## ============================================================
cat("\n=== TEST 14: Expected survival ===\n")
# Use sf01 (default) – predict with all failure times
sf14 <- survival_forest(X, Y, D, num.trees = 500, seed = 42)
pred14 <- predict(sf14)$predictions
ft14 <- sf14$failure.times
K14 <- length(ft14)

esurv14 <- expected_survival_r(sf14)
cat("E[T|X] range:", range(esurv14), "\n")
cat("E[T|X] mean:", mean(esurv14), "\n")

# Save ALL failure-time columns for integration comparison
df14 <- as.data.frame(X)
df14$time <- Y; df14$status <- D
for (j in 1:K14) df14[[paste0("r_s", j)]] <- pred14[, j]
df14$r_esurv <- esurv14
write_test(df14, "test14_expected_survival")
cat("Test 14: failure times=", K14, "\n")

## ============================================================
## TEST 15: Expected survival consistency check
## E[T|X] should be lower when X1 is large (higher hazard)
## ============================================================
cat("\n=== TEST 15: Expected survival consistency ===\n")
sf15 <- survival_forest(X, Y, D, num.trees = 500, seed = 42)
esurv15 <- expected_survival_r(sf15)

# Correlation with X1 should be negative (high X1 -> high hazard -> low survival)
cor_x1 <- cor(X[, 1], esurv15)
cat("Correlation E[T|X] with X1:", round(cor_x1, 4),
    "(should be negative)\n")

df15 <- as.data.frame(X)
df15$time <- Y; df15$status <- D
df15$r_esurv <- esurv15
write_test(df15, "test15_esurv_consistency")

## ============================================================
## TEST 16: numfailures(50) — R uses failure.times= (vector of specific times)
## Stata's numfailures(50) limits to 50 unique failure times.
## In R we subset to the 50 smallest unique event times.
## ============================================================
cat("\n=== TEST 16: numfailures(50) ===\n")
unique_ft <- sort(unique(Y[D == 1]))
ft_subset50 <- unique_ft[1:min(50, length(unique_ft))]
sf16 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        failure.times = ft_subset50)
pred16 <- predict(sf16)$predictions
ft16 <- sf16$failure.times
K16 <- min(20, length(ft16))
cat("failure times:", length(ft16), "| using first", K16, "\n")

df16 <- as.data.frame(X)
df16$time <- Y; df16$status <- D
for (j in 1:K16) df16[[paste0("r_s", j)]] <- pred16[, j]
write_test(df16, "test16_numfailures50")

## ============================================================
## TEST 17: Combined noutput=50 + predtype=0 (Nelson-Aalen) + cluster
## ============================================================
cat("\n=== TEST 17: Combined noutput50+NA+cluster ===\n")
set.seed(42)
cluster_ids17 <- rep(1:50, each = 10)
sf17 <- survival_forest(X, Y, D, num.trees = 500, seed = 42,
                        prediction.type = "Nelson-Aalen",
                        clusters = cluster_ids17)
pred17 <- predict(sf17)$predictions
ft17 <- sf17$failure.times
K17 <- min(50, length(ft17))
cat("failure times:", length(ft17), "| using first", K17, "\n")

df17 <- as.data.frame(X)
df17$time <- Y; df17$status <- D; df17$cluster_id <- cluster_ids17
for (j in 1:K17) df17[[paste0("r_s", j)]] <- pred17[, j]
write_test(df17, "test17_combined")

cat("\n=== ALL R TESTS COMPLETE ===\n")
