# Set the working directory to the extracted release folder before sourcing.
source("iRBC_functions.R")
source("iRBCrd.R")

set.seed(17041998)
n <- 400
x <- runif(n, -1, 1)
z <- 0.3 * x + rnorm(n)
d <- rbinom(n, 1, 0.2 + 0.6 * (x >= 0))
y <- 1 + x + x^2 + 0.8 * d + 0.5 * z + rnorm(n)
Z <- cbind(z = z)

# The bandwidth is illustrative, not an estimated optimum.
fits <- list(
  sharp = iRBCrd(y, x, h = 0.6),
  fuzzy = iRBCrd(y, x, h = 0.6, fuzzy = d),
  adjusted_sharp = iRBCrd(y, x, h = 0.6, covs = Z),
  adjusted_fuzzy = iRBCrd(y, x, h = 0.6, fuzzy = d, covs = Z)
)
results <- do.call(rbind, lapply(fits, function(f) {
  c(estimate = f$estimate_irbc, se = f$se_irbc,
    lower = f$ci_irbc[1], upper = f$ci_irbc[2])
}))
print(results)

# B = 99 keeps this example quick; use a larger B for substantive inference.
boot_fit <- iRBCrd(y, x, h = 0.6, fuzzy = d, covs = Z,
  res = "irbc-hc3", q.method = "sample", fast = FALSE,
  boot.ci = TRUE, boot.B = 99, boot.wild = "rademacher",
  boot.seed = 17041998, boot.keep = TRUE)
print(rbind(analytic = boot_fit$ci_irbc, bootstrap = boot_fit$ci_irbc_boot))
stopifnot(identical(boot_fit$estimate_irbc, boot_fit$estimate_irbc_boot),
          boot_fit$boot.success == 99)

# Optional externally estimated bandwidth, requiring rdrobust:
# install.packages("rdrobust")
# bw <- rdrobust::rdbwselect(y, x, fuzzy = d, covs = Z,
#   p = 1, q = 2, kernel = "triangular", vce = "hc3",
#   bwselect = "mserd")
# h <- as.numeric(bw$bws[1, 1:2])
# fit <- iRBCrd(y, x, h = h, fuzzy = d, covs = Z, bwcheck = NULL)
# comparison <- rdrobust::rdrobust(y, x, fuzzy = d, covs = Z,
#   p = 1, q = 2, kernel = "triangular", vce = "hc3", h = h, b = h)
# Setting b = h in this comparison imposes rho = 1.

# Kink designs: continuous outcome and treatment curves with slope changes.
set.seed(17041998)
x <- runif(500, -1, 1)
Z <- cbind(z = .3 * x + rnorm(length(x)))
d <- .3 * x + 2 * pmax(x, 0) + .2 * Z[, 1] + rnorm(length(x), sd = .1)
y <- 1 + x + x^2 + .8 * d + .4 * Z[, 1] + rnorm(length(x), sd = .15)
kink_fits <- list(
  sharp = iRBCrd(y, x, h = .7, deriv = 1),
  fuzzy = iRBCrd(y, x, h = .7, deriv = 1, fuzzy = d),
  adjusted_sharp = iRBCrd(y, x, h = .7, deriv = 1, covs = Z),
  adjusted_fuzzy = iRBCrd(y, x, h = .7, deriv = 1, fuzzy = d, covs = Z,
    q.method = "sample", res = "cct-hc3", boot.ci = TRUE, boot.B = 99,
    boot.wild = "rademacher", boot.seed = 17041998)
)
print(do.call(rbind, lapply(kink_fits, function(f) {
  c(estimate = f$estimate_irbc, se = f$se_irbc,
    lower = f$ci_irbc[1], upper = f$ci_irbc[2])
})))
print(kink_fits$adjusted_fuzzy$ci_irbc_boot)

# Population kink Q, retaining the full sandwich variance:
asymptotic_kink <- iRBCrd(y, x, h = .7, deriv = 1, p = 2,
                          q.method = "asymptotic", fast = FALSE)
print(c(Q_minus = asymptotic_kink$Q.minus, Q_plus = asymptotic_kink$Q.plus))

# Optional estimated kink bandwidths: use the same deriv, p, fuzzy, and covs
# settings in the selector and both estimators. Use "cerrd" for CE selection.
# bw <- rdrobust::rdbwselect(y, x, deriv = 1, p = 2, q = 3,
#   fuzzy = d, covs = Z, kernel = "triangular", vce = "hc3", bwselect = "mserd")
# h <- as.numeric(bw$bws[1, 1:2])
# kink <- iRBCrd(y, x, h = h, deriv = 1, p = 2, fuzzy = d, covs = Z,
#   res = "cct-hc3", bwcheck = NULL, boot.ci = TRUE, boot.B = 999,
#   boot.seed = 17041998)
# rbc <- rdrobust::rdrobust(y, x, deriv = 1, p = 2, q = 3,
#   fuzzy = d, covs = Z, h = h, b = h, vce = "hc3", kernel = "triangular")
