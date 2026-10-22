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
