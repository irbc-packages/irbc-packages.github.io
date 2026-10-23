iRBCrd: IMPROVED ROBUST BIAS CORRECTION FOR RD DESIGNS
Version 2.0 -- October 23, 2026
Website: https://irbc-packages.github.io/

AUTHORS AND CONTACTS
--------------------
Giuseppe Cavaliere
  University of Bologna and Exeter Business School
  giuseppe.cavaliere@unibo.it
Silvia Goncalves
  McGill University
  silvia.goncalves@mcgill.ca
Morten Orregaard Nielsen
  Aarhus University
  mon@econ.au.dk
Edoardo Zanelli
  Aarhus University
  ez@econ.au.dk

1. OVERVIEW
-----------
iRBCrd provides bias-corrected point estimates, analytic confidence intervals,
and plug-in wild bootstrap confidence intervals for:
  - Sharp regression discontinuity (RD).
  - Fuzzy RD.
  - Covariate-adjusted sharp RD.
  - Covariate-adjusted fuzzy RD, adjusting both outcome and treatment equations.

Use iRBCrd() for every design. Supply fuzzy for a treatment variable and covs
for covariate adjustment. The function selects the appropriate implementation
from these inputs.

This is a distribution of R functions loaded with source(), not an installable
R package or a CRAN release. No license file is included. Bootstrap inference
uses the plug-in wild bootstrap.

2. FILES AND LOADING
--------------------
iRBC_functions.R  Kernels, local-polynomial fits, variance calculations, and
                  bandwidth regularization helpers.
iRBCrd.R          Unified entry point and the four RD implementations.
README.txt        This manual.
examples.R        Reproducible examples for the four designs.
CHANGELOG.txt     Release overview.
smoke_test.R      Standalone interface and bootstrap checks.
VERIFICATION.txt  Summary of the development regression checks.
CONTACT_SOURCES.txt  Sources used to verify the institutional contacts.

Extract the ZIP into a folder. With R's working directory set to that folder:

  source("iRBC_functions.R")
  source("iRBCrd.R")

Both files are needed, in that order. No additional R packages are required
for estimation. rdrobust is optional for bandwidth selection or comparisons.
Loading these files does not attach or install any additional packages.

To run the examples or standalone checks after setting the working directory:

  source("examples.R")
  source("smoke_test.R")

The functions can be used in R, RStudio, or via Rscript. They do not require
RStudio or a compiler. This release was tested in R 4.5.2 on macOS.

3. QUICK START
--------------
  set.seed(17041998)
  n <- 400
  x <- runif(n, -1, 1)
  z <- 0.3 * x + rnorm(n)
  d <- rbinom(n, 1, 0.2 + 0.6 * (x >= 0))
  y <- 1 + x + x^2 + 0.8 * d + 0.5 * z + rnorm(n)
  Z <- cbind(z = z)

  sharp <- iRBCrd(y, x, h = 0.6)
  fuzzy_fit <- iRBCrd(y, x, h = 0.6, fuzzy = d)
  adjusted_sharp <- iRBCrd(y, x, h = 0.6, covs = Z)
  adjusted_fuzzy <- iRBCrd(y, x, h = 0.6, fuzzy = d, covs = Z)

  fit <- iRBCrd(y, x, h = 0.6, fuzzy = d, covs = Z,
                res = "irbc-hc3", q.method = "sample", fast = FALSE,
                boot.ci = TRUE, boot.B = 999,
                boot.wild = "rademacher", boot.seed = 17041998)
  fit$estimate_irbc
  fit$se_irbc
  fit$ci_irbc
  fit$estimate_irbc_boot
  fit$ci_irbc_boot
  fit$boot.success

The bandwidth above is illustrative. It is not an estimated optimal bandwidth.
The bootstrap and analytic reportable estimates are identical. Their intervals
can differ, and the bootstrap interval need not be symmetric around that
bias-corrected estimate.

4. FUNCTION INTERFACE
---------------------
  iRBCrd(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
         alpha = NULL, fast = NULL, g.loo = NULL,
         boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
         boot.wild = NULL, boot.keep = NULL,
         bwcheck = .bwcheck_default,
         h.plus = NULL, h.minus = NULL,
         kersel.plus = NULL, kersel.minus = NULL,
         q.method = NULL,
         fuzzy = NULL, covs = NULL,
         fuzzy.boot.weights = NULL, gamma.boot.weights = NULL,
         d = NULL)

Design and estimation:
  y, x       Numeric vectors in the same observation order. Prepare finite,
             complete-case data before calling the function.
  h          Required positive bandwidth, scalar or a pair ordered left, right.
  c          Cutoff; default 0. The right side includes x = c.
  p          Local-polynomial order; default 1. Use a nonnegative integer.
  fuzzy      Treatment vector. NULL selects sharp RD. This is not a Boolean.
             Binary, discrete, and continuous treatments are accepted.
  d          Alias for fuzzy. Do not supply both fuzzy and d.
  covs       Numeric covariate vector, matrix, or data frame; default NULL.
             One row per observation. Both equations are adjusted in fuzzy RD.
  kersel     "tri" (default), "epa", or "uni". Long aliases "triangular",
             "epanechnikov", and "uniform" are accepted. A pair is allowed.
  res        "irbc-hc0" through "irbc-hc3", or "cct-hc0" through "cct-hc3".
             The default is pointwise iRBC-HC3 ("irbc-hc3").
  alpha      Significance level; default 0.05 gives a 95% confidence interval.
  q.method   "sample" or "asymptotic". Default "sample" when fast = FALSE;
             "asymptotic" when fast = TRUE.
  fast       Standard-error shortcut; default FALSE. See Section 7.
  g.loo      Leave-one-out auxiliary fits for bias estimation; default FALSE.
  bwcheck    Minimum-neighbor regularization rule; default 21. NULL disables it.
  h.minus, h.plus
             Optional left/right bandwidth overrides.
  kersel.minus, kersel.plus
             Optional left/right kernel overrides.

Bootstrap:
  boot.ci    TRUE to add a plug-in bootstrap interval; default FALSE.
  boot.B     Number of draws; default 999 when boot.ci = TRUE. A supplied
             value must be a positive integer. At least two finite draws are
             needed for an interval; use substantially more for inference.
  boot.seed  Optional numeric seed. Default NULL uses the current R RNG state.
             A supplied seed is passed to set.seed() after integer conversion.
             The prior random-number state is not restored.
  boot.wild  "rademacher" (default), "mammen", or "normal".
  boot.keep  Retain bootstrap statistic draws; default TRUE.
  fuzzy.boot.weights
             Plain fuzzy RD: "fixed" (default) or "bootstrap".
             Covariate-adjusted fuzzy RD: "fixed" only.
  gamma.boot.weights
             Covariate-adjusted fuzzy RD: "fixed" only (default).

5. BANDWIDTHS AND DATA
----------------------
The function does not select an optimal bandwidth. Supply h, for example from
rdrobust::rdbwselect(), with the appropriate fuzzy and covariate inputs in the
selector. An unnamed bandwidth pair is left, right. Named pairs accept
left/right or minus/plus. A scalar applies to both sides. Explicit side
arguments override those choices. Kernels follow the same convention.

  fit <- iRBCrd(y, x, h = c(left = 0.5, right = 0.7),
                kersel = c(left = "tri", right = "epa"))

When using rdrobust for a comparison with rho = 1, explicitly use the same
selected h and set b = h there. For an exact bandwidth comparison, consider
bwcheck = NULL here, provided all needed local fits remain nonsingular.

With bwcheck = k, a local bandwidth is enlarged if necessary to reach the kth
closest observation, or the farthest observation if fewer than k are present.
This counts observations, not distinct running-variable values. The initial
cutoff regularization uses pooled observations for scalar h without side
overrides, and side-specific observations for paired bandwidth inputs.
Auxiliary fits can enlarge their own bandwidths as well.

Returned h.requested, h.minus, and h.plus describe the requested and effective
cutoff bandwidths, not the bandwidth at every auxiliary fitting point. This
rule does not guarantee full rank with ties, collinearity, or zero endpoint
kernel weights. Do not trim the data to the cutoff h window before calling
iRBCrd(): auxiliary fits need observations outside that window. Bias and
variance calculations can involve observations up to 2h from the cutoff when
regularization is disabled; pointwise residual fits use their own neighbors.

6. RESIDUALS
------------
iRBC-HCk:
  Fit a p-th order local polynomial separately on each cutoff side, centered
  at each X_i. The raw residual is Y_i - ghat(X_i). There is no additional
  Q-weighted bias subtraction in this fitted value. Apply the HC adjustment
  using the corresponding pointwise fit's leverage and local sample size.

CCT-HCk:
  Fit a (p+1)-th order local polynomial at the cutoff on each side. Evaluate
  it at X_i, subtract it from Y_i, and apply the HC adjustment using that
  cutoff fit's leverage and local sample size.

For raw residual e_i and leverage ell_i:
  HC0: e_i
  HC1: e_i / sqrt((n_local - df) / n_local)
  HC2: e_i / sqrt(1 - ell_i)
  HC3: e_i / (1 - ell_i)

Here df = p+1 for iRBC residuals and df = p+2 for CCT residuals; n_local counts
nonzero-kernel observations. The same residuals are used in the variance
estimator and bootstrap data-generating process. In joint designs, residuals
are evaluated componentwise and their cross-products enter the covariance.

"hc0" through "hc3" select iRBC-HCk; "rbc-hc0" through "rbc-hc3" select
CCT-HCk. The returned res field uses "irbc-hck" or "cct-hck", respectively,
with k replaced by 0, 1, 2, or 3.

res = "loo" uses unscaled residuals from the auxiliary fits. g.loo = TRUE
also changes the bias correction's auxiliary fits. res = "loo" by itself
does not activate leave-one-out fitting, and g.loo does not switch the separate
iRBC-HCk residual construction to leave-one-out fitting.

7. Q AND ANALYTIC INFERENCE
---------------------------
q.method = "sample", fast = FALSE uses the finite-sample Q_n values.
q.method = "asymptotic", fast = FALSE substitutes kernel-specific Q constants
in the bias, variance, and bootstrap calculations, retaining the full variance
calculation. Asymptotic Q is available only for p = 1:
  triangular: 1.4082; Epanechnikov: 1.3571; uniform: 1.2000.

fast = TRUE also requires p = 1 and asymptotic Q. In sharp RD it approximates
the side-specific standard errors by 0.84, 0.83, or 0.86 times the corresponding
RBC standard errors for the three kernels. Joint designs instead use their
full component covariance matrix for the final standard error. fast = TRUE
cannot be combined with q.method = "sample".

Let tau_y, tau_d, and tau_z denote uncorrected discontinuity estimates, and
B_y, B_d, and B_z their iRBC bias corrections on the estimator scale.
  Sharp RD: estimate = tau_y - B_y.
  Fuzzy RD: estimate = tau_y/tau_d - B_y/tau_d + tau_y*B_d/tau_d^2.
The fuzzy correction is for the linearized ratio; it is not a ratio of two
separately bias-corrected estimates.

Covariate adjustment estimates common slopes while allowing separate p-th
order polynomial terms on each side. For sharp RD it adjusts tau_y and B_y
by gamma' tau_z and gamma' B_z. Fuzzy RD adjusts both outcome and treatment
discontinuities and their bias corrections before applying the linearized
ratio correction. Covariances between components are retained.

The analytic interval is the bias-corrected estimate plus/minus
qnorm(1-alpha/2) * se_irbc.

8. PLUG-IN BOOTSTRAP
--------------------
The running variable, bandwidths, and Q adjustments are held fixed. For local
observations, bootstrap outcomes equal pointwise local-polynomial fitted
values plus HC-adjusted residuals multiplied by wild weights. The same draw
is used for all components of an observation in joint designs, preserving
their dependence. Observations outside the cutoff working window are left
unchanged; the bandwidth is not re-estimated in each bootstrap draw.

For each sharp component j, the bootstrap statistic is
  T_j* = Q_plus [ghat_j,plus*(c) - ghat_j,plus(c)]
       - Q_minus [ghat_j,minus*(c) - ghat_j,minus(c)].

The sharp, fuzzy, and covariate-adjusted statistics use the relevant component
linear combination. By default, fuzzy and covariate weights remain fixed.
Plain fuzzy RD additionally allows bootstrapped fuzzy weights as specified
by fuzzy.boot.weights. Covariate-adjusted designs keep gamma fixed; adjusted
fuzzy RD also keeps its fuzzy linearization weights fixed.

Let se_us denote the standard error of the uncorrected Q-weighted component
combination, and mhat = se_irbc/se_us. The interval uses empirical quantiles
of the bootstrap statistic at
  pnorm(mhat * qnorm(c(alpha/2, 1-alpha/2))).
Quantiles use R's type = 7 interpolation and only finite draws. They are
inverted around the corresponding uncorrected discontinuity or ratio estimate.
The REPORTABLE bootstrap point estimate is nevertheless the same bias-corrected
estimate as for analytic iRBC. This distinction is intentional.

Rademacher weights take -1 or 1 with equal probability. Mammen weights use
the usual mean-zero, variance-one asymmetric two-point distribution. Normal
weights are also available. Bootstrap quantiles have simulation error,
especially at high confidence levels. B = 999 is a default, not a guarantee
of accurate tail quantiles. Inspect boot.success and stability across seeds
or larger B when differences between intervals are small.

9. RESULTS
-----------
The primary iRBCrd() result fields are:
  estimate_irbc       Bias-corrected point estimate.
  se_irbc             Analytic standard error.
  ci_irbc             Analytic interval, lower then upper.
  estimate_irbc_boot  Same point estimate, when boot.ci = TRUE.
  ci_irbc_boot        Plug-in bootstrap interval, lower then upper.
  h.requested        Requested bandwidths.
  h.minus, h.plus     Effective cutoff bandwidths.
  Q.minus, Q.plus     Side-specific Q values.
  res, q.method      Residual and Q choices.
  boot.B, boot.success, boot.wild, boot.probs, boot.mhat, boot.t
                     Bootstrap details; boot.t requires boot.keep = TRUE.

Design-specific details include uncorrected discontinuities, bias components,
first-stage estimates, covariate slopes, covariance matrices, and sharp.*
component objects. Residual arrays are computational inputs over their local
support, not a full-sample diagnostic for fitted values.

10. REFERENCES
---------------
Cavaliere, G., Goncalves, S., Nielsen, M. O., and Zanelli, E.:
  Improved Inference for Regression Discontinuity Designs (current manuscript).
  Earlier public version: Improved Inference for Nonparametric Regression
  and Regression-Discontinuity Designs, https://arxiv.org/abs/2512.00566.
Calonico, S., Cattaneo, M. D., Farrell, M. H., and Titiunik, R. (2019):
  Regression Discontinuity Designs Using Covariates, Review of Economics
  and Statistics, 101(3), 442-451.

The code is supplied as R functions. No install.packages("iRBCrd") or
library(iRBCrd) step applies to this distribution. See CHANGELOG.txt for the
release overview and VERIFICATION.txt for tests performed before preparing
the ZIP.
