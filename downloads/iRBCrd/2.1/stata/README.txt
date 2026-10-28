iRBCrd: IMPROVED ROBUST BIAS CORRECTION FOR RD DESIGNS
Stata version 2.1 -- October 28, 2026
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
irbcrd provides bias-corrected point estimates, analytic confidence intervals,
and bootstrap confidence intervals for:
  - Sharp regression discontinuity (RD).
  - Fuzzy RD.
  - Covariate-adjusted sharp RD.
  - Covariate-adjusted fuzzy RD, adjusting both outcome and treatment equations.
  - Sharp and fuzzy kink RD, with or without covariate adjustment.

Use irbcrd for every design. Supply fuzzy() for a treatment variable and covs()
for covariate adjustment. Add kink for a kink design. The command selects the
appropriate implementation from these inputs; omitting kink selects level RD.

This is a distribution of native Stata/Mata source files. No license file is
included. Bootstrap inference uses a fixed-regressor wild bootstrap.

2. FILES AND LOADING
--------------------
irbcrd.mata       Local-polynomial fits, bias and variance calculations,
                  bandwidth regularization, and bootstrap implementations.
irbcrd.ado        Unified command interface, defaults, and results display.
README.txt        This manual.
irbcrd.sthlp      Stata help, accessed with help irbcrd.
examples.do      Reproducible level and kink examples.
CHANGELOG.txt    Release overview.
irbcrd.pkg       Local installation manifest.
stata.toc        Local installation catalog.
tests/           Standalone interface checks and R comparison scripts.
tests/VALIDATION.txt  Summary of the native regression checks.

Extract the ZIP into a folder. With Stata's working directory set to that folder:

  adopath ++ "."

Keep irbcrd.ado and irbcrd.mata together. Mata functions are compiled on first
use. No additional Stata packages are required for estimation. rdrobust is
optional for bandwidth selection or comparisons. R is not required.

Alternatively, install the extracted files on Stata's ado-path:

  net install irbcrd, from("/absolute/path/to/extracted/folder") replace

To run the examples or standalone checks after setting the working directory:

  do examples.do
  do tests/run_tests.do

These do-files replace the dataset in memory. The command requires Stata 16
or later and can be used interactively or in batch mode. This release was
tested in Stata/SE 18.0 on macOS. Start a new session after changing versions.

3. QUICK START
--------------
  clear
  set seed 17041998
  set obs 400
  generate double x = 2*runiform()-1
  generate double z = .3*x+rnormal()
  generate double d = rbinomial(1,.2+.6*(x>=0))
  generate double y = 1+x+x^2+.8*d+.5*z+rnormal()

  irbcrd y x, h(.6)
  irbcrd y x, h(.6) fuzzy(d)
  irbcrd y x, h(.6) covs(z)
  irbcrd y x, h(.6) fuzzy(d) covs(z)

The bandwidth above is illustrative. It is not an estimated optimal bandwidth.
The bootstrap and analytic reportable estimates are identical. Their intervals
can differ, and the bootstrap interval need not be symmetric around that
bias-corrected estimate.

4. COMMAND INTERFACE
---------------------
  irbcrd y x [if] [in], h(# [#])
      [c(#) p(#) kernel(name [name]) fuzzy(treatment) kink covs(varlist)
       res(name) alpha(#) bwcheck(#) fast
       boot_ci boot_b(#) boot_seed(#) boot_wild(name) boot_keep
       q_method(name) fuzzy_boot_weights(name) gamma_boot_weights(name)
       nodisplay]

Design and estimation:
  y, x       Numeric outcome and running-variable names. if/in restrict the
             sample; observations missing any analysis variable are omitted.
  h          Required positive bandwidth, one value or a pair ordered left, right.
  c          Cutoff; default 0. The right side includes x = c.
  p          Local-polynomial order; default 1 for levels and 2 for kinks.
             Use a nonnegative integer for levels, and p >= 1 for kinks.
  kernel     tri (default), epa, or uni. Long aliases triangular,
             epanechnikov, and uniform are accepted. A pair is allowed,
             ordered left, right.
  fuzzy      Treatment variable. Omission selects sharp RD. This is not a flag.
             Binary, discrete, and continuous treatments are accepted.
  kink       Add this flag for slope changes; omit it for level discontinuities.
  covs       Numeric covariate variable list; omitted by default.
             Both equations are adjusted in fuzzy RD.
  res        irbc-hc0 through irbc-hc3, or cct-hc0 through cct-hc3.
             The default is pointwise iRBC-HC3 (irbc-hc3).
  alpha      Significance level; default .05 gives a 95% confidence interval.
             This default does not depend on Stata's set level setting.
  q_method   sample or asymptotic. Default sample when fast is omitted;
             asymptotic when fast is supplied.
  fast       Standard-error shortcut; default off. See Section 8.
  bwcheck    Minimum-neighbor regularization rule; default 21. Zero disables it.
  nodisplay  Suppress the results display without changing calculations.

Bootstrap:
  boot_ci    Add this flag for a bootstrap interval; default off.
  boot_b     Number of draws; default 999 when boot_ci is supplied. A supplied
             value must be a positive integer. At least two finite draws are
             needed for an interval; use substantially more for inference.
  boot_seed  Optional numeric seed in [0,2147483647], converted to an integer.
             Omission uses the current Stata RNG state. A supplied seed is
             passed to set seed; the prior random-number state is not restored.
  boot_wild  rademacher (default), mammen, or normal.
  boot_keep  Retain bootstrap statistic draws and weight diagnostics; default off.
  fuzzy_boot_weights
             fixed (default) or bootstrap in every fuzzy design,
             including kink and covariate-adjusted designs.
  gamma_boot_weights
             fixed (default) or bootstrap in every covariate-adjusted
             design, including sharp, fuzzy, and kink designs.

Option names follow R's iRBCrd() interface with dots changed to underscores
and names written in lowercase. For example, boot.ci becomes boot_ci and
q.method becomes q_method. Stata flags replace R logical inputs: presence
means TRUE and omission means FALSE. Use bwcheck(0) for R's bwcheck = NULL.

5. BANDWIDTHS AND DATA
----------------------
The command does not select an optimal bandwidth. Supply h(), for example
from rdbwselect, with the appropriate fuzzy and covariate inputs in the
selector. For a kink, also set deriv(1) and the same p() in the selector.
A bandwidth pair is left, right; one value applies to both sides. Supply
bandwidths only through h(). Kernels follow the same convention through kernel().

  irbcrd y x, h(.5 .7) kernel(tri epa)

When using rdrobust for a comparison with rho = 1, explicitly use the same
selected h and set b = h there. For an exact bandwidth comparison, consider
bwcheck(0) here, provided all needed local fits remain nonsingular.

With bwcheck(k), a local bandwidth is enlarged if necessary to reach the kth
closest observation, or the farthest observation if fewer than k are present.
This counts observations, not distinct running-variable values. The initial
cutoff regularization uses pooled observations for a single h value,
and side-specific observations for paired bandwidth inputs.
Auxiliary fits can enlarge their own bandwidths as well.

For kink RD, cutoff regularization is side-specific even for a single h value.
The reported e(h_minus) and e(h_plus) are the actual slope-fitting bandwidths.
These effective bandwidths also start the side's auxiliary fits.
e(pilot_h_minus) and e(pilot_h_plus) record any further bandwidth enlargement
in the level pilots used for the bias. Kink bias, variance, and bootstrap
calculations use the same coefficient maps, including enlarged pilot fits.

Returned e(h) contains the requested minus/plus bandwidth pair.
e(h_minus) and e(h_plus) are output fields, not inputs; they report the effective
cutoff bandwidths, not the bandwidth at every auxiliary fitting point. This
rule does not guarantee full rank with ties, collinearity, or zero endpoint
kernel weights. Do not trim the data to the cutoff h window before calling
irbcrd: auxiliary fits need observations outside that window. Bias and
variance calculations can involve observations up to 2h from the cutoff when
regularization is disabled; pointwise residual fits use their own neighbors.

6. EXTENSIONS
-------------
All extensions use irbcrd and support analytic and bootstrap inference.
The design inputs can be combined; bandwidth and kernel conventions remain
as in Section 5.

Kink RD:
  Add kink to estimate the right slope minus the left slope at the cutoff.
  The default is local quadratic estimation (p = 2); p = 1 or a higher order
  can also be supplied. Sharp kink reports the slope difference in units
  of y/x, without dividing by a known policy-schedule slope change. For that
  target, rescale y by the known change before estimation. Bias correction
  uses slope-specific Q adjustments; residuals and bootstrap outcomes
  remain in levels.

    irbcrd y x, h(.6) kink

Fuzzy RD:
  Supply the treatment variable through fuzzy(). The target is the outcome
  discontinuity divided by the treatment discontinuity, or the corresponding
  ratio of slope changes when kink is supplied. iRBC corrects the linearized
  ratio and retains the covariance between the outcome and treatment
  components. A near-zero first stage prevents ordinary ratio inference;
  a nonzero estimate alone does not establish a strong first stage.
  Bootstrap weighting choices are described in Section 7.

    irbcrd y x, h(.6) fuzzy(d)

Covariate adjustment:
  Supply covs() to add common covariate slopes to the side-specific order-p
  weighted regressions. The outcome discontinuity is adjusted by the
  estimated covariate coefficients times the covariate discontinuities.
  Fuzzy designs adjust both outcome and treatment equations, using separate
  coefficient vectors and retaining their covariance. In kink designs, the
  adjustment uses covariate slope changes. Preserving the unadjusted target
  requires covariate balance at the cutoff: continuity for level RD and no
  slope change for kink RD, alongside the design's identification conditions.

    irbcrd y x, h(.6) covs(z)

7. BOOTSTRAP
-------------
The running variable, bandwidths, and Q adjustments are held fixed. For local
observations, bootstrap outcomes equal pointwise local-polynomial fitted
values plus HC-adjusted residuals multiplied by wild weights. The same draw
is used for all components of an observation in joint designs, preserving
their dependence. Observations outside the cutoff working window are left
unchanged; the bandwidth is not re-estimated in each bootstrap draw.

  irbcrd y x, h(.6) boot_ci boot_b(999) boot_seed(17041998)

For each sharp component j, the bootstrap statistic is
  T_j* = Q_plus [ghat_j,plus*(c) - ghat_j,plus(c)]
       - Q_minus [ghat_j,minus*(c) - ghat_j,minus(c)].

For kink designs, replace the cutoff level fits in this expression with
cutoff slope fits. Bootstrap outcomes and residuals are still constructed
in levels, and Q uses the slope-specific adjustment.

The sharp, fuzzy, and covariate-adjusted statistics use the relevant component
linear combination. By default, fuzzy and covariate weights remain fixed.
The two choices can be changed independently in level and kink designs:
  fuzzy_boot_weights(bootstrap) uses the bootstrap first-stage discontinuity
  (or slope change) in the denominator of the fuzzy weights. The original
  outcome-to-treatment ratio remains fixed, as in the linearized statistic.
  gamma_boot_weights(bootstrap) refits the common covariate coefficients
  using bootstrap outcomes and covariates. Fuzzy designs also refit the
  treatment equation. The bootstrap first stage is adjusted using whichever
  covariate coefficients are selected by gamma_boot_weights().
Thus, covariate-adjusted fuzzy designs allow all four fixed/bootstrap
combinations. These options change the bootstrap statistic's combination
weights, not the original analytic estimate, variance, or Q adjustments.
Draws with singular covariate refits or undefined fuzzy weights are omitted
from the interval calculation.

Let se_us denote the standard error of the uncorrected Q-weighted component
combination, and mhat = se_irbc/se_us. The interval uses empirical quantiles
of the bootstrap statistic at
  normal(mhat * invnormal(alpha/2))
  normal(mhat * invnormal(1-alpha/2)).
Quantiles use linear interpolation equivalent to R's type = 7 and only finite
draws. They are inverted around the corresponding uncorrected discontinuity
or ratio estimate. The REPORTABLE bootstrap point estimate is nevertheless
the same bias-corrected estimate as for analytic iRBC. This distinction
is intentional.

Rademacher weights take -1 or 1 with equal probability. Mammen weights use
the usual mean-zero, variance-one asymmetric two-point distribution. Normal
weights are also available. Bootstrap quantiles have simulation error,
especially at high confidence levels. B = 999 is a default, not a guarantee
of accurate tail quantiles. Check stability across seeds or larger B when
differences between intervals are small. R and Stata use different RNGs;
an identical numeric seed does not produce identical draws across languages.

8. OTHER DESIGN CHOICES
-----------------------
Residuals:
  iRBC-HCk fits a p-th order local polynomial separately on each cutoff side,
  centered at each X_i. The raw residual is Y_i - ghat(X_i). There is no
  additional Q-weighted bias subtraction in this fitted value. The HC
  adjustment uses the pointwise fit's leverage and local sample size.

  CCT-HCk fits a (p+1)-th order local polynomial at the cutoff on each side.
  Evaluate it at X_i, subtract it from Y_i, and apply the HC adjustment using
  that cutoff fit's leverage and local sample size. Both residual choices
  use fits in levels, including for kink designs.

For raw residual e_i and leverage ell_i:
  HC0: e_i
  HC1: e_i / sqrt((n_local - df) / n_local)
  HC2: e_i / sqrt(1 - ell_i)
  HC3: e_i / (1 - ell_i)

Here df = p+1 for iRBC residuals and df = p+2 for CCT residuals; n_local counts
nonzero-kernel observations. The same residuals are used in the variance
estimator and bootstrap data-generating process. In joint designs, residuals
are evaluated componentwise and their cross-products enter the covariance.

hc0 through hc3 select iRBC-HCk; rbc-hc0 through rbc-hc3 select CCT-HCk.
The returned e(res) field uses irbc-hck or cct-hck, respectively, with k
replaced by 0, 1, 2, or 3.

res(loo) uses unscaled residuals from the ordinary pointwise auxiliary fits,
without an HC adjustment.

Q and analytic inference:
  q_method(sample) uses finite-sample Q_n values (the default).
  q_method(asymptotic) uses kernel-specific population values; omit fast
  to retain the full variance calculation. Asymptotic Q requires p = 1 for
  level RD and also supports higher orders for kink RD.
  fast uses asymptotic Q and a standard-error shortcut for sharp designs,
  with slope-specific constants for kink RD; joint designs retain full
  covariance matrices. It cannot use sample Q.

The analytic interval is estimate_irbc +/- invnormal(1-alpha/2) * se_irbc.

9. RESULTS
-----------
The primary irbcrd result fields are:
  e(estimate_irbc)       Bias-corrected point estimate.
  e(se_irbc)             Analytic standard error.
  e(ci_irbc)             Analytic interval, lower then upper.
  e(estimate_irbc_boot)  Same point estimate, when boot_ci is supplied.
  e(ci_irbc_boot)        Bootstrap interval, lower then upper.
  e(h)                  Requested bandwidths, as a minus/plus row vector.
  e(h_minus), e(h_plus)  Effective cutoff bandwidths.
  e(kernel), e(kernel_minus), e(kernel_plus)
                        Kernel specifications used on the two sides.
  e(res)                Residual specification.
  e(boot_B), e(boot_wild), e(boot_t)
                        Bootstrap details; e(boot_t) is omitted by default.

Supply boot_keep to retain bootstrap statistics and, where applicable,
component statistics, covariate coefficients, and combination weights for
individual draws. This does not change point estimates or confidence intervals.
The retained arrays include e(boot_t_components), e(boot_linear_weights),
e(boot_fuzzy_weights), e(boot_adjusted_first_stage), and e(boot_gamma) or
e(boot_gamma_y)/e(boot_gamma_d), depending on the design.

Design-specific details include uncorrected discontinuities, bias components,
first-stage estimates, covariate slopes, covariance matrices, and component
estimates. Sharp designs report e(tau_hat); fuzzy designs report e(zeta_hat)
and e(first_stage). Adjusted designs report e(gamma_hat), or e(gamma_y_hat)
and e(gamma_d_hat) when fuzzy. Component estimates, corrections, and covariance
matrices are returned as labeled matrices. e(epshat_minus) and e(epshat_plus)
are computational inputs over their local support, not full-sample diagnostics
for fitted values.

Kink results also report e(kink), e(p), e(design), e(linear_weights),
e(kink_weights_minus), e(kink_weights_plus), e(pilot_h_minus), and
e(pilot_h_plus). Component estimates are slope changes. Coefficient maps and
residuals follow the selected observation order within each side.

e(b), e(V), and e(sample) follow Stata's estimation-command conventions.
Excluded observations are counted in e(N_omitted). Large returned matrices
can be accessed using Mata's st_matrix().

10. REFERENCES
---------------
Cavaliere, G., Goncalves, S., Nielsen, M. O., and Zanelli, E. (2026):
  Improved Inference for Nonparametric Regression and Regression-Discontinuity
  Designs. arXiv:2512.00566. https://arxiv.org/abs/2512.00566.
Calonico, S., Cattaneo, M. D., Farrell, M. H., and Titiunik, R. (2019):
  Regression Discontinuity Designs Using Covariates, Review of Economics
  and Statistics, 101(3), 442-451.

The code is supplied as native Stata/Mata source files. See CHANGELOG.txt for
the release overview and tests/VALIDATION.txt for checks performed before
preparing the ZIP.
