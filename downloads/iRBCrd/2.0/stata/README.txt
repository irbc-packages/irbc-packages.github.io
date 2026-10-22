iRBCrd for Stata / Mata
Version 2.0 DEVELOPMENT BUILD -- October 22, 2026

IMPORTANT: THIS IS NOT A VALIDATED PUBLIC RELEASE
------------------------------------------------
The native Stata/Mata source has been written against the R 2.0 implementation.
No licensed Stata installation was available during development. The code has
NOT been compiled or executed in Stata. Do not rely on it for research results
or publish it as validated until tests/run_tests.do passes in Stata.

Target: Stata 16 or later. This is a compatibility target, not a tested claim.
The estimation command has no R, Python, compiler, or rdrobust dependency.
R is used only by the optional developer tools that generate test fixtures.
This build is separate from the R release. The R implementation and existing
empirical results have not been changed.

AUTHORS
-------
Giuseppe Cavaliere, Silvia Goncalves, Morten Orregaard Nielsen, and
Edoardo Zanelli.
Contacts: giuseppe.cavaliere@unibo.it; silvia.goncalves@mcgill.ca;
          mon@econ.au.dk; ez@econ.au.dk.
No license has been added to this development distribution.

FILES
-----
irbcrd.ado       Unified Stata command, option parsing, and stored results.
irbcrd.mata      Native numerical implementation; compiled on first use.
irbcrd.sthlp     Stata help, accessible with help irbcrd.
examples.do     Illustrative examples for the four designs.
irbcrd.pkg      Local net-install manifest.
stata.toc       Local installation catalog.
tests/          R reference snapshots, numeric fixtures, formula audit, and
                the native Stata validation script.
VALIDATION.txt  What has and has not been tested.

LOCAL USE
---------
Set Stata's working directory to this folder, then type:

  adopath ++ "."
  help irbcrd
  do examples.do

For installation from this folder rather than a temporary search path:

  net install irbcrd, from("/absolute/path/to/this/folder")

The website download is a ZIP archive; use its extracted folder for local
net install. Keep irbcrd.ado and irbcrd.mata together on the ado-path.

QUICK START
-----------
  irbcrd y x, h(.6)
  irbcrd y x, h(.6) fuzzy(d)
  irbcrd y x, h(.6) covs(z1 z2)
  irbcrd y x, h(.6) fuzzy(d) covs(z1 z2)

  irbcrd y x, h(.5 .7) fuzzy(d) covs(z1 z2) ///
      residuals(cct-hc3) qmethod(sample) bootstrap ///
      reps(999) wild(rademacher) seed(17041998)

The bandwidths in these examples are illustrative, not estimated optima.
The command does not estimate a bandwidth: supply one selected externally,
for example with rdbwselect and the same fuzzy/covariate design inputs.

OPTIONS
-------
All designs use irbcrd y x [if] [in], h(...).

h(# [#])              Required bandwidth(s), ordered left, right.
c(#)                  Cutoff, default 0; x=c belongs to the right side.
p(#)                  Polynomial order, default 1, nonnegative integer.
fuzzy(d)              Treatment variable; omitting it selects sharp RD.
covs(z1 z2 ...)        Numeric covariates, with a common slope for each
                      component and side-specific polynomial terms.
                      Fuzzy RD adjusts both outcome and treatment equations.
kernel(tri [epa])      One or two kernels, left then right. Choices: tri,
                      epa, uni; triangular, epanechnikov, uniform accepted.
hleft(#), hright(#)    Override the respective h() entry.
kleft(...), kright(...)  Override the respective kernel entry.
residuals(...)        irbc-hc0..irbc-hc3 or cct-hc0..cct-hc3; default irbc-hc3.
                      hc0..hc3 aliases select iRBC residuals; rbc-hc0..rbc-hc3
                      aliases select CCT residuals. loo is also available.
qmethod(sample)       Finite-sample Qn, default unless fast is requested.
qmethod(asymptotic)   Kernel-specific Q: tri 1.4082, epa 1.3571, uni 1.2.
                      Requires p(1); does not request fast standard errors.
fast                 Asymptotic Q and the R fast-SE convention. In plain
                      sharp RD use .84/.83/.86 times side-specific RBC SEs;
                      joint designs retain full component covariances.
level(#)             Confidence level, default Stata's c(level), normally 95.
bwcheck(#)           Minimum-observation bandwidth regularization; default 21.
noregularize         Disable regularization; equivalent to bwcheck(0).
gloo                 Leave-one-out auxiliary bias fits; off by default.
bootstrap            Add the plug-in wild bootstrap CI; analytic CI remains.
reps(#)              Bootstrap draws, default 999. At least 2 finite draws
                      are required for a CI. More are needed for reliable tails.
wild(...)            rademacher (default), mammen, or normal; rad and gaussian
                      aliases are accepted.
seed(#)              Optional Stata seed. Changes the current RNG state, as
                      R's boot.seed does. Requires bootstrap.
fuzzyweights(fixed)  Default fuzzy linearization weights.
fuzzyweights(bootstrap)  Available for fuzzy RD WITHOUT covariates only,
                      matching R 2.0's bootstrapped first-stage denominator.
gammaweights(fixed)  Only available with fuzzy() and covs(); fixed is the only
                      supported choice. Covariate slopes are never resampled.
nokeep               Do not retain individual bootstrap statistic draws.
nodisplay            Suppress the results display, but still return estimates.
drawmatrix(name)     Advanced reproducibility/testing input. A Stata matrix
                      with N rows (complete selected sample in its original
                      order) and reps columns supplies the wild multipliers.
                      Only local rows are used. Requires bootstrap; do not
                      combine with seed(). Reported wild distribution is
                      'supplied'. This option is subject to Stata matrix limits.

NUMERICAL CONVENTIONS
---------------------
The analytic estimate is bias corrected. Fuzzy correction uses the linearized
ratio, not a ratio of separately corrected estimates. Covariate coefficients
are estimated using pooled common-slope weighted regressions after removing
side-specific polynomial terms. Joint covariances are preserved.

iRBC-HCk uses Y_i minus the order-p local fit centered at X_i. CCT-HCk uses
the order-(p+1) cutoff fit evaluated at X_i. HC adjustments use the matching
local leverage and positive-kernel sample size. Residuals are used in both
variance estimation and the bootstrap DGP.

Bootstrap X, h, and Q are fixed. An observation receives the same wild weight
in every component (outcome, treatment, and covariates). Observations outside
the R routine's working window remain unchanged. The statistic is the
appropriate linear combination of Q-weighted bootstrap component changes.
Empirical quantiles use R type-7 interpolation at
  normal((se_irbc/se_us) * invnormal((alpha/2, 1-alpha/2))).
The interval is inverted around the UNCORRECTED jump or ratio; the REPORTABLE
bootstrap estimate equals the analytic bias-corrected estimate. The interval
therefore need not be centered on the displayed point estimate.

Matching seed numbers in R and Stata do NOT imply matching bootstrap draws.
Use common supplied multiplier matrices to test numerical agreement.

Regularization counts observations, not distinct x values. A scalar h uses
pooled cutoff regularization; a two-element h or a side override uses separate
cutoff regularization. Auxiliary fits may enlarge further. All these details,
including the sharp/joint variance and bootstrap-support conventions after a
side recheck, are preserved from R 2.0. e(h) is the reported bandwidth pair;
e(h_fit) is the actual pair in the cutoff local fits, after a side recheck.
They can differ under regularization. This port does not silently change the
R implementation. Do not trim the input sample to the selected h windows.

KNOWN RESTRICTION
-----------------
R 2.0's g.loo path removes zero x and zero outcome/component entries separately.
That can misalign rows, especially for a binary treatment. This development
port rejects gloo when any selected x or component equals zero, instead of
silently returning different results. Without zeros it uses index-based
leave-one-out fits matching R. The default gloo=off works with binary treatment
and observations exactly at the cutoff. No change was made to the R code.

STORED RESULTS
--------------
e(b), e(V)                 Analytic bias-corrected coefficient and variance.
e(sample)                  Complete selected sample, NOT just observations in h.
e(estimate_irbc), e(se_irbc), e(ci_irbc)
e(estimate_irbc_boot), e(ci_irbc_boot)  Present with bootstrap.
e(tau_uncorrected), e(bias_irbc), e(first_stage)
e(h_requested), e(h), e(h_fit), e(Q), e(N_h)  Side pairs, left then right.
e(components), e(component_bias), e(component_estimates), e(linear_weights)
e(cov_irbc), e(cov_us)       Joint component covariance matrices.
e(gamma)                   Common slopes, rows=covariates; columns=outcome,
                           then treatment if fuzzy. Absent without covariates.
e(se_us), e(boot_probs), e(boot_mhat), e(boot_B), e(boot_success)
e(table_irbc)              Compact result matrix. Bootstrap SE is blank because
                           inference uses prepivoted quantiles, not a second SE.
e(N_selected), e(N_omitted), e(N)  Missing-data accounting.
e(design), e(residuals), e(qmethod), e(cmdline), e(validation)

Components are ordered outcome, treatment (if present), covariates in covs()
order. Bootstrap draws and component residuals are retained in the Mata
structure irbc_last, avoiding Stata's smaller matrix-size limit:

  mata: irbc_last.boot_t
  mata: irbc_last.boot_components
  mata: irbc_last.boot_first_stage
  mata: irbc_last.residuals

These diagnostics describe the last successful call, are overwritten by the
next successful call, and are not saved by estimates store. nokeep removes
bootstrap draws, not the reported bootstrap interval. Residual arrays are
local computational inputs, not full-sample fitted-value diagnostics.

VALIDATION
----------
tests/fixtures contains data, supplied wild weights, and expected values
generated by an unchanged snapshot of R 2.0. In a fresh Stata session, cd to
this folder and run:

  do tests/run_tests.do

The script replaces the current data. It compares estimates, residuals, Qs,
slopes, covariances, every bootstrap statistic, and interval endpoints, then
checks native RNG reproducibility, missing-data handling, option validation,
and interface behavior. A passing log is required before release.

Developer-only regeneration (base R, no extra R packages):
  Rscript tests/generate_fixtures.R
  Rscript tests/check_bundle.R

The generator uses the local R 2.0 source when available and otherwise uses
the included reference snapshot. The expected source hashes are checked.

This also compares an independent R evaluation of the port's matrix formulas
with the R reference. A PASS here is an algebra check, NOT a Stata validation.
See tests/FORMULA_AUDIT.txt and VALIDATION.txt for the actual status.

OUT OF SCOPE
------------
No direct/debiased bootstrap, clustered SEs, kink designs, internal optimal
bandwidth selection, weak-first-stage-robust confidence sets, or factor-variable
expansion. Missing analysis variables are excluded in Stata's usual manner;
the R function instead expects callers to prepare complete cases. Singular
fits fail explicitly; the command does not silently drop collinear covariates.

REFERENCES
----------
Cavaliere, G., Goncalves, S., Nielsen, M. O., and Zanelli, E.:
  Improved Inference for Regression Discontinuity Designs (current manuscript).
  Earlier public version: https://arxiv.org/abs/2512.00566.
Calonico, S., Cattaneo, M. D., Farrell, M. H., and Titiunik, R. (2019):
  Regression Discontinuity Designs Using Covariates, Review of Economics
  and Statistics 101(3), 442-451.

Implementation references (official Stata documentation):
  https://www.stata.com/manuals/m-1ado.pdf
  https://www.stata.com/manuals/m-2struct.pdf
  https://www.stata.com/manuals/m-5cholinv.pdf
  https://www.stata.com/manuals/preturn.pdf
