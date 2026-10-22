* Development examples. Stata execution has not yet been validated.
* Run in a new session from this folder. This replaces the data in memory.
version 16.0
clear all
set more off
adopath ++ "."
set seed 17041998
set obs 400
generate double x = 2*runiform()-1
generate double z1 = .3*x+rnormal()
generate double z2 = -.2*x+rnormal()
generate double d = runiform()<(.2+.6*(x>=0))
generate double y = 1+x+x^2+.8*d+.5*z1-.2*z2+rnormal()

* h is illustrative and supplied explicitly, not selected by these examples.
irbcrd y x, h(.6)
irbcrd y x, h(.6) fuzzy(d)
irbcrd y x, h(.6) covs(z1 z2)
irbcrd y x, h(.6) fuzzy(d) covs(z1 z2)

irbcrd y x, h(.6) fuzzy(d) covs(z1 z2) residuals(cct-hc3) ///
    qmethod(sample) bootstrap reps(399) wild(rademacher) seed(17041998)
matrix list e(ci_irbc)
matrix list e(ci_irbc_boot)

irbcrd y x, h(.5 .7) kernel(triangular epanechnikov) ///
    qmethod(asymptotic) bootstrap reps(399) wild(mammen) seed(17041998)

* The bootstrap first-stage denominator option is only for plain fuzzy RD.
irbcrd y x, h(.6) fuzzy(d) bootstrap fuzzyweights(bootstrap) ///
    reps(399) residuals(cct-hc3) seed(17041998)
