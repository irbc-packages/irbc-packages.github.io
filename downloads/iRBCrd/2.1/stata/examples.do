* With the extracted package folder as the working directory.
* These examples replace the dataset in memory; h values are illustrative.
version 16.0
clear all
set more off
set seed 17041998
set obs 500
adopath ++ "."
generate double x = 2*runiform()-1
generate double z1 = .3*x+rnormal()
generate double z2 = -.2*x+rnormal()
generate double d = .3+.6*(x>=0)+.8*x*(x>=0)+.15*z1+.15*rnormal()
generate double y = 1+x+x^2+.8*d+.4*z1-.2*z2+.3*rnormal()

irbcrd y x, h(.6)
irbcrd y x, h(.6) fuzzy(d)
irbcrd y x, h(.6) covs(z1 z2)
irbcrd y x, h(.6) fuzzy(d) covs(z1 z2)

irbcrd y x, h(.5 .7) kernel(tri epa) res(cct-hc3) q_method(asymptotic) bwcheck(0)
irbcrd y x, h(.6) alpha(.01) boot_ci boot_b(999) boot_seed(17041998)
irbcrd y x, h(.6) fuzzy(d) covs(z1 z2) boot_ci boot_b(999) ///
    boot_seed(17041998) boot_wild(mammen) boot_keep ///
    fuzzy_boot_weights(bootstrap) gamma_boot_weights(bootstrap)
matrix list e(ci_irbc_boot)
mata: st_matrix("e(boot_t)")[1..5,.]

* Every design also supports kink, with p(2) by default.
irbcrd y x, h(.7) kink
irbcrd y x, h(.7) kink fuzzy(d)
irbcrd y x, h(.7) kink covs(z1 z2)
irbcrd y x, h(.7) kink fuzzy(d) covs(z1 z2) ///
    res(cct-hc3) q_method(asymptotic) boot_ci boot_b(999) ///
    boot_seed(17041998) fuzzy_boot_weights(bootstrap) gamma_boot_weights(bootstrap)
