######################################################################################
## iRBCrd: improved robust bias correction for RD designs
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 2.0 (22.10.2026)
## Contacts: giuseppe.cavaliere@unibo.it; silvia.goncalves@mcgill.ca;
##           mon@econ.au.dk; ez@econ.au.dk
######################################################################################


iRBCrd <- function(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
                   alpha = NULL, fast = NULL, g.loo = NULL,
                   boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                   boot.wild = NULL, boot.keep = NULL,
                   bwcheck = .bwcheck_default,
                   h.plus = NULL, h.minus = NULL,
                   kersel.plus = NULL, kersel.minus = NULL,
                   q.method = NULL,
                   fuzzy = NULL, covs = NULL,
                   fuzzy.boot.weights = NULL, gamma.boot.weights = NULL,
                   d = NULL) {

  if (!is.null(fuzzy) && !is.null(d)) {
    stop("Supply the treatment through either 'fuzzy' or its alias 'd', not both.")
  }
  if (is.null(fuzzy)) fuzzy <- d
  has.fuzzy <- !is.null(fuzzy)
  has.covs <- !is.null(covs)

  if (!has.fuzzy && !is.null(fuzzy.boot.weights)) {
    stop("'fuzzy.boot.weights' requires a fuzzy RD design: supply 'fuzzy' or 'd'.")
  }
  if (!(has.fuzzy && has.covs) && !is.null(gamma.boot.weights)) {
    stop("'gamma.boot.weights' is supported only for covariate-adjusted fuzzy RD: supply treatment and 'covs'.")
  }

  args <- list(
    y = y, x = x, h = h, c = c, p = p, kersel = kersel, res = res,
    alpha = alpha, fast = fast, g.loo = g.loo,
    boot.ci = boot.ci, boot.B = boot.B, boot.seed = boot.seed,
    boot.wild = boot.wild, boot.keep = boot.keep, bwcheck = bwcheck,
    h.plus = h.plus, h.minus = h.minus,
    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
    q.method = q.method
  )

  if (has.fuzzy) {
    args$d <- fuzzy
    args$fuzzy.boot.weights <- fuzzy.boot.weights
    if (has.covs) {
      args$covs <- covs
      args$gamma.boot.weights <- gamma.boot.weights
      return(do.call(iRBCcfrd, args))
    }
    return(do.call(iRBCfrd, args))
  }
  if (has.covs) {
    args$covs <- covs
    return(do.call(iRBCcrd, args))
  }
  do.call(.irbcrd_sharp, args)
}


.irbcrd_sharp <- function(y, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
                   alpha = NULL, fast = NULL, g.loo = NULL,
                   boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                   boot.wild = NULL, boot.keep = NULL,
                   bwcheck = .bwcheck_default,
                   h.plus = NULL, h.minus = NULL,
                   kersel.plus = NULL, kersel.minus = NULL,
                   q.method = NULL) {

  # y       <- dependent variable
  # x       <- regressor
  # h       <- bandwidth; scalar or c(left, right), optionally named
  # c       <- cutoff
  # p       <- polynomial order                 (default = 1)
  # kersel  <- kernel function used; scalar or c(left, right) (default = "tri")
  # res     <- residuals class used             (default = "irbc-hc3")
  # alpha   <- significance level               (default = 0.05)
  # fast    <- fast std.error computation       (default = FALSE)
  # g.loo   <- leave-one-out Bias estimation    (default = FALSE)
  # boot.ci <- fixed-regressor wild bootstrap CI (default = FALSE)
  # boot.B  <- bootstrap repetitions             (default = 999 if boot.ci = TRUE)
  # q.method <- Q construction: "sample" or "asymptotic"
  #             (default: "sample" if fast = FALSE, "asymptotic" if fast = TRUE)
  # bwcheck <- minimum observations used for bandwidth regularization;
  #            set to NULL to disable regularization (default = 21)


  ############################# Error Checks #############################

  if (!is.numeric(h) || !(length(h) %in% c(1, 2)) || any(!is.finite(h)) || any(h <= 0)) {
    stop("'h' must be a positive number or a length-2 vector of positive numbers")
  }


  if (!is.numeric(c) & !is.null(c)) {
    stop("'c' must be numeric")
  } else if (is.null(c)) {
    c <- 0
  }


  if (length(p) == 0) {
    p  <- 1
  } else if (length(p) > 1) {
    stop("Polynomial order p incorrectly specified.\n")
  }


  if (!is.null(kersel) && (!is.character(kersel) || !(length(kersel) %in% c(1, 2)) ||
                           any(!(kersel %in% c("uni", "uniform", "tri", "triangular", "epa", "epanechnikov"))))) {
    stop("Kernel function incorrectly specified.\n")
  } else if (is.null(kersel)) {
    kersel <- "tri"
  }

  res <- .irbcrd_normalize_res(res)


  if (is.null(alpha))  {
    alpha  <- 0.05
  } else if (alpha <= 0 | alpha >= 1) {
    stop("significance level incorrectly specified.\n")
  }


  if (is.null(fast)) {
    fast <- FALSE
  } else if (!is.logical(fast) || length(fast) > 1) {
    stop("'fast' must be a single TRUE/FALSE")
  }


  if (is.null(g.loo)) {
    g.loo <- FALSE
  } else if (!is.logical(g.loo) || length(g.loo) > 1) {
    stop("'g.loo' must be a single TRUE/FALSE")
  }


  if (is.null(boot.ci)) {
    boot.ci <- FALSE
  } else if (!is.logical(boot.ci) || length(boot.ci) > 1) {
    stop("'boot.ci' must be a single TRUE/FALSE")
  }


  if (is.null(boot.B)) {
    boot.B <- if (boot.ci) 999 else 0
  } else if (!is.numeric(boot.B) || length(boot.B) > 1 || boot.B < 1 || boot.B != floor(boot.B)) {
    stop("'boot.B' must be a positive integer")
  }
  boot.B <- as.integer(boot.B)


  if (!is.null(boot.seed) && (!is.numeric(boot.seed) || length(boot.seed) > 1)) {
    stop("'boot.seed' must be a single numeric value")
  }


  boot.wild <- .irbcrd_normalize_boot_wild(boot.wild)


  if (is.null(boot.keep)) {
    boot.keep <- TRUE
  } else if (!is.logical(boot.keep) || length(boot.keep) > 1) {
    stop("'boot.keep' must be a single TRUE/FALSE")
  }


  bwcheck <- .validate_bwcheck(bwcheck)
  q.method <- .irbcrd_normalize_q_method(q.method, fast = fast)


  parse_two_sided <- function(value, value.name) {
    if (length(value) == 1) {
      return(list(minus = value[[1]], plus = value[[1]]))
    }

    nms <- names(value)
    if (!is.null(nms) && any(nzchar(nms))) {
      nms <- tolower(nms)
      plus.index <- which(nms %in% c("plus", "right", "p", "r"))
      minus.index <- which(nms %in% c("minus", "left", "m", "l"))
      if (length(plus.index) != 1 || length(minus.index) != 1) {
        stop("Named length-2 '", value.name, "' must identify one left/minus and one right/plus value.")
      }
      return(list(minus = value[[minus.index]], plus = value[[plus.index]]))
    }

    list(minus = value[[1]], plus = value[[2]])
  }

  h.scalar.default <- length(h) == 1 && is.null(h.minus) && is.null(h.plus)
  h.side <- parse_two_sided(h, "h")
  if (!is.null(h.minus)) {
    if (!is.numeric(h.minus) || length(h.minus) != 1 || !is.finite(h.minus) || h.minus <= 0) {
      stop("'h.minus' must be a positive number")
    }
    h.side$minus <- h.minus
  }
  if (!is.null(h.plus)) {
    if (!is.numeric(h.plus) || length(h.plus) != 1 || !is.finite(h.plus) || h.plus <= 0) {
      stop("'h.plus' must be a positive number")
    }
    h.side$plus <- h.plus
  }

  k.side <- parse_two_sided(kersel, "kersel")
  if (!is.null(kersel.minus)) k.side$minus <- kersel.minus
  if (!is.null(kersel.plus)) k.side$plus <- kersel.plus
  k.side$minus <- .canonical_kernel(k.side$minus)
  k.side$plus <- .canonical_kernel(k.side$plus)


  ##########################  Initialization  ##########################

  # Sample size
  n <- length(x)

  # Specify variables at each side of the cutoff
  indp <- x >= c
  indm <- x < c
  xp <- x[indp]
  yp <- y[indp]
  xm <- x[indm]
  ym <- y[indm]


  # Side-specific bandwidths and kernels
  if (h.scalar.default) {
    hp <- hm <- .regularize_h(x = x, x0 = c, h = h.side$plus, bwcheck = bwcheck)
  } else {
    hp <- .regularize_h(x = xp, x0 = c, h = h.side$plus, bwcheck = bwcheck)
    hm <- .regularize_h(x = xm, x0 = c, h = h.side$minus, bwcheck = bwcheck)
  }
  kerselp <- k.side$plus
  kerselm <- k.side$minus
  h.out <- if (isTRUE(all.equal(hp, hm))) hp else c(minus = hm, plus = hp)
  kersel.out <- if (identical(kerselp, kerselm)) kerselp else c(minus = kerselm, plus = kerselp)


  ##########################  ATE and Bias estimation    ##########################

  # ATE estimation
  ghatp   <- LP.est(y = yp, x = xp, x0 = c, h = hp, p = p, v = 0,
                    kersel = kerselp, bwcheck = bwcheck)
  ghatm   <- LP.est(y = ym, x = xm, x0 = c, h = hm, p = p, v = 0,
                    kersel = kerselm, bwcheck = bwcheck)
  tau.hat <- ghatp - ghatm

  # Initialization of bias estimation
  # The iRBC bias correction uses auxiliary local fits over the 2h window.
  U.Z <- (indp & abs((x - c) / (2 * hp)) < 1) |
    (indm & abs((x - c) / (2 * hm)) < 1)
  nonzero_indices <- which(U.Z != 0)
  ghat_vec <- numeric(n)
  if (g.loo == TRUE) ghat_vec.loo <- numeric(n)
  if (q.method == "sample") C.xi <- numeric(n)

  # Loop for bias estimation
  for (ix in seq_along(nonzero_indices)) {
    i <- nonzero_indices[ix]
    x0i <- x[i]
    if (x0i >= c) {

      ghat_vec[i] <- LP.est(y = yp, x = xp, x0 = x0i, h = hp, p = p, v = 0,
                            kersel = kerselp, bwcheck = bwcheck)

      if (q.method == "sample") C.xi[i] <- LP.est(y = ((xp - x0i) / hp)^(p + 1), x = xp,
                                                  x0 = x0i, h = hp, p = p, v = 0,
                                                  kersel = kerselp, bwcheck = bwcheck)

      if (g.loo == TRUE) {

        xp.loo <- x * (x >= c)
        yp.loo <- y * (x >= c)
        xp.loo <- xp.loo[-i]
        yp.loo <- yp.loo[-i]
        xp.loo <- xp.loo[xp.loo != 0]
        yp.loo <- yp.loo[yp.loo != 0]

        ghat_vec.loo[i] <- LP.est(y = yp.loo, x = xp.loo, x0 = x0i, h = hp,
                                  p = p, v = 0, kersel = kerselp,
                                  bwcheck = bwcheck)

      }

    } else if (x0i < c) {

      ghat_vec[i] <- LP.est(y = ym, x = xm, x0 = x0i, h = hm, p = p, v = 0,
                            kersel = kerselm, bwcheck = bwcheck)

      if (q.method == "sample") C.xi[i] <- LP.est(y = ((xm - x0i) / hm)^(p + 1), x = xm,
                                                  x0 = x0i, h = hm, p = p, v = 0,
                                                  kersel = kerselm, bwcheck = bwcheck)

      if (g.loo == TRUE) {

        xm.loo <- x * (x < c)
        ym.loo <- y * (x < c)
        xm.loo <- xm.loo[-i]
        ym.loo <- ym.loo[-i]
        xm.loo <- xm.loo[xm.loo != 0]
        ym.loo <- ym.loo[ym.loo != 0]

        ghat_vec.loo[i] <- LP.est(y = ym.loo, x = xm.loo, x0 = x0i, h = hm,
                                  p = p, v = 0, kersel = kerselm,
                                  bwcheck = bwcheck)

      }
    }
  }

  # Estimation of Q at each side of the cutoff
  if (q.method == "sample") {
    Cp.LP <- LP.est(y = C.xi[indp], x = xp, x0 = c, h = hp, p = p, v = 0,
                    kersel = kerselp, bwcheck = bwcheck)
    Cm.LP <- LP.est(y = C.xi[indm], x = xm, x0 = c, h = hm, p = p, v = 0,
                    kersel = kerselm, bwcheck = bwcheck)
    Cp    <- LP.est(y = ((xp - c) / hp)^(p + 1), x = xp, x0 = c, h = hp,
                    p = p, v = 0, kersel = kerselp, bwcheck = bwcheck)
    Cm    <- LP.est(y = ((xm - c) / hm)^(p + 1), x = xm, x0 = c, h = hm,
                    p = p, v = 0, kersel = kerselm, bwcheck = bwcheck)
    Qp    <- Cp / Cp.LP
    Qm    <- Cm / Cm.LP
  } else if (q.method == "asymptotic") {
    Qp <- .irbcrd_Q_fast(kerselp, p = p)
    Qm <- .irbcrd_Q_fast(kerselm, p = p)
  }

  # Bias Estimation
  if (g.loo == TRUE) {
    residualsp_gs <- ghat_vec.loo[indp] - c(ghatp)
    residualsm_gs <- ghat_vec.loo[indm] - c(ghatm)
  } else if (g.loo == FALSE) {
    residualsp_gs <- ghat_vec[indp] - c(ghatp)
    residualsm_gs <- ghat_vec[indm] - c(ghatm)
  }

  Bhatp_lp      <- LP.est(y = residualsp_gs, x = xp, x0 = c, h = hp, p = p, v = 0,
                          kersel = kerselp, bwcheck = bwcheck)
  Bhatm_lp      <- LP.est(y = residualsm_gs, x = xm, x0 = c, h = hm, p = p, v = 0,
                          kersel = kerselm, bwcheck = bwcheck)
  Bhat_irbc      <- Bhatp_lp * Qp - Bhatm_lp * Qm

  # De-biased ATE
  tau.hat_irbc <- tau.hat - Bhat_irbc


  #####################  Standard errors and Confidence Intervals  #####################

  # Residuals
  if (res == "loo") {
    if (g.loo == TRUE) {
      epshat <- y - ghat_vec.loo
    } else {
      epshat <- y - ghat_vec
    }

    epshatp <- epshat[indp]
    epshatm <- epshat[indm]
  } else if (res %in% c("cct-hc0", "cct-hc1", "cct-hc2", "cct-hc3")) {
    rp.pp1 <- .poly_basis(xp - c, p + 1)
    rm.pp1 <- .poly_basis(xm - c, p + 1)
    Kp.X <- K(u = (xp - c) / hp, kersel = kerselp) / hp
    Km.X <- K(u = (xm - c) / hm, kersel = kerselm) / hm
    invGp.q <- qrXXinv(sqrt(Kp.X) * rp.pp1)
    invGm.q <- qrXXinv(sqrt(Km.X) * rm.pp1)
    betap.q <- invGp.q %*% crossprod(rp.pp1 * Kp.X, yp)
    betam.q <- invGm.q %*% crossprod(rm.pp1 * Km.X, ym)

    if (res == "cct-hc0") {
      epshatp <- yp - rp.pp1 %*% betap.q
      epshatm <- ym - rm.pp1 %*% betam.q
    } else {
      Qp.vec <- rowSums((rp.pp1 %*% invGp.q) * (rp.pp1 * Kp.X))
      Qm.vec <- rowSums((rm.pp1 %*% invGm.q) * (rm.pp1 * Km.X))
      if (res == "cct-hc1") {
        nhp <- sum(Kp.X != 0)
        nhm <- sum(Km.X != 0)
        if (nhp <= p + 2 || nhm <= p + 2) {
          stop("HC1 residual correction requires more nonzero-kernel observations than q + 1 on each side.")
        }
        epshatp <- (yp - rp.pp1 %*% betap.q) / (((nhp - (p + 2)) / nhp)^(0.5))
        epshatm <- (ym - rm.pp1 %*% betam.q) / (((nhm - (p + 2)) / nhm)^(0.5))
      } else if (res == "cct-hc2") {
        epshatp <- (yp - rp.pp1 %*% betap.q) / ((1 - Qp.vec)^(0.5))
        epshatm <- (ym - rm.pp1 %*% betam.q) / ((1 - Qm.vec)^(0.5))
      } else if (res == "cct-hc3") {
        epshatp <- (yp - rp.pp1 %*% betap.q) / (1 - Qp.vec)
        epshatm <- (ym - rm.pp1 %*% betam.q) / (1 - Qm.vec)
      }
    }
  } else if (res %in% c("irbc-hc0", "irbc-hc1", "irbc-hc2", "irbc-hc3")) {
    irbc.res.p <- .irbcrd_residuals_side(y = yp, x = xp, x.cut = c,
                                              h = hp, p = p,
                                              kersel = kerselp, res = res,
                                              q.method = q.method,
                                              bwcheck = bwcheck)
    irbc.res.m <- .irbcrd_residuals_side(y = ym, x = xm, x.cut = c,
                                              h = hm, p = p,
                                              kersel = kerselm, res = res,
                                              q.method = q.method,
                                              bwcheck = bwcheck)
    epshatp <- irbc.res.p$epshat
    epshatm <- irbc.res.m$epshat
  }

  # iRBC Standard errors
  if (fast == 0) {
    sep_irbc <- viRBC(x = xp, x0 = c, p = p, h = hp, kersel = kerselp,
                      epshat = epshatp, Q = Qp, bwcheck = bwcheck) /
      sqrt(length(xp) * hp)
    sem_irbc <- viRBC(x = xm, x0 = c, p = p, h = hm, kersel = kerselm,
                      epshat = epshatm, Q = Qm, bwcheck = bwcheck) /
      sqrt(length(xm) * hm)
  } else if (fast == 1) {
    Kratio.fast <- function(kersel) {
      if (kersel == "tri") {
        0.84
      } else if (kersel == "epa") {
        0.83
      } else if (kersel == "uni") {
        0.86
      }
    }
    sep_irbc <- Kratio.fast(kerselp) *
      (CCT.se(x = xp, x0 = c, h = hp, p = p, kersel = kerselp,
              epshat = epshatp, bwcheck = bwcheck)$se.rb)
    sem_irbc <- Kratio.fast(kerselm) *
      (CCT.se(x = xm, x0 = c, h = hm, p = p, kersel = kerselm,
              epshat = epshatm, bwcheck = bwcheck)$se.rb)
  }

  se_irbc <- sqrt(sep_irbc^2 + sem_irbc^2)


  # Confidence Intervals
  ci_irbc <- c(tau.hat_irbc - qnorm(1 - alpha / 2) * se_irbc, tau.hat_irbc - qnorm(alpha / 2) * se_irbc)


  #####################  Fixed-regressor wild bootstrap CI  #####################

  if (boot.ci == TRUE) {
    if (!is.null(boot.seed)) set.seed(as.integer(boot.seed))

    epshat <- numeric(n)
    epshat[indp] <- as.numeric(epshatp)
    epshat[indm] <- as.numeric(epshatm)

    sep_us <- Qp * LP.se(x = xp, x0 = c, h = hp, p = p, kersel = kerselp,
                         epshat = epshatp, bwcheck = bwcheck)
    sem_us <- Qm * LP.se(x = xm, x0 = c, h = hm, p = p, kersel = kerselm,
                         epshat = epshatm, bwcheck = bwcheck)
    se_us <- sqrt(sep_us^2 + sem_us^2)

    mhat <- as.numeric(se_irbc / se_us)
    if (!is.finite(mhat) || mhat <= 0) {
      stop("Bootstrap CI cannot be computed because the non-bias-corrected standard error is not positive.")
    }

    H.inv <- function(u) pnorm(mhat * qnorm(u))
    boot.probs <- H.inv(c(alpha / 2, 1 - alpha / 2))

    h.boot.p <- .regularize_h(x = xp, x0 = c, h = hp, bwcheck = bwcheck)
    h.boot.m <- .regularize_h(x = xm, x0 = c, h = hm, bwcheck = bwcheck)
    boot.local <- (indp & abs(x - c) <= h.boot.p) |
      (indm & abs(x - c) <= h.boot.m)
    if (!any(boot.local & indp) || !any(boot.local & indm)) {
      stop("Bootstrap CI requires observations on both sides of the cutoff inside the bandwidth.")
    }

    ghat.boot.mean <- numeric(n)
    boot.local.indices <- which(boot.local)
    for (ix in seq_along(boot.local.indices)) {
      i <- boot.local.indices[ix]
      if (x[i] >= c) {
        ghat.boot.mean[i] <- LP.est(y = yp, x = xp, x0 = x[i], h = hp,
                                    p = p, v = 0, kersel = kerselp,
                                    bwcheck = bwcheck)
      } else {
        ghat.boot.mean[i] <- LP.est(y = ym, x = xm, x0 = x[i], h = hm,
                                    p = p, v = 0, kersel = kerselm,
                                    bwcheck = bwcheck)
      }
    }

    boot.t <- rep(NA_real_, boot.B)

    for (bb in seq_len(boot.B)) {
      e.boot <- .irbcrd_wild_draws(length(boot.local.indices), boot.wild)

      y.boot <- y
      y.boot[boot.local.indices] <- ghat.boot.mean[boot.local.indices] +
        epshat[boot.local.indices] * e.boot

      ghatp.boot <- LP.est(y = y.boot[indp], x = xp, x0 = c, h = hp,
                           p = p, v = 0, kersel = kerselp,
                           bwcheck = bwcheck)
      ghatm.boot <- LP.est(y = y.boot[indm], x = xm, x0 = c, h = hm,
                           p = p, v = 0, kersel = kerselm,
                           bwcheck = bwcheck)

      boot.t[bb] <- as.numeric(Qp * (ghatp.boot - ghatp) - Qm * (ghatm.boot - ghatm))
    }

    boot.t.ok <- boot.t[is.finite(boot.t)]
    boot.success <- length(boot.t.ok)
    if (boot.success < 2) {
      warning("Fewer than two successful bootstrap repetitions; bootstrap CI is NA.")
      ci_irbc_boot <- c(NA_real_, NA_real_)
    } else {
      q.boot <- as.numeric(quantile(boot.t.ok, probs = boot.probs, names = FALSE, type = 7))
      ci_irbc_boot <- c(as.numeric(tau.hat) - q.boot[2], as.numeric(tau.hat) - q.boot[1])
    }
  }

  
  # Output
  out <- list(n = n, h = h.out, h.requested = h.side, bwcheck = bwcheck,
              h.plus = hp, h.minus = hm,
              kersel = kersel.out, kersel.plus = kerselp, kersel.minus = kerselm,
              res = res, q.method = q.method,
              Q.plus = Qp, Q.minus = Qm, tau.hat = tau.hat, tau.hat_irbc = tau.hat_irbc,
              estimate_irbc = as.numeric(tau.hat_irbc),
              Bhat.plus.lp = Bhatp_lp, Bhat.minus.lp = Bhatm_lp,
              Bhat.irbc = Bhat_irbc, se_irbc = se_irbc,
              ci_irbc = ci_irbc, epshat.plus = epshatp, epshat.minus = epshatm)

  if (boot.ci == TRUE) {
    out <- c(out, list(estimate_irbc_boot = as.numeric(tau.hat_irbc),
                       ci_irbc_boot = ci_irbc_boot,
                       boot.B = boot.B,
                       boot.success = boot.success,
                       boot.impl = "plug-in",
                       boot.probs = boot.probs,
                       boot.mhat = mhat,
                       boot.wild = boot.wild,
                       se_us = se_us))
    if (boot.keep == TRUE) out$boot.t <- boot.t
  }
  
  return(out)

}


.irbcrd_normalize_res <- function(res) {
  if (is.null(res)) {
    return("irbc-hc3")
  }
  if (!is.character(res) || length(res) != 1) {
    stop("Residuals incorrectly specified.\n")
  }

  key <- tolower(trimws(res))
  key <- gsub("[[:space:]_]", "-", key)

  if (key == "loo") {
    return("loo")
  }
  if (key %in% c("hc0", "hc1", "hc2", "hc3")) {
    return(paste0("irbc-", key))
  }
  if (key %in% c("irbc-hc0", "irbc-hc1", "irbc-hc2", "irbc-hc3")) {
    return(key)
  }
  if (key %in% c("cct-hc0", "cct-hc1", "cct-hc2", "cct-hc3",
                 "rbc-hc0", "rbc-hc1", "rbc-hc2", "rbc-hc3")) {
    return(sub("^rbc-", "cct-", key))
  }

  stop("Residuals incorrectly specified.\n")
}


.irbcrd_normalize_q_method <- function(q.method, fast) {
  if (is.null(q.method)) {
    return(if (isTRUE(fast)) "asymptotic" else "sample")
  }
  if (!is.character(q.method) || length(q.method) != 1) {
    stop("'q.method' must be NULL, 'sample', or 'asymptotic'")
  }

  key <- tolower(trimws(q.method))
  key <- gsub("[[:space:]_\\.-]", "", key)

  if (key %in% c("sample", "samples", "finite", "finitesample", "qn")) {
    if (isTRUE(fast)) {
      stop("'q.method = \"sample\"' cannot be combined with fast = TRUE.")
    }
    return("sample")
  }
  if (key %in% c("asymptotic", "asym", "population", "limit", "fast")) {
    return("asymptotic")
  }

  stop("'q.method' must be NULL, 'sample', or 'asymptotic'")
}


.irbcrd_Q_fast <- function(kersel, p = 1) {
  kersel <- .canonical_kernel(kersel)
  if (!isTRUE(all.equal(p, 1))) {
    stop("Asymptotic Q constants are currently implemented for p = 1 only.")
  }
  if (kersel == "tri") {
    return(1.4082)
  }
  if (kersel == "epa") {
    return(1.3571)
  }
  if (kersel == "uni") {
    return(1.2)
  }
  stop("Invalid kernel specification.")
}


.irbcrd_residuals_side <- function(y, x, x.cut, h, p, kersel, res,
                                        q.method,
                                        bwcheck = .bwcheck_default) {
  hc <- sub("^irbc-", "", res)
  n <- length(x)
  epshat <- numeric(n)
  fitted <- rep(NA_real_, n)
  leverage <- rep(NA_real_, n)
  local.n <- rep(NA_integer_, n)
  Q.point <- rep(NA_real_, n)

  h.cut <- .regularize_h(x = x, x0 = x.cut, h = h, bwcheck = bwcheck)
  active <- which(abs((x - x.cut) / (2 * h.cut)) <= 1)
  if (length(active) == 0L) {
    return(list(epshat = epshat, fitted = fitted, leverage = leverage,
                Q.point = Q.point, active = active))
  }

  for (jj in active) {
    h.j <- .regularize_h(x = x, x0 = x[jj], h = h, bwcheck = bwcheck)
    fitted[jj] <- LP.est(y = y, x = x, x0 = x[jj], h = h,
                         p = p, v = 0, kersel = kersel,
                         bwcheck = bwcheck)
    local.n[jj] <- sum(K(u = (x - x[jj]) / h.j, kersel = kersel) != 0)
    leverage[jj] <- as.numeric(wix0(x = x, xi = x[jj], x0 = x[jj],
                                    h = h, p = p, v = 0,
                                    kersel = kersel,
                                    bwcheck = bwcheck)) / (n * h.j)
  }

  raw <- as.numeric(y[active] - fitted[active])
  if (hc == "hc0") {
    epshat[active] <- raw
  } else if (hc == "hc1") {
    df.factor <- (local.n[active] - (p + 1)) / local.n[active]
    if (any(!is.finite(df.factor)) || any(df.factor <= 0)) {
      stop("iRBC HC1 residual correction requires more local observations than p + 1.")
    }
    epshat[active] <- raw / sqrt(df.factor)
  } else if (hc %in% c("hc2", "hc3")) {
    adj <- 1 - leverage[active]
    if (any(!is.finite(adj)) || any(adj <= 0)) {
      stop("iRBC HC2/HC3 residual correction encountered non-positive leverage adjustment.")
    }
    if (hc == "hc2") {
      epshat[active] <- raw / sqrt(adj)
    } else {
      epshat[active] <- raw / adj
    }
  } else {
    stop("Unknown iRBC residual correction.")
  }

  list(epshat = epshat, fitted = fitted, leverage = leverage,
       Q.point = Q.point, active = active)
}


.irbcrd_normalize_boot_wild <- function(boot.wild) {
  if (is.null(boot.wild)) {
    return("rademacher")
  }
  if (!is.character(boot.wild) || length(boot.wild) > 1) {
    stop("'boot.wild' must be 'normal', 'rademacher', or 'mammen'")
  }

  key <- tolower(trimws(boot.wild))
  key <- gsub("[[:space:]_\\.-]", "", key)
  key <- gsub("'", "", key, fixed = TRUE)

  if (key %in% c("normal", "gaussian", "standardnormal")) {
    return("normal")
  }
  if (key %in% c("rademacher", "rad")) {
    return("rademacher")
  }
  if (key %in% c("mammen", "mammens", "mammen1993", "mammentwopoint",
                 "mammenwild")) {
    return("mammen")
  }

  stop("'boot.wild' must be 'normal', 'rademacher', or 'mammen'")
}


.irbcrd_normalize_fuzzy_boot_weights <- function(fuzzy.boot.weights) {
  if (is.null(fuzzy.boot.weights)) {
    return("fixed")
  }
  if (!is.character(fuzzy.boot.weights) || length(fuzzy.boot.weights) > 1) {
    stop("'fuzzy.boot.weights' must be 'fixed' or 'bootstrap'")
  }

  key <- tolower(trimws(fuzzy.boot.weights))
  key <- gsub("[[:space:]_\\.-]", "", key)

  if (key %in% c("fixed", "eq615", "equation615", "linear", "paper")) {
    return("fixed")
  }
  if (key %in% c("bootstrap", "boot", "star", "starred",
                 "bootstrapweights", "bootweights")) {
    return("bootstrap")
  }

  stop("'fuzzy.boot.weights' must be 'fixed' or 'bootstrap'")
}


.irbcrd_normalize_gamma_boot_weights <- function(gamma.boot.weights) {
  if (is.null(gamma.boot.weights)) {
    return("fixed")
  }
  if (!is.character(gamma.boot.weights) || length(gamma.boot.weights) != 1L ||
      is.na(gamma.boot.weights)) {
    stop("'gamma.boot.weights' must be NULL or 'fixed'; covariate-adjusted fuzzy RD keeps the original-sample covariate coefficients fixed.")
  }

  key <- tolower(trimws(gamma.boot.weights))
  key <- gsub("[[:space:]_\\.-]", "", key)

  if (key %in% c("fixed", "original", "sample", "hat")) {
    return("fixed")
  }
  stop("'gamma.boot.weights' must be NULL or 'fixed'; covariate-adjusted fuzzy RD keeps the original-sample covariate coefficients fixed.")
}


.irbcrd_wild_draws <- function(n, wild) {
  n <- as.integer(n)
  if (n < 0) stop("'n' must be non-negative")

  wild <- .irbcrd_normalize_boot_wild(wild)
  if (wild == "normal") {
    return(rnorm(n))
  }
  if (wild == "rademacher") {
    return(sample(c(-1, 1), n, replace = TRUE))
  }

  sqrt5 <- sqrt(5)
  values <- c(-(sqrt5 - 1) / 2, (sqrt5 + 1) / 2)
  probs <- c((sqrt5 + 1) / (2 * sqrt5),
             (sqrt5 - 1) / (2 * sqrt5))
  sample(values, n, replace = TRUE, prob = probs)
}


.frd_side_weights <- function(x, x0, h, p, kersel, Q,
                              weights = c("irbc", "lp"),
                              bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)
  weights <- match.arg(weights)
  Q <- as.numeric(Q)
  n <- length(x)
  if (n == 0) {
    return(list(active = integer(0), omega = numeric(0)))
  }

  if (weights == "irbc") {
    active <- which(abs((x - x0) / (2 * h)) <= 1)
    if (length(active) == 0) {
      return(list(active = integer(0), omega = numeric(0)))
    }
    W <- .wiRBC_vec(x = x, xi = x[active], x0 = x0, h = h,
                    p = p, kersel = kersel, Q = Q,
                    bwcheck = bwcheck)
  } else {
    active <- which(abs((x - x0) / h) <= 1)
    if (length(active) == 0) {
      return(list(active = integer(0), omega = numeric(0)))
    }
    W <- Q * wix0(x = x, xi = x[active], x0 = x0, h = h,
                  p = p, v = 0, kersel = kersel,
                  bwcheck = bwcheck)
  }

  list(active = active, omega = as.numeric(W) / (n * h))
}


.irbcrd_joint_sandwich <- function(x.plus, x.minus, x0,
                                   h.plus, h.minus, p,
                                   kersel.plus, kersel.minus,
                                   Q.plus, Q.minus,
                                   eps.plus, eps.minus,
                                   weights = c("irbc", "lp"),
                                   component.names = NULL,
                                   bwcheck = .bwcheck_default) {
  weights <- match.arg(weights)
  eps.plus <- as.matrix(eps.plus)
  eps.minus <- as.matrix(eps.minus)

  if (nrow(eps.plus) != length(x.plus) ||
      nrow(eps.minus) != length(x.minus)) {
    stop("Residual matrices must have one row per observation on their respective cutoff side.")
  }
  if (ncol(eps.plus) < 1L || ncol(eps.plus) != ncol(eps.minus)) {
    stop("Residual matrices must have the same positive number of columns.")
  }
  if (any(!is.finite(eps.plus)) || any(!is.finite(eps.minus))) {
    stop("Residual matrices must contain only finite values.")
  }

  d <- ncol(eps.plus)
  if (is.null(component.names)) {
    component.names <- colnames(eps.plus)
  }
  if (is.null(component.names)) {
    component.names <- paste0("component", seq_len(d))
  }
  if (length(component.names) != d) {
    stop("'component.names' must have one entry per residual column.")
  }

  cov.mat <- matrix(0, nrow = d, ncol = d)

  add_side <- function(x, h, kersel, Q, eps) {
    side.w <- .frd_side_weights(x = x, x0 = x0, h = h, p = p,
                                kersel = kersel, Q = Q,
                                weights = weights,
                                bwcheck = bwcheck)
    if (length(side.w$active) == 0) {
      return(matrix(0, nrow = d, ncol = d))
    }
    eps.mat <- eps[side.w$active, , drop = FALSE]
    crossprod(eps.mat * side.w$omega)
  }

  cov.mat <- cov.mat +
    add_side(x = x.plus, h = h.plus, kersel = kersel.plus,
             Q = Q.plus, eps = eps.plus)
  cov.mat <- cov.mat +
    add_side(x = x.minus, h = h.minus, kersel = kersel.minus,
             Q = Q.minus, eps = eps.minus)

  dimnames(cov.mat) <- list(component.names, component.names)
  cov.mat
}


.irbcrd_linear_se <- function(cov.mat, weights, design = "Linear-combination") {
  cov.mat <- as.matrix(cov.mat)
  weights <- as.numeric(weights)
  if (nrow(cov.mat) != ncol(cov.mat) ||
      nrow(cov.mat) != length(weights)) {
    stop("The covariance matrix and linear-combination weights are not conformable.")
  }
  v <- as.numeric(t(weights) %*% cov.mat %*% weights)
  if (v < 0 && v > -1e-12) v <- 0
  if (!is.finite(v) || v < 0) {
    stop(design, " standard error could not be computed.")
  }
  sqrt(v)
}


.frd_joint_sandwich <- function(x.plus, x.minus, x0, h.plus, h.minus, p,
                                kersel.plus, kersel.minus, Q.plus, Q.minus,
                                eps.y.plus, eps.y.minus,
                                eps.d.plus, eps.d.minus,
                                weights = c("irbc", "lp"),
                                bwcheck = .bwcheck_default) {
  .irbcrd_joint_sandwich(
    x.plus = x.plus, x.minus = x.minus, x0 = x0,
    h.plus = h.plus, h.minus = h.minus, p = p,
    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
    Q.plus = Q.plus, Q.minus = Q.minus,
    eps.plus = cbind(outcome = as.numeric(eps.y.plus),
                     treatment = as.numeric(eps.d.plus)),
    eps.minus = cbind(outcome = as.numeric(eps.y.minus),
                      treatment = as.numeric(eps.d.minus)),
    weights = weights,
    component.names = c("outcome", "treatment"),
    bwcheck = bwcheck
  )
}


.frd_linear_se <- function(cov.mat, weights) {
  .irbcrd_linear_se(cov.mat = cov.mat, weights = weights,
                    design = "Fuzzy RD")
}


iRBCfrd <- function(y, d, x, h, c = NULL, p = NULL, kersel = NULL, res = NULL,
                    alpha = NULL, fast = NULL, g.loo = NULL,
                    boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                    boot.wild = NULL, boot.keep = NULL,
                    bwcheck = .bwcheck_default,
                    h.plus = NULL, h.minus = NULL,
                    kersel.plus = NULL, kersel.minus = NULL,
                    q.method = NULL,
                    fuzzy.boot.weights = NULL) {

  if (!is.numeric(y) || !is.numeric(d) || !is.numeric(x)) {
    stop("'y', 'd', and 'x' must be numeric vectors")
  }
  if (length(y) != length(x) || length(d) != length(x)) {
    stop("'y', 'd', and 'x' must have the same length")
  }
  if (any(!is.finite(y)) || any(!is.finite(d)) || any(!is.finite(x))) {
    stop("'y', 'd', and 'x' must contain only finite values")
  }

  if (is.null(c)) c <- 0
  if (is.null(p)) p <- 1
  if (is.null(alpha)) alpha <- 0.05
  if (is.null(fast)) fast <- FALSE
  if (is.null(g.loo)) g.loo <- FALSE
  if (is.null(boot.ci)) boot.ci <- FALSE
  if (is.null(boot.keep)) boot.keep <- TRUE
  bwcheck <- .validate_bwcheck(bwcheck)
  boot.wild <- .irbcrd_normalize_boot_wild(boot.wild)
  fuzzy.boot.weights <- .irbcrd_normalize_fuzzy_boot_weights(fuzzy.boot.weights)

  if (is.null(boot.B)) {
    boot.B <- if (boot.ci) 999 else 0
  } else if (!is.numeric(boot.B) || length(boot.B) > 1 || boot.B < 1 || boot.B != floor(boot.B)) {
    stop("'boot.B' must be a positive integer")
  }
  boot.B <- as.integer(boot.B)

  sharp.y <- iRBCrd(y = y, x = x, h = h, c = c, p = p, kersel = kersel,
                    res = res, alpha = alpha, fast = fast, g.loo = g.loo,
                    boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
                    boot.wild = boot.wild, boot.keep = FALSE,
                    bwcheck = bwcheck,
                    h.plus = h.plus, h.minus = h.minus,
                    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
                    q.method = q.method)
  sharp.d <- iRBCrd(y = d, x = x, h = h, c = c, p = p, kersel = kersel,
                    res = res, alpha = alpha, fast = fast, g.loo = g.loo,
                    boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
                    boot.wild = boot.wild, boot.keep = FALSE,
                    bwcheck = bwcheck,
                    h.plus = h.plus, h.minus = h.minus,
                    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
                    q.method = q.method)

  tau.y <- as.numeric(sharp.y$tau.hat)
  tau.d <- as.numeric(sharp.d$tau.hat)
  if (!is.finite(tau.d) || abs(tau.d) <= sqrt(.Machine$double.eps)) {
    stop("The estimated first-stage discontinuity is too close to zero for fuzzy RD inference.")
  }

  frd.weights <- c(1 / tau.d, -tau.y / (tau.d^2))
  zeta.hat <- tau.y / tau.d

  Bhat.y <- tau.y - as.numeric(sharp.y$tau.hat_irbc)
  Bhat.d <- tau.d - as.numeric(sharp.d$tau.hat_irbc)
  Bhat.frd <- as.numeric(c(Bhat.y, Bhat.d) %*% frd.weights)
  zeta.hat_irbc <- zeta.hat - Bhat.frd

  indp <- x >= c
  indm <- x < c
  xp <- x[indp]
  xm <- x[indm]

  cov.irbc <- .frd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.y.plus = sharp.y$epshat.plus,
    eps.y.minus = sharp.y$epshat.minus,
    eps.d.plus = sharp.d$epshat.plus,
    eps.d.minus = sharp.d$epshat.minus,
    weights = "irbc",
    bwcheck = bwcheck)
  se_irbc <- .frd_linear_se(cov.irbc, frd.weights)

  ci_irbc <- c(zeta.hat_irbc - qnorm(1 - alpha / 2) * se_irbc,
               zeta.hat_irbc - qnorm(alpha / 2) * se_irbc)

  cov.us <- .frd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.y.plus = sharp.y$epshat.plus,
    eps.y.minus = sharp.y$epshat.minus,
    eps.d.plus = sharp.d$epshat.plus,
    eps.d.minus = sharp.d$epshat.minus,
    weights = "lp",
    bwcheck = bwcheck)
  se_us <- .frd_linear_se(cov.us, frd.weights)

  out <- list(n = length(x), h = sharp.y$h, h.requested = sharp.y$h.requested,
              bwcheck = bwcheck, h.plus = sharp.y$h.plus,
              h.minus = sharp.y$h.minus, kersel = sharp.y$kersel,
              kersel.plus = sharp.y$kersel.plus,
              kersel.minus = sharp.y$kersel.minus,
              res = sharp.y$res, q.method = sharp.y$q.method,
              Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
              tau.y.hat = tau.y, tau.d.hat = tau.d,
              tau.y.hat_irbc = as.numeric(sharp.y$tau.hat_irbc),
              tau.d.hat_irbc = as.numeric(sharp.d$tau.hat_irbc),
              zeta.hat = zeta.hat, zeta.hat_irbc = zeta.hat_irbc,
              estimate_irbc = as.numeric(zeta.hat_irbc),
              Bhat.y = Bhat.y, Bhat.d = Bhat.d, Bhat.frd = Bhat.frd,
              frd.weights = frd.weights, cov.irbc = cov.irbc,
              cov.us = cov.us, se_irbc = se_irbc, se_us = se_us,
              ci_irbc = ci_irbc, sharp.outcome = sharp.y,
              sharp.treatment = sharp.d)

  if (boot.ci == TRUE) {
    if (!is.null(boot.seed)) set.seed(as.integer(boot.seed))

    eps.y <- numeric(length(x))
    eps.d <- numeric(length(x))
    eps.y[indp] <- as.numeric(sharp.y$epshat.plus)
    eps.y[indm] <- as.numeric(sharp.y$epshat.minus)
    eps.d[indp] <- as.numeric(sharp.d$epshat.plus)
    eps.d[indm] <- as.numeric(sharp.d$epshat.minus)

    boot.local <- (indp & abs(x - c) <= sharp.y$h.plus) |
      (indm & abs(x - c) <= sharp.y$h.minus)
    if (!any(boot.local & indp) || !any(boot.local & indm)) {
      stop("Bootstrap CI requires observations on both sides of the cutoff inside the bandwidth.")
    }

    ghat.y.boot.mean <- numeric(length(x))
    ghat.d.boot.mean <- numeric(length(x))
    boot.local.indices <- which(boot.local)
    yp <- y[indp]
    ym <- y[indm]
    dp <- d[indp]
    dm <- d[indm]
    ghatp.y <- LP.est(y = yp, x = xp, x0 = c, h = sharp.y$h.plus,
                      p = p, v = 0, kersel = sharp.y$kersel.plus,
                      bwcheck = bwcheck)
    ghatm.y <- LP.est(y = ym, x = xm, x0 = c, h = sharp.y$h.minus,
                      p = p, v = 0, kersel = sharp.y$kersel.minus,
                      bwcheck = bwcheck)
    ghatp.d <- LP.est(y = dp, x = xp, x0 = c, h = sharp.y$h.plus,
                      p = p, v = 0, kersel = sharp.y$kersel.plus,
                      bwcheck = bwcheck)
    ghatm.d <- LP.est(y = dm, x = xm, x0 = c, h = sharp.y$h.minus,
                      p = p, v = 0, kersel = sharp.y$kersel.minus,
                      bwcheck = bwcheck)

    for (ix in seq_along(boot.local.indices)) {
      i <- boot.local.indices[ix]
      if (x[i] >= c) {
        ghat.y.boot.mean[i] <- LP.est(y = yp, x = xp, x0 = x[i],
                                      h = sharp.y$h.plus, p = p, v = 0,
                                      kersel = sharp.y$kersel.plus,
                                      bwcheck = bwcheck)
        ghat.d.boot.mean[i] <- LP.est(y = dp, x = xp, x0 = x[i],
                                      h = sharp.y$h.plus, p = p, v = 0,
                                      kersel = sharp.y$kersel.plus,
                                      bwcheck = bwcheck)
      } else {
        ghat.y.boot.mean[i] <- LP.est(y = ym, x = xm, x0 = x[i],
                                      h = sharp.y$h.minus, p = p, v = 0,
                                      kersel = sharp.y$kersel.minus,
                                      bwcheck = bwcheck)
        ghat.d.boot.mean[i] <- LP.est(y = dm, x = xm, x0 = x[i],
                                      h = sharp.y$h.minus, p = p, v = 0,
                                      kersel = sharp.y$kersel.minus,
                                      bwcheck = bwcheck)
      }
    }

    boot.t <- rep(NA_real_, boot.B)
    for (bb in seq_len(boot.B)) {
      e.boot <- .irbcrd_wild_draws(length(boot.local.indices), boot.wild)

      y.boot <- y
      d.boot <- d
      y.boot[boot.local.indices] <- ghat.y.boot.mean[boot.local.indices] +
        eps.y[boot.local.indices] * e.boot
      d.boot[boot.local.indices] <- ghat.d.boot.mean[boot.local.indices] +
        eps.d[boot.local.indices] * e.boot

      ghatp.y.boot <- LP.est(y = y.boot[indp], x = xp, x0 = c,
                             h = sharp.y$h.plus, p = p, v = 0,
                             kersel = sharp.y$kersel.plus,
                             bwcheck = bwcheck)
      ghatm.y.boot <- LP.est(y = y.boot[indm], x = xm, x0 = c,
                             h = sharp.y$h.minus, p = p, v = 0,
                             kersel = sharp.y$kersel.minus,
                             bwcheck = bwcheck)
      ghatp.d.boot <- LP.est(y = d.boot[indp], x = xp, x0 = c,
                             h = sharp.y$h.plus, p = p, v = 0,
                             kersel = sharp.y$kersel.plus,
                             bwcheck = bwcheck)
      ghatm.d.boot <- LP.est(y = d.boot[indm], x = xm, x0 = c,
                             h = sharp.y$h.minus, p = p, v = 0,
                             kersel = sharp.y$kersel.minus,
                             bwcheck = bwcheck)

      T.rd.y <- as.numeric(sharp.y$Q.plus * (ghatp.y.boot - ghatp.y) -
                             sharp.y$Q.minus * (ghatm.y.boot - ghatm.y))
      T.rd.d <- as.numeric(sharp.y$Q.plus * (ghatp.d.boot - ghatp.d) -
                             sharp.y$Q.minus * (ghatm.d.boot - ghatm.d))
      if (fuzzy.boot.weights == "fixed") {
        boot.t[bb] <- as.numeric(c(T.rd.y, T.rd.d) %*% frd.weights)
      } else {
        tau.d.boot <- as.numeric(ghatp.d.boot - ghatm.d.boot)
        if (is.finite(tau.d.boot) &&
            abs(tau.d.boot) > sqrt(.Machine$double.eps)) {
          frd.weights.boot <- c(1 / tau.d.boot,
                                -tau.y / (tau.d * tau.d.boot))
          boot.t[bb] <- as.numeric(c(T.rd.y, T.rd.d) %*% frd.weights.boot)
        } else {
          boot.t[bb] <- NA_real_
        }
      }
    }

    boot.t.ok <- boot.t[is.finite(boot.t)]
    boot.success <- length(boot.t.ok)
    if (!is.finite(se_us) || se_us <= 0) {
      stop("Fuzzy RD bootstrap calibration cannot be computed because the non-bias-corrected standard error is not positive.")
    }
    mhat <- as.numeric(se_irbc / se_us)
    H.inv <- function(u) pnorm(mhat * qnorm(u))
    boot.probs <- H.inv(c(alpha / 2, 1 - alpha / 2))

    if (boot.success < 2) {
      warning("Fewer than two successful bootstrap repetitions; bootstrap CI is NA.")
      ci_irbc_boot <- c(NA_real_, NA_real_)
    } else {
      q.boot <- as.numeric(quantile(boot.t.ok, probs = boot.probs,
                                    names = FALSE, type = 7))
      ci_irbc_boot <- c(as.numeric(zeta.hat) - q.boot[2],
                        as.numeric(zeta.hat) - q.boot[1])
    }

    out <- c(out, list(estimate_irbc_boot = as.numeric(zeta.hat_irbc),
                       ci_irbc_boot = ci_irbc_boot,
                       boot.B = boot.B,
                       boot.success = boot.success,
                       boot.probs = boot.probs,
                       boot.mhat = mhat,
                       boot.wild = boot.wild,
                       fuzzy.boot.weights = fuzzy.boot.weights))
    if (boot.keep == TRUE) out$boot.t <- boot.t
  }

  return(out)
}


.irbcrd_validate_covs <- function(covs, n) {
  if (is.data.frame(covs)) {
    numeric.columns <- vapply(covs, is.numeric, logical(1))
    if (!all(numeric.columns)) {
      stop("All columns of 'covs' must be numeric.")
    }
    covs <- as.matrix(covs)
  } else if (is.numeric(covs) && is.null(dim(covs))) {
    covs <- matrix(covs, ncol = 1L)
  } else if (!is.matrix(covs) || !is.numeric(covs)) {
    stop("'covs' must be a numeric vector, matrix, or data frame.")
  }

  if (nrow(covs) != n) {
    stop("'covs' must have one row per observation in 'y' and 'x'.")
  }
  if (ncol(covs) < 1L) {
    stop("'covs' must contain at least one covariate.")
  }
  if (any(!is.finite(covs))) {
    stop("'covs' must contain only finite values.")
  }

  storage.mode(covs) <- "double"
  covariate.names <- colnames(covs)
  if (is.null(covariate.names)) {
    covariate.names <- paste0("z", seq_len(ncol(covs)))
  } else {
    empty <- !nzchar(covariate.names)
    covariate.names[empty] <- paste0("z", which(empty))
    covariate.names <- make.unique(covariate.names)
  }
  colnames(covs) <- covariate.names
  covs
}


.irbcrd_common_covariate_slope <- function(y, covs, x, x0,
                                           h.plus, h.minus, p,
                                           kersel.plus, kersel.minus) {
  d <- ncol(covs)

  side_terms <- function(index, h, kersel) {
    x.side <- x[index]
    y.side <- y[index]
    z.side <- covs[index, , drop = FALSE]
    u <- (x.side - x0) / h
    w <- K(u = u, kersel = kersel) / h
    active <- which(w > 0)
    if (length(active) <= p + 1L) {
      stop("The common-slope covariate regression requires more than p + 1 positive-weight observations on each side.")
    }

    R <- .poly_basis(u[active], p)
    w <- w[active]
    y.side <- y.side[active]
    z.side <- z.side[active, , drop = FALSE]

    A <- crossprod(R * sqrt(w))
    A.inv <- try(chol2inv(chol(A)), silent = TRUE)
    if (inherits(A.inv, "try-error")) {
      stop("The local-polynomial terms are singular in the common-slope covariate regression.")
    }

    U.z <- crossprod(R * w, z.side)
    U.y <- crossprod(R * w, y.side)
    ZWZ <- crossprod(z.side * w, z.side) -
      crossprod(U.z, A.inv %*% U.z)
    ZWY <- crossprod(z.side * w, y.side) -
      crossprod(U.z, A.inv %*% U.y)

    list(ZWZ = ZWZ, ZWY = ZWY, n.active = length(active))
  }

  plus <- side_terms(x >= x0, h = h.plus, kersel = kersel.plus)
  minus <- side_terms(x < x0, h = h.minus, kersel = kersel.minus)
  ZWZ <- plus$ZWZ + minus$ZWZ
  ZWY <- plus$ZWY + minus$ZWY
  ZWZ <- (ZWZ + t(ZWZ)) / 2

  ZWZ.inv <- try(chol2inv(chol(ZWZ)), silent = TRUE)
  if (inherits(ZWZ.inv, "try-error")) {
    stop("The covariates are locally collinear after partialling out the side-specific polynomial terms.")
  }

  gamma.hat <- as.numeric(ZWZ.inv %*% ZWY)
  names(gamma.hat) <- colnames(covs)
  if (length(gamma.hat) != d || any(!is.finite(gamma.hat))) {
    stop("The common covariate coefficient could not be computed.")
  }

  list(gamma.hat = gamma.hat, ZWZ = ZWZ, ZWY = ZWY,
       n.active = c(minus = minus$n.active, plus = plus$n.active))
}


.irbcrd_lp_component_weights <- function(x, x0, h, p, kersel,
                                         bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)
  as.numeric(wix0(x = x, xi = x, x0 = x0, h = h, p = p, v = 0,
                  kersel = kersel, bwcheck = bwcheck)) /
    (length(x) * h)
}


.irbcrd_lp_component_estimate <- function(values, x, x0, h, p, kersel,
                                          bwcheck = .bwcheck_default) {
  values <- as.matrix(values)
  if (nrow(values) != length(x)) {
    stop("Component data must have one row per running-variable observation.")
  }
  omega <- .irbcrd_lp_component_weights(
    x = x, x0 = x0, h = h, p = p, kersel = kersel,
    bwcheck = bwcheck
  )
  as.numeric(crossprod(omega, values))
}


## Covariate-adjusted sharp RD with a CCF&T common covariate coefficient.
iRBCcrd <- function(y, covs, x, h, c = NULL, p = NULL,
                    kersel = NULL, res = NULL,
                    alpha = NULL, fast = NULL, g.loo = NULL,
                    boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                    boot.wild = NULL, boot.keep = NULL,
                    bwcheck = .bwcheck_default,
                    h.plus = NULL, h.minus = NULL,
                    kersel.plus = NULL, kersel.minus = NULL,
                    q.method = NULL) {

  if (!is.numeric(y) || !is.numeric(x)) {
    stop("'y' and 'x' must be numeric vectors.")
  }
  if (length(y) != length(x)) {
    stop("'y' and 'x' must have the same length.")
  }
  if (length(y) < 1L || any(!is.finite(y)) || any(!is.finite(x))) {
    stop("'y' and 'x' must be nonempty and contain only finite values.")
  }
  covs <- .irbcrd_validate_covs(covs, n = length(x))

  if (is.null(c)) c <- 0
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) {
    stop("'c' must be a finite numeric scalar.")
  }
  if (!any(x < c) || !any(x >= c)) {
    stop("Observations are required on both sides of the cutoff.")
  }
  if (is.null(p)) p <- 1
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      p < 0 || p != floor(p)) {
    stop("'p' must be a non-negative integer.")
  }
  p <- as.integer(p)
  if (is.null(alpha)) alpha <- 0.05
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be strictly between zero and one.")
  }
  if (is.null(fast)) fast <- FALSE
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("'fast' must be a single TRUE/FALSE.")
  }
  if (is.null(g.loo)) g.loo <- FALSE
  if (!is.logical(g.loo) || length(g.loo) != 1L || is.na(g.loo)) {
    stop("'g.loo' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.ci)) boot.ci <- FALSE
  if (!is.logical(boot.ci) || length(boot.ci) != 1L || is.na(boot.ci)) {
    stop("'boot.ci' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.keep)) boot.keep <- TRUE
  if (!is.logical(boot.keep) || length(boot.keep) != 1L || is.na(boot.keep)) {
    stop("'boot.keep' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.B)) {
    boot.B <- if (boot.ci) 999L else 0L
  } else if (!is.numeric(boot.B) || length(boot.B) != 1L ||
             !is.finite(boot.B) || boot.B < 1 || boot.B != floor(boot.B)) {
    stop("'boot.B' must be a positive integer.")
  }
  boot.B <- as.integer(boot.B)
  if (!is.null(boot.seed) &&
      (!is.numeric(boot.seed) || length(boot.seed) != 1L ||
       !is.finite(boot.seed))) {
    stop("'boot.seed' must be a finite numeric scalar.")
  }

  bwcheck <- .validate_bwcheck(bwcheck)
  boot.wild <- .irbcrd_normalize_boot_wild(boot.wild)

  sharp.y <- iRBCrd(
    y = y, x = x, h = h, c = c, p = p, kersel = kersel,
    res = res, alpha = alpha, fast = fast, g.loo = g.loo,
    boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
    boot.wild = boot.wild, boot.keep = FALSE,
    bwcheck = bwcheck,
    h.plus = h.plus, h.minus = h.minus,
    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
    q.method = q.method
  )

  covariate.names <- colnames(covs)
  sharp.covs <- lapply(seq_len(ncol(covs)), function(jj) {
    iRBCrd(
      y = covs[, jj], x = x, h = h, c = c, p = p, kersel = kersel,
      res = sharp.y$res, alpha = alpha, fast = fast, g.loo = g.loo,
      boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
      boot.wild = boot.wild, boot.keep = FALSE,
      bwcheck = bwcheck,
      h.plus = h.plus, h.minus = h.minus,
      kersel.plus = kersel.plus, kersel.minus = kersel.minus,
      q.method = sharp.y$q.method
    )
  })
  names(sharp.covs) <- covariate.names

  same_design <- vapply(sharp.covs, function(fit) {
    isTRUE(all.equal(fit$h.plus, sharp.y$h.plus, tolerance = 1e-12)) &&
      isTRUE(all.equal(fit$h.minus, sharp.y$h.minus, tolerance = 1e-12)) &&
      identical(fit$kersel.plus, sharp.y$kersel.plus) &&
      identical(fit$kersel.minus, sharp.y$kersel.minus) &&
      isTRUE(all.equal(fit$Q.plus, sharp.y$Q.plus, tolerance = 1e-12)) &&
      isTRUE(all.equal(fit$Q.minus, sharp.y$Q.minus, tolerance = 1e-12))
  }, logical(1))
  if (!all(same_design)) {
    stop("Outcome and covariate components did not use identical effective RD specifications.")
  }

  gamma.fit <- .irbcrd_common_covariate_slope(
    y = y, covs = covs, x = x, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus
  )
  gamma.hat <- gamma.fit$gamma.hat

  component.names <- c("outcome", covariate.names)
  tau.components <- c(
    outcome = as.numeric(sharp.y$tau.hat),
    vapply(sharp.covs, function(fit) as.numeric(fit$tau.hat), numeric(1))
  )
  tau.components.irbc <- c(
    outcome = as.numeric(sharp.y$tau.hat_irbc),
    vapply(sharp.covs,
           function(fit) as.numeric(fit$tau.hat_irbc), numeric(1))
  )
  names(tau.components) <- component.names
  names(tau.components.irbc) <- component.names

  crd.weights <- c(outcome = 1, -gamma.hat)
  names(crd.weights) <- component.names
  Bhat.components <- tau.components - tau.components.irbc
  tau.hat <- as.numeric(crossprod(crd.weights, tau.components))
  Bhat.crd <- as.numeric(crossprod(crd.weights, Bhat.components))
  tau.hat_irbc <- tau.hat - Bhat.crd

  indp <- x >= c
  indm <- x < c
  xp <- x[indp]
  xm <- x[indm]
  n.components <- length(component.names)
  eps.plus <- matrix(NA_real_, nrow = length(xp), ncol = n.components,
                     dimnames = list(NULL, component.names))
  eps.minus <- matrix(NA_real_, nrow = length(xm), ncol = n.components,
                      dimnames = list(NULL, component.names))
  eps.plus[, 1L] <- as.numeric(sharp.y$epshat.plus)
  eps.minus[, 1L] <- as.numeric(sharp.y$epshat.minus)
  for (jj in seq_along(sharp.covs)) {
    eps.plus[, jj + 1L] <- as.numeric(sharp.covs[[jj]]$epshat.plus)
    eps.minus[, jj + 1L] <- as.numeric(sharp.covs[[jj]]$epshat.minus)
  }

  cov.irbc <- .irbcrd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.plus = eps.plus, eps.minus = eps.minus,
    weights = "irbc", component.names = component.names,
    bwcheck = bwcheck
  )
  se_irbc <- .irbcrd_linear_se(
    cov.mat = cov.irbc, weights = crd.weights,
    design = "Covariate-adjusted RD"
  )

  cov.us <- .irbcrd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.plus = eps.plus, eps.minus = eps.minus,
    weights = "lp", component.names = component.names,
    bwcheck = bwcheck
  )
  se_us <- .irbcrd_linear_se(
    cov.mat = cov.us, weights = crd.weights,
    design = "Covariate-adjusted RD bootstrap calibration"
  )

  ci_irbc <- c(
    tau.hat_irbc - qnorm(1 - alpha / 2) * se_irbc,
    tau.hat_irbc - qnorm(alpha / 2) * se_irbc
  )

  out <- list(
    n = length(x), h = sharp.y$h, h.requested = sharp.y$h.requested,
    bwcheck = bwcheck, h.plus = sharp.y$h.plus,
    h.minus = sharp.y$h.minus, kersel = sharp.y$kersel,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    res = sharp.y$res, q.method = sharp.y$q.method,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    covariate.names = covariate.names,
    gamma.hat = gamma.hat, beta.covs = gamma.hat,
    gamma.normal.matrix = gamma.fit$ZWZ,
    gamma.normal.vector = gamma.fit$ZWY,
    gamma.n = gamma.fit$n.active,
    tau.y.hat = tau.components[["outcome"]],
    tau.z.hat = tau.components[-1L],
    tau.y.hat_irbc = tau.components.irbc[["outcome"]],
    tau.z.hat_irbc = tau.components.irbc[-1L],
    tau.hat = tau.hat, tau.hat_irbc = tau.hat_irbc,
    estimate_irbc = as.numeric(tau.hat_irbc),
    Bhat.y = Bhat.components[["outcome"]],
    Bhat.z = Bhat.components[-1L],
    Bhat.crd = Bhat.crd,
    crd.weights = crd.weights,
    component.estimates = tau.components,
    component.estimates.irbc = tau.components.irbc,
    component.bias = Bhat.components,
    cov.irbc = cov.irbc, cov.us = cov.us,
    se_irbc = se_irbc, se_us = se_us, ci_irbc = ci_irbc,
    epshat.plus = eps.plus, epshat.minus = eps.minus,
    sharp.outcome = sharp.y, sharp.covariates = sharp.covs
  )

  if (boot.ci == TRUE) {
    if (!is.null(boot.seed)) set.seed(as.integer(boot.seed))
    if (!is.finite(se_us) || se_us <= 0) {
      stop("Covariate-adjusted RD bootstrap calibration requires a positive non-bias-corrected standard error.")
    }
    mhat <- as.numeric(se_irbc / se_us)
    if (!is.finite(mhat) || mhat <= 0) {
      stop("Covariate-adjusted RD bootstrap calibration could not be computed.")
    }
    H.inv <- function(u) pnorm(mhat * qnorm(u))
    boot.probs <- H.inv(c(alpha / 2, 1 - alpha / 2))

    component.data <- cbind(outcome = as.numeric(y), covs)
    eps.full <- matrix(0, nrow = length(x), ncol = n.components,
                       dimnames = list(NULL, component.names))
    eps.full[indp, ] <- eps.plus
    eps.full[indm, ] <- eps.minus

    boot.local <- (indp & abs(x - c) <= sharp.y$h.plus) |
      (indm & abs(x - c) <= sharp.y$h.minus)
    if (!any(boot.local & indp) || !any(boot.local & indm)) {
      stop("Bootstrap CI requires observations on both sides of the cutoff inside the bandwidth.")
    }
    boot.local.indices <- which(boot.local)
    ghat.boot.mean <- matrix(0, nrow = length(x), ncol = n.components,
                             dimnames = list(NULL, component.names))
    for (ii in boot.local.indices) {
      if (x[ii] >= c) {
        ghat.boot.mean[ii, ] <- .irbcrd_lp_component_estimate(
          values = component.data[indp, , drop = FALSE],
          x = xp, x0 = x[ii], h = sharp.y$h.plus, p = p,
          kersel = sharp.y$kersel.plus, bwcheck = bwcheck
        )
      } else {
        ghat.boot.mean[ii, ] <- .irbcrd_lp_component_estimate(
          values = component.data[indm, , drop = FALSE],
          x = xm, x0 = x[ii], h = sharp.y$h.minus, p = p,
          kersel = sharp.y$kersel.minus, bwcheck = bwcheck
        )
      }
    }

    omega.plus <- .irbcrd_lp_component_weights(
      x = xp, x0 = c, h = sharp.y$h.plus, p = p,
      kersel = sharp.y$kersel.plus, bwcheck = bwcheck
    )
    omega.minus <- .irbcrd_lp_component_weights(
      x = xm, x0 = c, h = sharp.y$h.minus, p = p,
      kersel = sharp.y$kersel.minus, bwcheck = bwcheck
    )
    ghat.plus <- as.numeric(crossprod(
      omega.plus, component.data[indp, , drop = FALSE]
    ))
    ghat.minus <- as.numeric(crossprod(
      omega.minus, component.data[indm, , drop = FALSE]
    ))

    boot.t <- rep(NA_real_, boot.B)
    boot.t.components <- matrix(
      NA_real_, nrow = boot.B, ncol = n.components,
      dimnames = list(NULL, component.names)
    )
    for (bb in seq_len(boot.B)) {
      e.boot <- .irbcrd_wild_draws(length(boot.local.indices), boot.wild)
      component.boot <- component.data
      component.boot[boot.local.indices, ] <-
        ghat.boot.mean[boot.local.indices, , drop = FALSE] +
        sweep(eps.full[boot.local.indices, , drop = FALSE],
              MARGIN = 1L, STATS = e.boot, FUN = "*")

      ghat.plus.boot <- as.numeric(crossprod(
        omega.plus, component.boot[indp, , drop = FALSE]
      ))
      ghat.minus.boot <- as.numeric(crossprod(
        omega.minus, component.boot[indm, , drop = FALSE]
      ))
      T.components <- as.numeric(
        as.numeric(sharp.y$Q.plus) * (ghat.plus.boot - ghat.plus) -
          as.numeric(sharp.y$Q.minus) * (ghat.minus.boot - ghat.minus)
      )
      boot.t.components[bb, ] <- T.components
      boot.t[bb] <- as.numeric(crossprod(crd.weights, T.components))
    }

    boot.t.ok <- boot.t[is.finite(boot.t)]
    boot.success <- length(boot.t.ok)
    if (boot.success < 2L) {
      warning("Fewer than two successful bootstrap repetitions; bootstrap CI is NA.")
      ci_irbc_boot <- c(NA_real_, NA_real_)
    } else {
      q.boot <- as.numeric(quantile(
        boot.t.ok, probs = boot.probs, names = FALSE, type = 7
      ))
      ci_irbc_boot <- c(tau.hat - q.boot[2L],
                        tau.hat - q.boot[1L])
    }

    out <- c(out, list(
      estimate_irbc_boot = as.numeric(tau.hat_irbc),
      ci_irbc_boot = ci_irbc_boot,
      boot.B = boot.B, boot.success = boot.success,
      boot.impl = "plug-in", boot.probs = boot.probs,
      boot.mhat = mhat, boot.wild = boot.wild,
      bootstrap.weights = "fixed-gamma"
    ))
    if (boot.keep == TRUE) {
      out$boot.t <- boot.t
      out$boot.t.components <- boot.t.components
    }
  }

  out
}


## Covariate-adjusted fuzzy RD with separate common covariate coefficients
## in the reduced-form and first-stage equations.
iRBCcfrd <- function(y, d, covs, x, h, c = NULL, p = NULL,
                     kersel = NULL, res = NULL,
                     alpha = NULL, fast = NULL, g.loo = NULL,
                     boot.ci = NULL, boot.B = NULL, boot.seed = NULL,
                     boot.wild = NULL, boot.keep = NULL,
                     bwcheck = .bwcheck_default,
                     h.plus = NULL, h.minus = NULL,
                     kersel.plus = NULL, kersel.minus = NULL,
                     q.method = NULL,
                     fuzzy.boot.weights = NULL,
                     gamma.boot.weights = NULL) {

  if (!is.numeric(y) || !is.numeric(d) || !is.numeric(x)) {
    stop("'y', 'd', and 'x' must be numeric vectors.")
  }
  if (length(y) != length(x) || length(d) != length(x)) {
    stop("'y', 'd', and 'x' must have the same length.")
  }
  if (length(y) < 1L || any(!is.finite(y)) || any(!is.finite(d)) ||
      any(!is.finite(x))) {
    stop("'y', 'd', and 'x' must be nonempty and contain only finite values.")
  }
  covs <- .irbcrd_validate_covs(covs, n = length(x))
  covariate.names <- make.unique(
    c("outcome", "treatment", colnames(covs))
  )[-c(1L, 2L)]
  colnames(covs) <- covariate.names

  if (is.null(c)) c <- 0
  if (!is.numeric(c) || length(c) != 1L || !is.finite(c)) {
    stop("'c' must be a finite numeric scalar.")
  }
  if (!any(x < c) || !any(x >= c)) {
    stop("Observations are required on both sides of the cutoff.")
  }
  if (is.null(p)) p <- 1
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) ||
      p < 0 || p != floor(p)) {
    stop("'p' must be a non-negative integer.")
  }
  p <- as.integer(p)
  if (is.null(alpha)) alpha <- 0.05
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be strictly between zero and one.")
  }
  if (is.null(fast)) fast <- FALSE
  if (!is.logical(fast) || length(fast) != 1L || is.na(fast)) {
    stop("'fast' must be a single TRUE/FALSE.")
  }
  if (is.null(g.loo)) g.loo <- FALSE
  if (!is.logical(g.loo) || length(g.loo) != 1L || is.na(g.loo)) {
    stop("'g.loo' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.ci)) boot.ci <- FALSE
  if (!is.logical(boot.ci) || length(boot.ci) != 1L || is.na(boot.ci)) {
    stop("'boot.ci' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.keep)) boot.keep <- TRUE
  if (!is.logical(boot.keep) || length(boot.keep) != 1L || is.na(boot.keep)) {
    stop("'boot.keep' must be a single TRUE/FALSE.")
  }
  if (is.null(boot.B)) {
    boot.B <- if (boot.ci) 999L else 0L
  } else if (!is.numeric(boot.B) || length(boot.B) != 1L ||
             !is.finite(boot.B) || boot.B < 1 || boot.B != floor(boot.B)) {
    stop("'boot.B' must be a positive integer.")
  }
  boot.B <- as.integer(boot.B)
  if (!is.null(boot.seed) &&
      (!is.numeric(boot.seed) || length(boot.seed) != 1L ||
       !is.finite(boot.seed))) {
    stop("'boot.seed' must be a finite numeric scalar.")
  }

  bwcheck <- .validate_bwcheck(bwcheck)
  boot.wild <- .irbcrd_normalize_boot_wild(boot.wild)
  fuzzy.boot.weights <-
    .irbcrd_normalize_fuzzy_boot_weights(fuzzy.boot.weights)
  if (fuzzy.boot.weights != "fixed") {
    stop("'fuzzy.boot.weights' must be NULL or 'fixed' for covariate-adjusted fuzzy RD; bootstrap fuzzy weights are no longer supported for this design.")
  }
  gamma.boot.weights <-
    .irbcrd_normalize_gamma_boot_weights(gamma.boot.weights)

  sharp.y <- iRBCrd(
    y = y, x = x, h = h, c = c, p = p, kersel = kersel,
    res = res, alpha = alpha, fast = fast, g.loo = g.loo,
    boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
    boot.wild = boot.wild, boot.keep = FALSE,
    bwcheck = bwcheck,
    h.plus = h.plus, h.minus = h.minus,
    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
    q.method = q.method
  )
  sharp.d <- iRBCrd(
    y = d, x = x, h = h, c = c, p = p, kersel = kersel,
    res = sharp.y$res, alpha = alpha, fast = fast, g.loo = g.loo,
    boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
    boot.wild = boot.wild, boot.keep = FALSE,
    bwcheck = bwcheck,
    h.plus = h.plus, h.minus = h.minus,
    kersel.plus = kersel.plus, kersel.minus = kersel.minus,
    q.method = sharp.y$q.method
  )

  sharp.covs <- lapply(seq_len(ncol(covs)), function(jj) {
    iRBCrd(
      y = covs[, jj], x = x, h = h, c = c, p = p, kersel = kersel,
      res = sharp.y$res, alpha = alpha, fast = fast, g.loo = g.loo,
      boot.ci = FALSE, boot.B = NULL, boot.seed = NULL,
      boot.wild = boot.wild, boot.keep = FALSE,
      bwcheck = bwcheck,
      h.plus = h.plus, h.minus = h.minus,
      kersel.plus = kersel.plus, kersel.minus = kersel.minus,
      q.method = sharp.y$q.method
    )
  })
  names(sharp.covs) <- covariate.names

  component.fits <- c(list(outcome = sharp.y, treatment = sharp.d),
                      sharp.covs)
  same.design <- vapply(component.fits[-1L], function(fit) {
    isTRUE(all.equal(fit$h.plus, sharp.y$h.plus, tolerance = 1e-12)) &&
      isTRUE(all.equal(fit$h.minus, sharp.y$h.minus, tolerance = 1e-12)) &&
      identical(fit$kersel.plus, sharp.y$kersel.plus) &&
      identical(fit$kersel.minus, sharp.y$kersel.minus) &&
      isTRUE(all.equal(fit$Q.plus, sharp.y$Q.plus, tolerance = 1e-12)) &&
      isTRUE(all.equal(fit$Q.minus, sharp.y$Q.minus, tolerance = 1e-12))
  }, logical(1))
  if (!all(same.design)) {
    stop("Outcome, treatment, and covariate components did not use identical effective RD specifications.")
  }

  gamma.y.fit <- .irbcrd_common_covariate_slope(
    y = y, covs = covs, x = x, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus
  )
  gamma.d.fit <- .irbcrd_common_covariate_slope(
    y = d, covs = covs, x = x, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus
  )
  gamma.y.hat <- gamma.y.fit$gamma.hat
  gamma.d.hat <- gamma.d.fit$gamma.hat

  component.names <- c("outcome", "treatment", covariate.names)
  tau.components <- c(
    outcome = as.numeric(sharp.y$tau.hat),
    treatment = as.numeric(sharp.d$tau.hat),
    vapply(sharp.covs, function(fit) as.numeric(fit$tau.hat), numeric(1))
  )
  tau.components.irbc <- c(
    outcome = as.numeric(sharp.y$tau.hat_irbc),
    treatment = as.numeric(sharp.d$tau.hat_irbc),
    vapply(sharp.covs,
           function(fit) as.numeric(fit$tau.hat_irbc), numeric(1))
  )
  names(tau.components) <- component.names
  names(tau.components.irbc) <- component.names
  Bhat.components <- tau.components - tau.components.irbc

  tau.z.hat <- tau.components[covariate.names]
  tau.z.hat_irbc <- tau.components.irbc[covariate.names]
  Bhat.z <- Bhat.components[covariate.names]
  tau.y.adj.hat <- as.numeric(
    tau.components[["outcome"]] - crossprod(gamma.y.hat, tau.z.hat)
  )
  tau.d.adj.hat <- as.numeric(
    tau.components[["treatment"]] - crossprod(gamma.d.hat, tau.z.hat)
  )
  if (!is.finite(tau.d.adj.hat) ||
      abs(tau.d.adj.hat) <= sqrt(.Machine$double.eps)) {
    stop("The covariate-adjusted first-stage discontinuity is too close to zero for fuzzy RD inference.")
  }

  tau.y.adj.hat_irbc <- as.numeric(
    tau.components.irbc[["outcome"]] -
      crossprod(gamma.y.hat, tau.z.hat_irbc)
  )
  tau.d.adj.hat_irbc <- as.numeric(
    tau.components.irbc[["treatment"]] -
      crossprod(gamma.d.hat, tau.z.hat_irbc)
  )
  Bhat.y.adj <- as.numeric(
    Bhat.components[["outcome"]] - crossprod(gamma.y.hat, Bhat.z)
  )
  Bhat.d.adj <- as.numeric(
    Bhat.components[["treatment"]] - crossprod(gamma.d.hat, Bhat.z)
  )

  frd.weights <- c(
    reduced.form = 1 / tau.d.adj.hat,
    first.stage = -tau.y.adj.hat / (tau.d.adj.hat^2)
  )
  cfrd.weights <- c(
    outcome = frd.weights[["reduced.form"]],
    treatment = frd.weights[["first.stage"]],
    -(frd.weights[["reduced.form"]] * gamma.y.hat +
        frd.weights[["first.stage"]] * gamma.d.hat)
  )
  names(cfrd.weights) <- component.names

  zeta.hat <- tau.y.adj.hat / tau.d.adj.hat
  Bhat.cfrd <- as.numeric(crossprod(cfrd.weights, Bhat.components))
  zeta.hat_irbc <- zeta.hat - Bhat.cfrd

  indp <- x >= c
  indm <- x < c
  xp <- x[indp]
  xm <- x[indm]
  n.components <- length(component.names)
  eps.plus <- matrix(NA_real_, nrow = length(xp), ncol = n.components,
                     dimnames = list(NULL, component.names))
  eps.minus <- matrix(NA_real_, nrow = length(xm), ncol = n.components,
                      dimnames = list(NULL, component.names))
  eps.plus[, 1L] <- as.numeric(sharp.y$epshat.plus)
  eps.minus[, 1L] <- as.numeric(sharp.y$epshat.minus)
  eps.plus[, 2L] <- as.numeric(sharp.d$epshat.plus)
  eps.minus[, 2L] <- as.numeric(sharp.d$epshat.minus)
  for (jj in seq_along(sharp.covs)) {
    eps.plus[, jj + 2L] <- as.numeric(sharp.covs[[jj]]$epshat.plus)
    eps.minus[, jj + 2L] <- as.numeric(sharp.covs[[jj]]$epshat.minus)
  }

  cov.irbc <- .irbcrd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.plus = eps.plus, eps.minus = eps.minus,
    weights = "irbc", component.names = component.names,
    bwcheck = bwcheck
  )
  se_irbc <- .irbcrd_linear_se(
    cov.mat = cov.irbc, weights = cfrd.weights,
    design = "Covariate-adjusted fuzzy RD"
  )

  cov.us <- .irbcrd_joint_sandwich(
    x.plus = xp, x.minus = xm, x0 = c,
    h.plus = sharp.y$h.plus, h.minus = sharp.y$h.minus, p = p,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    eps.plus = eps.plus, eps.minus = eps.minus,
    weights = "lp", component.names = component.names,
    bwcheck = bwcheck
  )
  se_us <- .irbcrd_linear_se(
    cov.mat = cov.us, weights = cfrd.weights,
    design = "Covariate-adjusted fuzzy RD bootstrap calibration"
  )

  ci_irbc <- c(
    zeta.hat_irbc - qnorm(1 - alpha / 2) * se_irbc,
    zeta.hat_irbc - qnorm(alpha / 2) * se_irbc
  )

  out <- list(
    n = length(x), h = sharp.y$h, h.requested = sharp.y$h.requested,
    bwcheck = bwcheck, h.plus = sharp.y$h.plus,
    h.minus = sharp.y$h.minus, kersel = sharp.y$kersel,
    kersel.plus = sharp.y$kersel.plus,
    kersel.minus = sharp.y$kersel.minus,
    res = sharp.y$res, q.method = sharp.y$q.method,
    Q.plus = sharp.y$Q.plus, Q.minus = sharp.y$Q.minus,
    covariate.names = covariate.names,
    gamma.y.hat = gamma.y.hat, gamma.d.hat = gamma.d.hat,
    gamma.y.normal.matrix = gamma.y.fit$ZWZ,
    gamma.y.normal.vector = gamma.y.fit$ZWY,
    gamma.y.n = gamma.y.fit$n.active,
    gamma.d.normal.matrix = gamma.d.fit$ZWZ,
    gamma.d.normal.vector = gamma.d.fit$ZWY,
    gamma.d.n = gamma.d.fit$n.active,
    tau.y.hat = tau.components[["outcome"]],
    tau.d.hat = tau.components[["treatment"]],
    tau.z.hat = tau.z.hat,
    tau.y.hat_irbc = tau.components.irbc[["outcome"]],
    tau.d.hat_irbc = tau.components.irbc[["treatment"]],
    tau.z.hat_irbc = tau.z.hat_irbc,
    tau.y.adj.hat = tau.y.adj.hat,
    tau.d.adj.hat = tau.d.adj.hat,
    tau.y.adj.hat_irbc = tau.y.adj.hat_irbc,
    tau.d.adj.hat_irbc = tau.d.adj.hat_irbc,
    zeta.hat = zeta.hat, zeta.hat_irbc = zeta.hat_irbc,
    estimate_irbc = as.numeric(zeta.hat_irbc),
    Bhat.y = Bhat.components[["outcome"]],
    Bhat.d = Bhat.components[["treatment"]],
    Bhat.z = Bhat.z,
    Bhat.y.adj = Bhat.y.adj, Bhat.d.adj = Bhat.d.adj,
    Bhat.cfrd = Bhat.cfrd,
    frd.weights = frd.weights, cfrd.weights = cfrd.weights,
    component.estimates = tau.components,
    component.estimates.irbc = tau.components.irbc,
    component.bias = Bhat.components,
    cov.irbc = cov.irbc, cov.us = cov.us,
    se_irbc = se_irbc, se_us = se_us, ci_irbc = ci_irbc,
    epshat.plus = eps.plus, epshat.minus = eps.minus,
    sharp.outcome = sharp.y, sharp.treatment = sharp.d,
    sharp.covariates = sharp.covs,
    fuzzy.boot.weights = fuzzy.boot.weights,
    gamma.boot.weights = gamma.boot.weights
  )

  if (boot.ci == TRUE) {
    if (!is.null(boot.seed)) set.seed(as.integer(boot.seed))
    if (!is.finite(se_us) || se_us <= 0) {
      stop("Covariate-adjusted fuzzy RD bootstrap calibration requires a positive non-bias-corrected standard error.")
    }
    mhat <- as.numeric(se_irbc / se_us)
    if (!is.finite(mhat) || mhat <= 0) {
      stop("Covariate-adjusted fuzzy RD bootstrap calibration could not be computed.")
    }
    H.inv <- function(u) pnorm(mhat * qnorm(u))
    boot.probs <- H.inv(c(alpha / 2, 1 - alpha / 2))

    component.data <- cbind(
      outcome = as.numeric(y), treatment = as.numeric(d), covs
    )
    eps.full <- matrix(0, nrow = length(x), ncol = n.components,
                       dimnames = list(NULL, component.names))
    eps.full[indp, ] <- eps.plus
    eps.full[indm, ] <- eps.minus

    boot.local <- (indp & abs(x - c) <= sharp.y$h.plus) |
      (indm & abs(x - c) <= sharp.y$h.minus)
    if (!any(boot.local & indp) || !any(boot.local & indm)) {
      stop("Bootstrap CI requires observations on both sides of the cutoff inside the bandwidth.")
    }
    boot.local.indices <- which(boot.local)
    ghat.boot.mean <- matrix(0, nrow = length(x), ncol = n.components,
                             dimnames = list(NULL, component.names))
    for (ii in boot.local.indices) {
      if (x[ii] >= c) {
        ghat.boot.mean[ii, ] <- .irbcrd_lp_component_estimate(
          values = component.data[indp, , drop = FALSE],
          x = xp, x0 = x[ii], h = sharp.y$h.plus, p = p,
          kersel = sharp.y$kersel.plus, bwcheck = bwcheck
        )
      } else {
        ghat.boot.mean[ii, ] <- .irbcrd_lp_component_estimate(
          values = component.data[indm, , drop = FALSE],
          x = xm, x0 = x[ii], h = sharp.y$h.minus, p = p,
          kersel = sharp.y$kersel.minus, bwcheck = bwcheck
        )
      }
    }

    omega.plus <- .irbcrd_lp_component_weights(
      x = xp, x0 = c, h = sharp.y$h.plus, p = p,
      kersel = sharp.y$kersel.plus, bwcheck = bwcheck
    )
    omega.minus <- .irbcrd_lp_component_weights(
      x = xm, x0 = c, h = sharp.y$h.minus, p = p,
      kersel = sharp.y$kersel.minus, bwcheck = bwcheck
    )
    ghat.plus <- as.numeric(crossprod(
      omega.plus, component.data[indp, , drop = FALSE]
    ))
    ghat.minus <- as.numeric(crossprod(
      omega.minus, component.data[indm, , drop = FALSE]
    ))

    boot.t <- rep(NA_real_, boot.B)
    boot.t.components <- matrix(
      NA_real_, nrow = boot.B, ncol = n.components,
      dimnames = list(NULL, component.names)
    )
    boot.linear.weights <- matrix(
      NA_real_, nrow = boot.B, ncol = n.components,
      dimnames = list(NULL, component.names)
    )
    boot.fuzzy.weights <- matrix(
      NA_real_, nrow = boot.B, ncol = 2L,
      dimnames = list(NULL, names(frd.weights))
    )
    boot.gamma.y <- matrix(
      NA_real_, nrow = boot.B, ncol = ncol(covs),
      dimnames = list(NULL, covariate.names)
    )
    boot.gamma.d <- matrix(
      NA_real_, nrow = boot.B, ncol = ncol(covs),
      dimnames = list(NULL, covariate.names)
    )
    boot.adjusted.first.stage <- rep(NA_real_, boot.B)

    for (bb in seq_len(boot.B)) {
      e.boot <- .irbcrd_wild_draws(length(boot.local.indices), boot.wild)
      component.boot <- component.data
      component.boot[boot.local.indices, ] <-
        ghat.boot.mean[boot.local.indices, , drop = FALSE] +
        sweep(eps.full[boot.local.indices, , drop = FALSE],
              MARGIN = 1L, STATS = e.boot, FUN = "*")

      ghat.plus.boot <- as.numeric(crossprod(
        omega.plus, component.boot[indp, , drop = FALSE]
      ))
      ghat.minus.boot <- as.numeric(crossprod(
        omega.minus, component.boot[indm, , drop = FALSE]
      ))
      tau.components.boot <- ghat.plus.boot - ghat.minus.boot
      names(tau.components.boot) <- component.names
      T.components <- as.numeric(
        as.numeric(sharp.y$Q.plus) * (ghat.plus.boot - ghat.plus) -
          as.numeric(sharp.y$Q.minus) * (ghat.minus.boot - ghat.minus)
      )
      names(T.components) <- component.names

      # Combine the paper's fixed fuzzy weights and fixed covariate slopes.
      gamma.y.use <- gamma.y.hat
      gamma.d.use <- gamma.d.hat

      tau.z.boot <- tau.components.boot[covariate.names]
      tau.d.adj.boot <- as.numeric(
        tau.components.boot[["treatment"]] -
          crossprod(gamma.d.use, tau.z.boot)
      )
      fuzzy.weights.use <- frd.weights

      linear.weights.use <- c(
        outcome = fuzzy.weights.use[["reduced.form"]],
        treatment = fuzzy.weights.use[["first.stage"]],
        -(fuzzy.weights.use[["reduced.form"]] * gamma.y.use +
            fuzzy.weights.use[["first.stage"]] * gamma.d.use)
      )
      names(linear.weights.use) <- component.names
      if (any(!is.finite(c(T.components, linear.weights.use)))) {
        next
      }

      boot.t.components[bb, ] <- T.components
      boot.fuzzy.weights[bb, ] <- fuzzy.weights.use
      boot.gamma.y[bb, ] <- gamma.y.use
      boot.gamma.d[bb, ] <- gamma.d.use
      boot.linear.weights[bb, ] <- linear.weights.use
      boot.adjusted.first.stage[bb] <- tau.d.adj.boot
      boot.t[bb] <- as.numeric(crossprod(linear.weights.use, T.components))
    }

    boot.t.ok <- boot.t[is.finite(boot.t)]
    boot.success <- length(boot.t.ok)
    if (boot.success < 2L) {
      warning("Fewer than two successful bootstrap repetitions; bootstrap CI is NA.")
      ci_irbc_boot <- c(NA_real_, NA_real_)
    } else {
      q.boot <- as.numeric(quantile(
        boot.t.ok, probs = boot.probs, names = FALSE, type = 7
      ))
      ci_irbc_boot <- c(zeta.hat - q.boot[2L],
                        zeta.hat - q.boot[1L])
    }

    out <- c(out, list(
      estimate_irbc_boot = as.numeric(zeta.hat_irbc),
      ci_irbc_boot = ci_irbc_boot,
      boot.B = boot.B, boot.success = boot.success,
      boot.impl = "plug-in", boot.probs = boot.probs,
      boot.mhat = mhat, boot.wild = boot.wild,
      bootstrap.weights = paste0(
        fuzzy.boot.weights, "-fuzzy/", gamma.boot.weights, "-gamma"
      )
    ))
    if (boot.keep == TRUE) {
      out$boot.t <- boot.t
      out$boot.t.components <- boot.t.components
      out$boot.linear.weights <- boot.linear.weights
      out$boot.fuzzy.weights <- boot.fuzzy.weights
      out$boot.gamma.y <- boot.gamma.y
      out$boot.gamma.d <- boot.gamma.d
      out$boot.adjusted.first.stage <- boot.adjusted.first.stage
    }
  }

  out
}
