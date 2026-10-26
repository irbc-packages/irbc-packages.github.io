######################################################################################
## Auxiliary functions for iRBCrd
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 2.1 (26.10.2026)
######################################################################################

.bwcheck_default <- 21

# Exact coefficient maps keep kink bias, variance, and bootstrap fits aligned.
.irbcrd_kink_lp <- function(x, x0, h, p, kersel, deriv = 0,
                            bwcheck = NULL) {
  if (!length(x)) stop("No observations are available for a kink local fit.")
  h <- .regularize_h(x, x0, h, bwcheck)
  u <- (x - x0) / h
  w <- K(u, kersel)
  active <- which(w > 0)
  if (length(active) < p + 1L) {
    stop("A kink local fit has too few positive-weight observations; increase h or bwcheck.")
  }
  R <- .poly_basis(u[active], p)
  sw <- sqrt(w[active])
  qr.fit <- qr(R * sw, tol = 1e-10)
  if (qr.fit$rank < p + 1L) {
    stop("A kink local-polynomial fit is rank deficient; inspect ties or increase h.")
  }
  Q <- qr.Q(qr.fit)
  map.pivot <- sweep(backsolve(qr.R(qr.fit), t(Q)), 2L, sw, "*")
  map <- map.pivot[order(qr.fit$pivot), , drop = FALSE]
  omega <- numeric(length(x))
  omega[active] <- map[deriv + 1L, ] * factorial(deriv) / h^deriv
  leverage <- numeric(length(x))
  leverage[active] <- rowSums(Q^2)
  list(omega = omega, h = h, active = active, map = map,
       leverage = leverage)
}

.irbcrd_kink_hc <- function(raw, leverage, n.local, df, hc) {
  if (hc == "hc0") return(raw)
  factor <- if (hc == "hc1") (n.local - df) / n.local else 1 - leverage
  if (any(!is.finite(factor)) || any(factor <= 0)) {
    stop("Kink HC residuals require positive degrees of freedom and leverage adjustments.")
  }
  if (hc != "hc3") factor <- sqrt(factor)
  raw / factor
}

# Population moments from the paper's kink auxiliary results. Only the outer
# equivalent kernel is differentiated; the inner pilot is a level smoother.
.irbcrd_kink_constants <- local({
  cache <- new.env(parent = emptyenv())
  function(kersel, p, nquad = 64L, length.ratio = FALSE) {
    kersel <- .canonical_kernel(kersel)
    key <- paste(kersel, p, nquad, length.ratio, sep = ":")
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache))
    coeff <- switch(kersel, tri = c(1, -1), uni = .5, epa = c(.75, 0, -.75))
    moments <- function(t, max.order) {
      ans <- matrix(0, length(t), max.order + 1L)
      for (j in 0:max.order) for (k in seq_along(coeff)) {
        power <- j + k
        ans[, j + 1L] <- ans[, j + 1L] +
          coeff[k] * (1 + (-1)^j * t^power) / power
      }
      ans
    }
    k <- seq_len(nquad - 1L)
    J <- matrix(0, nquad, nquad)
    J[cbind(k, k + 1L)] <- J[cbind(k + 1L, k)] <- k / sqrt(4 * k^2 - 1)
    eig <- eigen(J, symmetric = TRUE)
    quad <- function(a, b) list(
      x = (a + b) / 2 + (b - a) * eig$values / 2,
      w = (b - a) * eig$vectors[1L, ]^2)
    idx <- outer(0:p, 0:p, "+") + 1L
    m0 <- drop(moments(0, 2 * p + 2L))
    beta <- solve(matrix(m0[idx], p + 1L), diag(p + 1L)[, 2L])
    original <- function(u) {
      w <- drop(.poly_basis(u, p) %*% beta) * K(u, kersel)
      w[u < 0 | u > 1] <- 0
      w
    }
    pilot.beta <- function(t) {
      mm <- moments(t, 2 * p)
      e0 <- c(1, rep(0, p))
      t(vapply(seq_along(t), function(i)
        solve(matrix(mm[i, idx], p + 1L), e0), numeric(p + 1L)))
    }
    q01 <- quad(0, 1)
    C <- sum(beta * m0[(p + 2L):(2 * p + 2L)])
    pilot.C <- rowSums(pilot.beta(q01$x) *
      moments(q01$x, 2 * p + 1L)[, (p + 2L):(2 * p + 2L), drop = FALSE])
    C.lp <- sum(original(q01$x) * q01$w * pilot.C)
    if (!is.finite(C.lp) || abs(C.lp) < 1e-12) {
      stop("The population kink bias constant is numerically zero for this specification.")
    }
    out <- list(C = C, C.lp = C.lp, Q = C / C.lp)
    if (length.ratio) {
      q12 <- quad(1, 2)
      u <- c(q01$x, q12$x)
      uw <- c(q01$w, q12$w)
      convolution <- vapply(u, function(point) {
        cuts <- sort(unique(c(max(0, point - 1), min(1, point), 1)))
        ans <- 0
        for (j in seq_len(length(cuts) - 1L)) {
          rule <- quad(cuts[j], cuts[j + 1L])
          v <- point - rule$x
          pilot <- rowSums(pilot.beta(rule$x) * .poly_basis(v, p)) * K(v, kersel)
          ans <- ans + sum(original(rule$x) * rule$w * pilot)
        }
        ans
      }, numeric(1))
      corrected <- (1 + out$Q) * original(u) - out$Q * convolution
      q <- p + 1L
      idxq <- outer(0:q, 0:q, "+") + 1L
      betaq <- solve(matrix(m0[idxq], q + 1L), diag(q + 1L)[, 2L])
      rbc <- drop(.poly_basis(q01$x, q) %*% betaq) * K(q01$x, kersel)
      out$length.ratio <- sqrt(sum(corrected^2 * uw) / sum(rbc^2 * q01$w))
    }
    assign(key, out, cache)
    out
  }
})

.validate_bwcheck <- function(bwcheck) {
  if (is.null(bwcheck)) {
    return(NULL)
  }
  if (!is.numeric(bwcheck) || length(bwcheck) != 1 || !is.finite(bwcheck) ||
      bwcheck < 1 || bwcheck != floor(bwcheck)) {
    stop("'bwcheck' must be NULL or a positive integer")
  }
  as.integer(bwcheck)
}

.regularize_h <- function(x, x0, h, bwcheck = .bwcheck_default) {
  bwcheck <- .validate_bwcheck(bwcheck)
  if (!is.null(bwcheck)) {
    d <- abs(x - x0)
    if (length(d) > 0) {
      bw.index <- min(length(d), bwcheck)
      if (length(d) >= bwcheck) {
        bw.min <- sort(d, partial = bw.index)[bw.index]
      } else {
        bw.min <- sort(d)[bw.index]
      }
      h <- max(h, bw.min)
    }
  }
  h
}

.poly_basis <- function(u, degree) {
  u <- as.vector(u)
  if (degree == 0) {
    return(matrix(1, nrow = length(u), ncol = 1))
  }
  if (degree == 1) {
    return(cbind(1, u))
  }
  outer(u, 0:degree, `^`)
}

.canonical_kernel <- function(kersel) {
  if (!is.character(kersel) || length(kersel) != 1) {
    stop("Invalid kernel specification.")
  }

  switch(kersel,
         epanechnikov = "epa",
         epa = "epa",
         triangular = "tri",
         tri = "tri",
         uniform = "uni",
         uni = "uni",
         stop("Invalid kernel specification."))
}

# Kernel function
K <- function(u, kersel) {
  kersel <- .canonical_kernel(kersel)
  au <- abs(u)

  result <- if (kersel == "epa") {

    ifelse(au <= 1, 0.75 * (1 - u^2), 0)

  } else if (kersel == "tri") {

    ifelse(au <= 1, 1 - au, 0)

  } else if (kersel == "uni") {

    ifelse(au <= 1, 0.5, 0)

  } else {
    stop("Invalid kernel specification.")
  }

  return(as.numeric(result))
}


# Polynomial basis function
r_p <- function(p) {
  function(u) {
    .poly_basis(u, p)
  }
}


# Inverse of the Gamma_p matrix
Gamma.p.inv <- function(u, h, p, kersel) {
  n <- length(u)
  X <- .poly_basis(u, p)
  K.X <- K(u = u, kersel = kersel) / h
  nonzero_indices <- which(K.X != 0)
  X <- X[nonzero_indices, , drop = FALSE]
  K.X <- K.X[nonzero_indices]
  Gamma.inv <- qrXXinv(sqrt(K.X) * X) * n
}


# Local polynomial estimator
LP.est <- function(y, x, x0, h, p, v, kersel, bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)

  n <- length(x)
  u <- (x - x0) / h
  X <- .poly_basis(u, p)
  K.X <- K(u = u, kersel = kersel) / h
  nonzero_indices <- which(K.X != 0)
  X <- X[nonzero_indices, , drop = FALSE]
  K.X <- K.X[nonzero_indices]

  Gamma.inv <- qrXXinv(sqrt(K.X) * X) * n

  e1 <- numeric(p + 1)
  e1[v + 1] <- 1
  ghat <- factorial(v) * e1 %*% Gamma.inv %*% crossprod(X * K.X, y[nonzero_indices]) / (n * (h^v))
  return(ghat)
}


# Equivalent kernel of a local polynomial estimator evaluated at a single x_i
wix0 <- function(x, xi, x0, h, p, v, kersel, bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)

  e1 <- numeric(p + 1)
  e1[v + 1] <- 1
  u <- (x - x0) / h
  Gamma.inv <- Gamma.p.inv(u = u, h = h, p = p, kersel = kersel)
  coef <- as.numeric(factorial(v) * e1 %*% Gamma.inv)

  W <- numeric(length(xi))
  u_all <- (xi - x0) / h
  K_all <- K(u = u_all, kersel = kersel)
  nonzero_indices.2 <- which(K_all != 0)
  if (length(nonzero_indices.2) > 0) {
    R <- .poly_basis(u_all[nonzero_indices.2], p)
    W[nonzero_indices.2] <- as.numeric((R * K_all[nonzero_indices.2]) %*% coef)
  }
  return(W)
}

## -- Equivalent kernel of a (iRBC based) debiased LP estimator evaluated at x_i -- ##

.wiRBC_vec <- function(x, xi, x0, h, p, kersel, Q,
                       bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)
  n <- length(x)
  A <- c(1 + Q) * wix0(x = x, xi = xi, x0 = x0, h = h, p = p, v = 0,
                       kersel = kersel, bwcheck = bwcheck)
  B <- numeric(length(xi))

  U.Z <- abs((x - x0) / h) <= 1
  nonzero_indices <- which(U.Z != 0)
  if (length(nonzero_indices) > 0) {
    W0 <- wix0(x = x, xi = x[nonzero_indices], x0 = x0, h = h, p = p, v = 0,
               kersel = kersel, bwcheck = bwcheck)
    for (jj_pos in seq_along(nonzero_indices)) {
      jj <- nonzero_indices[jj_pos]
      B <- B + c(Q) * W0[jj_pos] *
        wix0(x = x, xi = xi, x0 = x[jj], h = h, p = p, v = 0,
             kersel = kersel, bwcheck = bwcheck) / (n * h)
    }
  }
  A - B
}

wiRBC <- function(x, xi, x0, h, p, kersel, Q, bwcheck = .bwcheck_default) {
  .wiRBC_vec(x = x, xi = xi, x0 = x0, h = h, p = p, kersel = kersel,
             Q = Q, bwcheck = bwcheck)
}


# Variance estimator of a (iRBC-based) bias corrected test statistic
viRBC <- function(x, x0, h, p, kersel, epshat, Q,
                  bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)
  n <- length(x)
  U.Z <- (abs((x - x0) / (2 * h)) <= 1)
  nonzero_indices <- which(U.Z != 0)
  if (length(nonzero_indices) == 0) {
    return(0)
  }

  W <- .wiRBC_vec(x = x, xi = x[nonzero_indices], x0 = x0, h = h,
                  p = p, kersel = kersel, Q = Q, bwcheck = bwcheck)
  C <- sum(((W * epshat[nonzero_indices])^2) / (n * h))
  sqrt(C)
}


# CCT (single bandwidth) standard errors
CCT.se <- function(x, x0, h, p, kersel, epshat,
                   bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)
  r.pp1 <- .poly_basis(x - x0, p + 1)
  r.p <- r.pp1[, 1:(p + 1), drop = FALSE]

  e.p1 <- matrix(0, (p + 2), 1)
  e.p1[p + 2] <- 1
  K.X <- K(u = (x - x0) / h, kersel = kersel) / h
  invG.q <- qrXXinv(sqrt(K.X) * r.pp1)
  invG.p <- qrXXinv(sqrt(K.X) * r.p)
  L <- crossprod(r.p * K.X, ((x - x0) / h)^(p + 1))

  Q.q <- t(t(r.p * K.X) - (h^(p + 1)) * (L %*% t(e.p1)) %*% t(t(invG.q %*% t(r.pp1)) * K.X))
  Q.p <- r.p * K.X

  MM.p <- crossprod(c(epshat) * Q.p)
  MM.q <- crossprod(c(epshat) * Q.q)
  MMM.p <- invG.p %*% MM.p %*% invG.p
  MMM.q <- invG.p %*% MM.q %*% invG.p

  M.p <- sqrt(MMM.p[1, 1])
  M.q <- sqrt(MMM.q[1, 1])

  ses <- list(se.us = M.p, se.rb = M.q)
  return(ses)
}


# Local-polynomial standard error without bias-correction adjustment
LP.se <- function(x, x0, h, p, kersel, epshat,
                  bwcheck = .bwcheck_default) {
  h <- .regularize_h(x = x, x0 = x0, h = h, bwcheck = bwcheck)

  r.p <- .poly_basis(x - x0, p)

  K.X <- K(u = (x - x0) / h, kersel = kersel) / h
  invG.p <- qrXXinv(sqrt(K.X) * r.p)
  Q.p <- r.p * K.X

  MM.p <- crossprod(c(epshat) * Q.p)
  MMM.p <- invG.p %*% MM.p %*% invG.p

  sqrt(MMM.p[1, 1])
}

qrXXinv <- function(x, ...) {
  chol2inv(chol(crossprod(x)))
}
