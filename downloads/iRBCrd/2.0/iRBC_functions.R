######################################################################################
## Auxiliary functions for iRBCrd
## Authors: Giuseppe Cavaliere, Sílvia Gonçalves, Morten Ø. Nielsen, & Edoardo Zanelli
## version: 2.0 (22.10.2026)
######################################################################################

.bwcheck_default <- 21

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
