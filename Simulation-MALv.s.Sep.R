###############################################################################
## Simulation study: comparison of Joint-MAL and Sep-AL
##
## Companion code for
##   "Bayesian Joint Quantile Regression of Longitudinal Ordinal-Continuous
##    Biomarkers and Time-to-Event Outcomes" (Niu and Zhao).
##
## This script reproduces the comparison with the separate-fitting baseline
## (Section 4.3 of the paper). It is self-contained and does not source any
## other file.
##
## Dependencies: base R, splines, parallel; GIGrvg is required only when the
## number of longitudinal outcomes p is not 1 or 3.
##
## ---------------------------------------------------------------------------
## Method 1: Joint-MAL (proposed method)
##   Joint quantile regression of three longitudinal outcomes (one ordinal,
##   two continuous) and a time-to-event outcome. The longitudinal errors
##   follow a multivariate asymmetric Laplace (MAL) distribution whose
##   correlation matrix Psi is estimated, and the p-dimensional random-effect
##   vector b_i enters the hazard through psi' b_i.
##
## Method 2: Sep-AL (separate-fitting baseline)
##   Each longitudinal outcome j = 1, ..., p is modeled with a univariate
##   asymmetric Laplace (AL) quantile regression, which is equivalent to
##   constraining Psi = I and ignores the cross-outcome dependence. The
##   random effects b_ij are scalars with priors N(0, sigma2_{b_j}); they are
##   updated within a single chain together with the shared survival
##   submodel.
##
## Assembly of the Sep-AL estimates:
##   beta_tau : row j is taken from the j-th submodel   (p * K = 9 entries)
##   psi      : entry j is taken from the j-th submodel (p     = 3 entries)
##   gamma    : estimated jointly in the single chain    (nq    = 2 entries)
##
## ---------------------------------------------------------------------------
## Aggregate accuracy measure over the 14 parameters (9 beta, 2 gamma, 3 psi):
##   APMSE = (1/S) sum_s [ ||beta_hat - beta||_F^2
##                       + ||gamma_hat - gamma||_2^2
##                       + ||psi_hat   - psi||_2^2 ]
##
## Settings: N = 500, low censoring (u0 = 6),
##           2 correlation levels (Psi off-diagonal 0.3 / 0.6)
##           x 3 quantile combinations x 100 replicates.
###############################################################################

library(splines)
library(parallel)


## =========================================================================
## 0. Utility functions
## =========================================================================

## Draw from N(mu, sd^2) truncated to (lo, hi), vectorized (inverse-CDF method).
rtnorm <- function(mu, sd, lo, hi) {
  n <- length(mu)
  a <- (lo - mu) / sd; b <- (hi - mu) / sd
  flip <- a > 0
  a2 <- ifelse(flip, -b, a); b2 <- ifelse(flip, -a, b)
  pa <- pnorm(a2); pb <- pnorm(b2)
  z  <- qnorm(pa + runif(n) * (pb - pa))
  z  <- ifelse(flip, -z, z)
  bad <- !is.finite(z) | z < a | z > b
  if (any(bad)) {
    z[bad] <- ifelse(is.finite(a[bad]) & (a[bad] > 0), a[bad],
                     ifelse(is.finite(b[bad]) & (b[bad] < 0), b[bad], 0))
  }
  mu + sd * z
}

## Draw from the inverse Gaussian distribution IG(mu, shape), vectorized
## (Michael, Schucany and Haas, 1976).
rinvgauss_vec <- function(mu, shape) {
  n <- length(mu)
  y <- rnorm(n)^2
  x <- mu + mu^2 * y / (2 * shape) -
    mu / (2 * shape) * sqrt(4 * mu * shape * y + mu^2 * y^2)
  u <- runif(n)
  pmax(ifelse(u <= mu / (mu + x), x, mu^2 / x), 1e-12)
}

## Draw from the generalized inverse Gaussian distribution GIG(lambda, chi, psi),
## with density proportional to x^(lambda-1) exp{-(chi/x + psi*x)/2}.
## lambda = +-1/2 (p = 1 or 3) uses exact inverse-Gaussian draws; other values use GIGrvg.
rgig_vec <- function(lambda, chi, psi) {
  chi <- pmax(chi, 1e-10)
  if (abs(lambda + 0.5) < 1e-12) {
    return(rinvgauss_vec(sqrt(chi / psi), chi))
  } else if (abs(lambda - 0.5) < 1e-12) {
    return(1 / rinvgauss_vec(sqrt(psi / chi), psi))
  } else {
    if (!requireNamespace("GIGrvg", quietly = TRUE))
      stop("p 不等于 1 或 3 时需要安装 GIGrvg 包")
    return(vapply(seq_along(chi), function(k)
      GIGrvg::rgig(1, lambda, chi[k], psi[k]), numeric(1)))
  }
}

## Draw from the inverse-Wishart distribution IW(df, S).
riwish <- function(df, S) {
  S <- (S + t(S)) / 2
  W <- rWishart(1, df, chol2inv(chol(S)))[, , 1]
  X <- chol2inv(chol(W))
  (X + t(X)) / 2
}

## Log density (up to a constant) of N(mu, S) evaluated at x.
ldmvn <- function(x, mu, S) {
  R <- chol(S)
  z <- backsolve(R, x - mu, transpose = TRUE)
  -sum(log(diag(R))) - 0.5 * sum(z^2)
}

## Nodes x and weights w of n-point Gauss-Legendre quadrature on [-1, 1].
gauss_legendre <- function(n) {
  i <- seq_len(n - 1)
  bb <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- bb; J[cbind(i + 1, i)] <- bb
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = e$values[o], w = 2 * e$vectors[1, o]^2)
}

## Map unconstrained values u to the Cholesky factor L of a correlation matrix
## through canonical partial correlations z = tanh(u) (Lewandowski et al., 2009).
u2cholcorr <- function(u, p) {
  Z <- matrix(0, p, p); Z[lower.tri(Z)] <- tanh(u)
  L <- matrix(0, p, p); L[1, 1] <- 1
  for (i in 2:p) {
    ss <- 0
    for (j in 1:(i - 1)) { L[i, j] <- Z[i, j] * sqrt(1 - ss); ss <- ss + L[i, j]^2 }
    L[i, i] <- sqrt(1 - ss)
  }
  list(L = L, z = Z[lower.tri(Z)], col = col(Z)[lower.tri(Z)])
}

## Log full conditional of Psi on the unconstrained scale u, including the LKJ(eta)
## prior, the likelihood term with scatter matrix S, and the Jacobian of tanh.
logpost_Psi <- function(u, S, n, eta, p) {
  cc <- u2cholcorr(u, p)
  Li <- forwardsolve(cc$L, diag(p))
  -n * sum(log(diag(cc$L))) - 0.5 * sum((Li %*% S) * Li) +
    sum((eta + (p - 1 - cc$col) / 2) * log(1 - cc$z^2))
}

## MAL parameters: skewness theta_j and squared scale delta_j^2 for quantile level tau_j.
theta_fun  <- function(tau) (1 - 2 * tau) / (tau * (1 - tau))
delta2_fun <- function(tau) 2 / (tau * (1 - tau))


## =========================================================================
## 1. Data generation
## =========================================================================
## Simulate one data set: longitudinal ordinal-continuous outcomes with MAL errors,
## survival times with baseline hazard exp{0.8 sin(0.5 t)}, uniform censoring on
## (0, u0), and irregular visit times (exponential gaps with rate visit_rate).
gen_data <- function(N = 300,
                     tau = c(0.5, 0.5, 0.5),
                     beta  = matrix(1, 3, 3),
                     Sigma_b = diag(3),
                     Psi = matrix(0.3, 3, 3) + diag(0.7, 3),
                     d1 = 0.2,
                     dvec_cont = NULL,
                     alpha = c(0, 1, 2, 3),
                     gamma = c(-0.5, -0.5),
                     psi   = c(-0.5, 0.5, -0.3),
                     visit_rate = 15,
                     u0 = 3,
                     tgrid_max = 30) {
  p <- nrow(beta); K <- ncol(beta)
  theta  <- theta_fun(tau)
  Lambda <- diag(sqrt(delta2_fun(tau)))
  if (is.null(dvec_cont)) dvec_cont <- rep(d1, p - 1)
  dvec <- c(d1, dvec_cont)
  D     <- diag(dvec)
  Sigma <- Lambda %*% Psi %*% Lambda
  Veps  <- D %*% Sigma %*% D
  Leps  <- t(chol(Veps))
  
  b  <- matrix(rnorm(N * p), N, p) %*% chol(Sigma_b)
  z  <- matrix(rnorm(N * length(gamma)), N)
  lp <- drop(z %*% gamma + b %*% psi)
  
  h0_true <- function(t) exp(0.8 * sin(0.5 * t))
  tg <- seq(0, tgrid_max, by = 0.001)
  hg <- h0_true(tg)
  Hg <- c(0, cumsum((hg[-1] + hg[-length(hg)]) / 2 * diff(tg)))
  E  <- -log(1 - runif(N)) / exp(lp)
  Tstar <- approx(Hg, tg, xout = E, rule = 2)$y
  Cens  <- runif(N, 0, u0)
  Tobs  <- pmin(Tstar, Cens, u0)
  Delta <- as.integer(Tstar <= pmin(Cens, u0))
  
  id <- integer(0); tt <- numeric(0)
  for (i in seq_len(N)) {
    ti <- 0
    repeat {
      nxt <- ti[length(ti)] + rexp(1, visit_rate)
      if (nxt > Tobs[i]) break
      ti <- c(ti, nxt)
    }
    id <- c(id, rep(i, length(ti))); tt <- c(tt, ti)
  }
  nobs <- length(id)
  X <- matrix(rnorm(nobs * K), nobs, K)
  
  w   <- rexp(nobs)
  eps <- sweep(matrix(0, nobs, p), 2, drop(D %*% theta), "+") * w +
    sqrt(w) * (matrix(rnorm(nobs * p), nobs, p) %*% t(Leps))
  
  Y <- X %*% t(beta) + b[id, , drop = FALSE] + eps
  y1star <- findInterval(Y[, 1], c(-Inf, alpha, Inf), left.open = TRUE)
  
  list(N = N, p = p, K = K, id = id, time = tt, X = X,
       Y2 = Y[, -1, drop = FALSE], y1star = y1star, R = length(alpha) + 1,
       Tobs = Tobs, Delta = Delta, Z = z, tau = tau, d1 = d1, dvec = dvec,
       truth = list(beta = beta, alpha = alpha, gamma = gamma, psi = psi,
                    Sigma_b = Sigma_b, Psi = Psi, b = b, Y1 = Y[, 1],
                    dvec = dvec))
}


## =========================================================================
## 2. B-spline baseline hazard
## =========================================================================
## B-spline basis of order ord with n_inner interior knots at quantiles of the observed
## times, two auxiliary knots, and K_GL-point Gauss-Legendre nodes for H0(T_i).
prep_spline <- function(Tobs, n_inner = 3, ord = 4, K_GL = 15) {
  tmax  <- max(Tobs) * 1.0001
  inner <- quantile(Tobs, probs = seq_len(n_inner) / (n_inner + 1), names = FALSE)
  aux   <- c(inner[1] / 2, (inner[n_inner] + tmax) / 2)
  knots <- c(rep(0, ord), sort(c(aux, inner)), rep(tmax, ord))
  bfun  <- function(x) splineDesign(knots, x, ord = ord, outer.ok = TRUE)
  gl <- gauss_legendre(K_GL)
  N  <- length(Tobs)
  S  <- outer(Tobs / 2, 1 + gl$x)
  Wt <- outer(Tobs / 2, gl$w)
  list(CT  = bfun(Tobs),
       Cgl = bfun(as.vector(S)),
       wgl = as.vector(Wt),
       N = N, K_GL = K_GL, L = length(knots) - ord, knots = knots)
}

## Cumulative baseline hazard H0(T_i) for all subjects via Gauss-Legendre quadrature.
H0_fun <- function(zeta, sp) {
  ex <- exp(drop(sp$Cgl %*% zeta))
  rowSums(matrix(sp$wgl * ex, sp$N, sp$K_GL))
}

## Survival log-likelihood given gamma, psi, zeta and random effects B.
surv_loglik <- function(gamma, psi, zeta, B, dat, sp, H0 = NULL) {
  if (is.null(H0)) H0 <- H0_fun(zeta, sp)
  eta <- drop(dat$Z %*% gamma + B %*% psi)
  sum(dat$Delta * (drop(sp$CT %*% zeta) + eta)) - sum(exp(eta) * H0)
}

## Metropolis-Hastings step with Gaussian proposal covariance scale^2 * H^{-1},
## where H is the negative Hessian at the current point; the Hastings correction
## accounts for the position-dependent proposal.
mh_hess <- function(cur, logpost, hess, scale) {
  S0   <- scale^2 * chol2inv(chol(hess(cur)))
  prop <- cur + drop(t(chol(S0)) %*% rnorm(length(cur)))
  H1   <- tryCatch(hess(prop), error = function(e) NULL)
  if (is.null(H1) || any(!is.finite(H1))) return(list(x = cur, acc = 0))
  S1 <- scale^2 * chol2inv(chol(H1))
  lr <- logpost(prop) - logpost(cur) + ldmvn(cur, prop, S1) - ldmvn(prop, cur, S0)
  if (is.finite(lr) && log(runif(1)) < lr) list(x = prop, acc = 1)
  else list(x = cur, acc = 0)
}


## =========================================================================
## 3. Joint-MAL: MH-within-Gibbs sampler (proposed method)
## =========================================================================
## Fit Joint-MAL by MH-within-Gibbs. Returns post-burn-in, thinned draws of
## beta, alpha, gamma, psi, zeta, Sigma_b, Psi and d.
fit_jointQR <- function(dat,
                        niter = 20000, nburn = 5000, thin = 10,
                        hyper = list(a0 = 0.5, b0 = 0.5,
                                     eta = 1,
                                     s2_d = 1.0,
                                     mb = NULL, Phib = NULL,
                                     s2_gamma = 100, s2_psi = 100, s2_zeta = 100,
                                     alpha_min = -10, alpha_max = 10),
                        n_inner = 3, K_GL = 15,
                        alpha_update = c("cowles", "gibbs"),
                        estimate_d = TRUE,
                        init = NULL, verbose = 500) {
  
  alpha_update <- match.arg(alpha_update)
  p <- dat$p; K <- dat$K; N <- dat$N; R <- dat$R
  id <- dat$id; X <- dat$X; nobs <- length(id)
  tau <- dat$tau
  theta  <- theta_fun(tau)
  Lambda <- diag(sqrt(delta2_fun(tau)))
  d1fix  <- dat$d1
  dvec   <- dat$dvec
  D      <- diag(dvec)
  Dtheta <- drop(D %*% theta)
  DL_inv <- 1 / (dvec * sqrt(delta2_fun(tau)))
  
  if (is.null(hyper$eta))  hyper$eta  <- 1
  if (is.null(hyper$mb))   hyper$mb   <- p + 1
  if (is.null(hyper$Phib)) hyper$Phib <- diag(p)
  nq <- ncol(dat$Z)
  
  sp <- prep_spline(dat$Tobs, n_inner = n_inner, K_GL = K_GL)
  L  <- sp$L
  Sig_g_inv <- diag(1 / hyper$s2_gamma, nq)
  Sig_p_inv <- diag(1 / hyper$s2_psi,   p)
  Sig_z_inv <- diag(1 / hyper$s2_zeta,  L)
  
  cp <- cumsum(table(factor(dat$y1star, levels = 1:R)))[-R] / nobs
  cp <- pmin(pmax(cp, 0.01), 0.99)
  s_y <- mean(apply(dat$Y2, 2, sd))
  alpha <- sort(qnorm(cp) * s_y)
  ab <- c(alpha[1] - s_y, alpha, alpha[R - 1] + s_y)
  y1 <- (ab[dat$y1star] + ab[dat$y1star + 1]) / 2
  Y  <- cbind(y1, dat$Y2)
  
  beta  <- t(qr.solve(X, Y))
  omega <- sqrt(abs(beta)) + 1
  lam   <- matrix(1, p, K)
  Psi   <- diag(p)
  uPsi  <- rep(0, p * (p - 1) / 2)
  sPsi  <- 1.2 / sqrt(nobs)
  acc_Psi <- 0
  log_d   <- log(dvec[-1])
  Sigma_b <- diag(p)
  B     <- matrix(0, N, p)
  w     <- rep(1, nobs)
  gamma <- rep(0, nq); psi <- rep(0, p)
  zeta  <- rep(log(sum(dat$Delta) / sum(dat$Tobs)), L)
  if (!is.null(init)) for (nm in names(init)) assign(nm, init[[nm]])
  
  sc  <- list(b = 1, gamma = 1, psi = 1, zeta = 0.5, alpha = 0.05, d = 0.1,
              beta1 = 0.05)
  acc <- list(b = 0, gamma = 0, psi = 0, zeta = 0, alpha = 0, d = 0,
              beta1 = 0); nacc <- 0
  amin <- hyper$alpha_min; amax <- hyper$alpha_max
  
  keep <- seq(nburn + thin, niter, by = thin)
  G <- length(keep); g <- 0
  out <- list(beta  = array(NA, c(G, p, K)),
              alpha = matrix(NA, G, R - 1),
              gamma = matrix(NA, G, nq),
              psi   = matrix(NA, G, p),
              zeta  = matrix(NA, G, L),
              Sigma_b = array(NA, c(G, p, p)),
              Psi   = array(NA, c(G, p, p)),
              dvec  = matrix(NA, G, p))
  
  for (it in seq_len(niter)) {
    
    if (estimate_d && p > 1) {
      dvec   <- c(d1fix, exp(log_d))
      D      <- diag(dvec); Dtheta <- drop(D %*% theta)
      DL_inv <- 1 / (dvec * sqrt(delta2_fun(tau)))
    }
    Sigma <- Lambda %*% Psi %*% Lambda
    V     <- D %*% Sigma %*% D
    Vinv  <- chol2inv(chol(V))
    Sinv  <- chol2inv(chol(Sigma))
    XB    <- X %*% t(beta)
    Bid   <- B[id, , drop = FALSE]
    
    V12 <- V[1, -1, drop = FALSE]; V22i <- solve(V[-1, -1])
    a_reg <- drop(V12 %*% V22i)
    s11   <- V[1, 1] - drop(V12 %*% V22i %*% t(V12))
    Mu  <- XB + Bid + outer(w, Dtheta)
    m1  <- Mu[, 1] + drop((Y[, -1, drop = FALSE] - Mu[, -1, drop = FALSE]) %*% a_reg)
    sd1 <- sqrt(w * s11)
    
    if (alpha_update == "cowles") {
      ll_ord <- function(al, m1_use) {
        cc <- c(-Inf, al, Inf)
        pr <- pnorm((cc[dat$y1star + 1] - m1_use) / sd1) -
          pnorm((cc[dat$y1star]     - m1_use) / sd1)
        sum(log(pmax(pr, 1e-300)))
      }
      g_cur  <- diff(c(amin, alpha))
      g_prop <- g_cur + rnorm(R - 1) * sc$alpha
      if (all(g_prop > 0)) {
        ap <- amin + cumsum(g_prop)
        if (ap[R - 1] < amax) {
          lr <- ll_ord(ap, m1) - ll_ord(alpha, m1)
          if (is.finite(lr) && log(runif(1)) < lr) {
            alpha <- ap; acc$alpha <- acc$alpha + 1
          }
        }
      }
      beta1_cur  <- beta[1, ]
      beta1_prop <- beta1_cur + rnorm(K) * sc$beta1
      if (all(abs(beta1_prop) < omega[1, ]^2)) {
        m1_prop <- m1 + drop(X %*% (beta1_prop - beta1_cur))
        lr <- ll_ord(alpha, m1_prop) - ll_ord(alpha, m1)
        if (is.finite(lr) && log(runif(1)) < lr) {
          beta[1, ] <- beta1_prop; m1 <- m1_prop; XB <- X %*% t(beta)
          acc$beta1 <- acc$beta1 + 1
        }
      }
      cuts <- c(-Inf, alpha, Inf)
      Y[, 1] <- rtnorm(m1, sd1, cuts[dat$y1star], cuts[dat$y1star + 1])
    }
    
    E    <- Y - XB - Bid
    chi  <- rowSums((E %*% Vinv) * E)
    psiG <- drop(t(theta) %*% Sinv %*% theta) + 2
    w    <- rgig_vec(1 - p / 2, chi, psiG)
    
    Cres  <- Y - XB - outer(w, Dtheta)
    s_i   <- as.vector(rowsum(1 / w, id))
    Gv    <- rowsum(Cres / w, id) %*% Vinv
    Sbinv <- chol2inv(chol(Sigma_b))
    H0    <- H0_fun(zeta, sp)
    etaz  <- drop(dat$Z %*% gamma)
    logtb <- function(Bm) {
      eb <- drop(Bm %*% psi)
      -0.5 * rowSums((Bm %*% Sbinv) * Bm) -
        0.5 * s_i * rowSums((Bm %*% Vinv) * Bm) +
        rowSums(Bm * Gv) +
        dat$Delta * eb - exp(etaz + eb) * H0
    }
    Rb  <- chol(Sbinv); Rbi <- backsolve(Rb, diag(p))
    eg  <- eigen(t(Rbi) %*% Vinv %*% Rbi, symmetric = TRUE)
    A   <- Rbi %*% eg$vectors
    dsd <- 1 / sqrt(1 + outer(s_i, eg$values))
    Bprop <- B + sc$b * ((matrix(rnorm(N * p), N, p) * dsd) %*% t(A))
    ok  <- log(runif(N)) < (logtb(Bprop) - logtb(B))
    B[ok, ] <- Bprop[ok, ]
    acc$b <- acc$b + mean(ok)
    Bid <- B[id, , drop = FALSE]
    
    Eta  <- Y - Bid - outer(w, Dtheta)
    Binv <- crossprod(X / w, X)
    Bmat <- chol2inv(chol(Binv))
    Mmat <- crossprod(Eta / w, X) %*% Bmat
    Q    <- kronecker(Binv, Vinv)
    mu_v <- as.vector(Mmat); bv <- as.vector(beta); om2 <- as.vector(omega^2)
    idx_skip <- if (alpha_update == "cowles") seq(1, p * K, by = p) else integer(0)
    for (l in setdiff(seq_along(bv), idx_skip)) {
      cm  <- mu_v[l] - sum(Q[l, -l] * (bv[-l] - mu_v[-l])) / Q[l, l]
      bv[l] <- rtnorm(cm, 1 / sqrt(Q[l, l]), -om2[l], om2[l])
    }
    beta <- matrix(bv, p, K); XB <- X %*% t(beta)
    
    Lres <- (Y - XB - Bid - outer(w, Dtheta)) / sqrt(w)
    U    <- sweep(Lres, 2, DL_inv, "*")
    S_r  <- crossprod(U)
    lp_cur <- logpost_Psi(uPsi, S_r, nobs, hyper$eta, p)
    for (rep_ in 1:5) {
      ucand  <- uPsi + sPsi * rnorm(length(uPsi))
      lp_new <- logpost_Psi(ucand, S_r, nobs, hyper$eta, p)
      if (log(runif(1)) < lp_new - lp_cur) {
        uPsi <- ucand; lp_cur <- lp_new; acc_Psi <- acc_Psi + 1
      }
    }
    Lc <- u2cholcorr(uPsi, p)$L; Psi <- Lc %*% t(Lc)
    Sigma <- Lambda %*% Psi %*% Lambda
    
    if (estimate_d && p > 1) {
      lp_d <- function(ld) {
        dv  <- c(d1fix, exp(ld)); Dv <- diag(dv)
        Vd  <- Dv %*% Sigma %*% Dv; Vdi <- chol2inv(chol(Vd))
        Er  <- Y - XB - Bid - outer(w, drop(Dv %*% theta))
        ll  <- -nobs * sum(log(dv[-1])) -
          0.5 * sum(rowSums((Er / w) * (Er %*% Vdi)))
        ll - 0.5 * sum(ld^2) / hyper$s2_d
      }
      hs_d <- function(ld) {
        eps <- 1e-4; H <- numeric(p - 1); lp0 <- lp_d(ld)
        for (j in seq_len(p - 1)) {
          ldp <- ld; ldp[j] <- ld[j] + eps
          ldm <- ld; ldm[j] <- ld[j] - eps
          H[j] <- -(lp_d(ldp) - 2 * lp0 + lp_d(ldm)) / eps^2
        }
        diag(pmax(H, 1e-6))
      }
      r_ <- mh_hess(log_d, lp_d, hs_d, sc$d)
      if (r_$acc == 1) {
        log_d  <- r_$x
        dvec   <- c(d1fix, exp(log_d)); D <- diag(dvec)
        Dtheta <- drop(D %*% theta); DL_inv <- 1 / (dvec * sqrt(delta2_fun(tau)))
        V <- D %*% Sigma %*% D; Vinv <- chol2inv(chol(V))
        acc$d <- acc$d + 1
      }
    }
    
    omega <- matrix(rexp(p * K, as.vector(lam)), p, K) + sqrt(abs(beta))
    lam   <- matrix(rgamma(p * K, 3 + hyper$a0, hyper$b0 + as.vector(omega)), p, K)
    
    if (alpha_update == "gibbs") for (r in seq_len(R - 1)) {
      lo <- max(c(Y[dat$y1star == r, 1], if (r > 1) alpha[r - 1], hyper$alpha_min))
      hi <- min(c(Y[dat$y1star == r + 1, 1], if (r < R - 1) alpha[r + 1], hyper$alpha_max))
      if (hi > lo) alpha[r] <- runif(1, lo, hi)
    }
    
    Sigma_b <- riwish(N + hyper$mb, crossprod(B) + hyper$Phib)
    
    eb   <- drop(B %*% psi)
    lp_g <- function(g) -0.5 * sum(g^2 * diag(Sig_g_inv)) +
      sum(dat$Delta * drop(dat$Z %*% g)) - sum(exp(drop(dat$Z %*% g) + eb) * H0)
    hs_g <- function(g) Sig_g_inv +
      crossprod(dat$Z * (exp(drop(dat$Z %*% g) + eb) * H0), dat$Z)
    r_ <- mh_hess(gamma, lp_g, hs_g, sc$gamma)
    gamma <- r_$x; acc$gamma <- acc$gamma + r_$acc
    etaz  <- drop(dat$Z %*% gamma)
    
    lp_p <- function(ps) -0.5 * sum(ps^2 * diag(Sig_p_inv)) +
      sum(dat$Delta * drop(B %*% ps)) - sum(exp(etaz + drop(B %*% ps)) * H0)
    hs_p <- function(ps) Sig_p_inv +
      crossprod(B * (exp(etaz + drop(B %*% ps)) * H0), B)
    r_ <- mh_hess(psi, lp_p, hs_p, sc$psi)
    psi <- r_$x; acc$psi <- acc$psi + r_$acc
    
    ee   <- exp(etaz + drop(B %*% psi))
    lp_z <- function(zt) -0.5 * sum(zt^2 * diag(Sig_z_inv)) +
      sum(dat$Delta * drop(sp$CT %*% zt)) - sum(ee * H0_fun(zt, sp))
    hs_z <- function(zt) {
      wv <- sp$wgl * exp(drop(sp$Cgl %*% zt)) * rep(ee, sp$K_GL)
      Sig_z_inv + crossprod(sp$Cgl * wv, sp$Cgl)
    }
    r_ <- mh_hess(zeta, lp_z, hs_z, sc$zeta)
    zeta <- r_$x; acc$zeta <- acc$zeta + r_$acc
    
    nacc <- nacc + 1
    if (it %% 50 == 0) {
      ar <- unlist(acc) / nacc
      acc[] <- 0; nacc <- 0
      if (it <= nburn) {
        for (nm in names(sc)) {
          if (nm == "alpha" && alpha_update == "gibbs") next
          if (nm == "beta1" && alpha_update == "gibbs") next
          if (nm == "d"     && !estimate_d)             next
          if (ar[nm] < 0.20)      sc[[nm]] <- sc[[nm]] * 0.70
          else if (ar[nm] < 0.25) sc[[nm]] <- sc[[nm]] * 0.85
          else if (ar[nm] > 0.50) sc[[nm]] <- sc[[nm]] * 1.30
          else if (ar[nm] > 0.45) sc[[nm]] <- sc[[nm]] * 1.15
          sc[[nm]] <- pmin(pmax(sc[[nm]], 1e-4), 10)
        }
        acc_Psi_rate <- acc_Psi / (5 * 50)
        if      (acc_Psi_rate < 0.20) sPsi <- sPsi * 0.70
        else if (acc_Psi_rate < 0.25) sPsi <- sPsi * 0.85
        else if (acc_Psi_rate > 0.50) sPsi <- sPsi * 1.30
        else if (acc_Psi_rate > 0.45) sPsi <- sPsi * 1.15
        sPsi <- pmin(pmax(sPsi, 1e-4), 5)
        acc_Psi <- 0
      }
    }
    
    if (it %in% keep) {
      g <- g + 1
      out$beta[g, , ]    <- beta
      out$alpha[g, ]     <- alpha
      out$gamma[g, ]     <- gamma
      out$psi[g, ]       <- psi
      out$zeta[g, ]      <- zeta
      out$Sigma_b[g, , ] <- Sigma_b
      out$Psi[g, , ]     <- Psi
      out$dvec[g, ]      <- dvec
    }
    if (verbose > 0 && it %% verbose == 0)
      cat(sprintf("  [Joint] iter %5d\n", it))
  }
  out$acc_Psi <- acc_Psi / (5 * (niter - nburn))
  out
}


## =========================================================================
## 4. Sep-AL: univariate AL submodels with a shared survival submodel
##
## Differences from Joint-MAL:
##   (1) Errors: each outcome j has its own univariate AL distribution
##       (no Psi, i.e. Psi = I).
##   (2) Random effects: b_ij are scalars with priors b_ij ~ N(0, sigma2_{b_j}).
##
## Components identical to Joint-MAL:
##   - b_ij updates include the survival contribution (through psi_j);
##   - gamma, psi and zeta are updated by the same Hessian-based MH steps;
##   - cut-points alpha and beta_{1.} are updated by Cowles collapsed MH;
##   - B-spline baseline hazard and Gauss-Legendre quadrature;
##   - prior specification and adaptive step-size rules.
## =========================================================================

## Fit Sep-AL by MH-within-Gibbs (same output structure as fit_jointQR, with the
## scalar random-effect variances s2b in place of Sigma_b and Psi).
fit_sepAL <- function(dat,
                      niter = 20000, nburn = 5000, thin = 10,
                      hyper = list(a0 = 0.5, b0 = 0.5,
                                   s2_d     = 1.0,
                                   mb       = 2,
                                   Phib     = 1,
                                   s2_gamma = 100,
                                   s2_psi   = 100,
                                   s2_zeta  = 100,
                                   alpha_min = -10,
                                   alpha_max =  10),
                      n_inner = 3, K_GL = 15,
                      estimate_d = TRUE,
                      verbose = 500) {
  
  ## ------- Dimensions -------
  p    <- dat$p; K <- dat$K; N <- dat$N; R <- dat$R
  id   <- dat$id; X <- dat$X; nobs <- length(id)
  tau  <- dat$tau
  nq   <- ncol(dat$Z)
  d1fix <- dat$d1
  
  ## AL parameters of each univariate submodel
  theta_v  <- theta_fun(tau)          # length p
  delta2_v <- delta2_fun(tau)         # length p
  
  ## ------- Default hyperparameters -------
  if (is.null(hyper$a0))        hyper$a0        <- 0.5
  if (is.null(hyper$b0))        hyper$b0        <- 0.5
  if (is.null(hyper$s2_d))      hyper$s2_d      <- 1.0
  if (is.null(hyper$mb))        hyper$mb        <- 2
  if (is.null(hyper$Phib))      hyper$Phib      <- 1
  if (is.null(hyper$s2_gamma))  hyper$s2_gamma  <- 100
  if (is.null(hyper$s2_psi))    hyper$s2_psi    <- 100
  if (is.null(hyper$s2_zeta))   hyper$s2_zeta   <- 100
  if (is.null(hyper$alpha_min)) hyper$alpha_min <- -10
  if (is.null(hyper$alpha_max)) hyper$alpha_max <-  10
  amin <- hyper$alpha_min; amax <- hyper$alpha_max
  
  ## ------- B-spline baseline hazard -------
  sp <- prep_spline(dat$Tobs, n_inner = n_inner, K_GL = K_GL)
  L  <- sp$L
  Sig_g_inv <- diag(1 / hyper$s2_gamma, nq)
  Sig_p_inv <- diag(1 / hyper$s2_psi,   p)
  Sig_z_inv <- diag(1 / hyper$s2_zeta,  L)
  
  ## ------- Initial values -------
  ## Cut-points and latent ordinal scores (j = 1)
  s_y   <- mean(apply(dat$Y2, 2, sd))
  cp    <- cumsum(table(factor(dat$y1star, levels = 1:R)))[-R] / nobs
  cp    <- pmin(pmax(cp, 0.01), 0.99)
  alpha <- sort(qnorm(cp) * s_y)
  ab    <- c(alpha[1] - s_y, alpha, alpha[R-1] + s_y)
  y1    <- (ab[dat$y1star] + ab[dat$y1star + 1]) / 2
  
  ## Response matrix Y (rows = observations, columns = outcomes);
  ## the first column holds the latent ordinal scores
  Y <- cbind(y1, dat$Y2)   # nobs x p
  
  ## Regression coefficients beta (p x K); row j belongs to submodel j
  beta  <- t(qr.solve(X, Y))
  omega <- sqrt(abs(beta)) + 1
  lam   <- matrix(1, p, K)
  
  ## Random effects: B is N x p; column j holds the scalar b_ij
  B     <- matrix(0, N, p)
  s2b   <- rep(1, p)          # scalar random-effect variance of each submodel
  
  ## Latent scale variables w (one per submodel), stored as an nobs x p matrix
  W <- matrix(1, nobs, p)
  
  ## Scale vector d (d_1 fixed at d1fix; d_j, j >= 2, estimated)
  dvec  <- dat$dvec
  log_d <- log(dvec[-1])      # length p - 1
  
  ## Survival parameters
  gamma <- rep(0, nq)
  psi   <- rep(0, p)
  zeta  <- rep(log(sum(dat$Delta) / sum(dat$Tobs)), L)
  
  ## ------- Step sizes and acceptance counters -------
  sc_alpha <- 0.05; sc_beta1 <- 0.05
  sc_b     <- rep(1,    p)
  sc_d     <- rep(0.1,  p - 1)
  sc_gamma <- 1; sc_psi <- 1; sc_zeta <- 0.5
  
  acc_alpha <- 0; acc_beta1 <- 0
  acc_b     <- rep(0, p)
  acc_d     <- rep(0, p - 1)
  acc_gamma <- 0; acc_psi <- 0; acc_zeta <- 0
  nacc      <- 0
  
  ## ------- Storage -------
  keep <- seq(nburn + thin, niter, by = thin)
  G    <- length(keep); g <- 0
  out  <- list(
    beta    = array(NA, c(G, p, K)),
    alpha   = matrix(NA, G, R - 1),
    gamma   = matrix(NA, G, nq),
    psi     = matrix(NA, G, p),
    zeta    = matrix(NA, G, L),
    s2b     = matrix(NA, G, p),
    dvec    = matrix(NA, G, p))
  
  ## =======================================================================
  ## Main MCMC loop
  ## =======================================================================
  for (it in seq_len(niter)) {
    
    ## Quantities shared by the updates in this iteration
    Bid  <- B[id, , drop = FALSE]   # nobs x p
    XB   <- X %*% t(beta)           # nobs x p
    etaz <- drop(dat$Z %*% gamma)
    H0   <- H0_fun(zeta, sp)
    
    ## -------------------------------------------------------------------
    ## (A) Cowles collapsed MH for the cut-points alpha and beta_{1.} (j = 1)
    ## -------------------------------------------------------------------
    d1   <- dvec[1]; th1 <- theta_v[1]; dl1 <- delta2_v[1]
    Vsc1 <- d1^2 * dl1
    ## Sep-AL: outcome 1 is modeled alone, without conditioning on outcomes 2..p
    m1   <- XB[, 1] + Bid[, 1] + W[, 1] * d1 * th1
    sd1  <- sqrt(W[, 1] * Vsc1)
    
    ll_ord <- function(al, m1_use) {
      cc <- c(-Inf, al, Inf)
      pr <- pnorm((cc[dat$y1star + 1] - m1_use) / sd1) -
        pnorm((cc[dat$y1star]     - m1_use) / sd1)
      sum(log(pmax(pr, 1e-300)))
    }
    
    ## Update alpha
    g_cur  <- diff(c(amin, alpha))
    g_prop <- g_cur + rnorm(R - 1) * sc_alpha
    if (all(g_prop > 0)) {
      ap <- amin + cumsum(g_prop)
      if (ap[R - 1] < amax) {
        lr <- ll_ord(ap, m1) - ll_ord(alpha, m1)
        if (is.finite(lr) && log(runif(1)) < lr) {
          alpha <- ap; acc_alpha <- acc_alpha + 1
        }
      }
    }
    
    ## Update beta_{1.}
    b1_cur  <- beta[1, ]
    b1_prop <- b1_cur + rnorm(K) * sc_beta1
    if (all(abs(b1_prop) < omega[1, ]^2)) {
      m1_prop <- m1 + drop(X %*% (b1_prop - b1_cur))
      lr <- ll_ord(alpha, m1_prop) - ll_ord(alpha, m1)
      if (is.finite(lr) && log(runif(1)) < lr) {
        beta[1, ] <- b1_prop; m1 <- m1_prop
        XB <- X %*% t(beta); acc_beta1 <- acc_beta1 + 1
      }
    }
    
    ## Draw the latent scores y_{1it} from truncated normals
    cuts   <- c(-Inf, alpha, Inf)
    Y[, 1] <- rtnorm(m1, sd1, cuts[dat$y1star], cuts[dat$y1star + 1])
    
    ## -------------------------------------------------------------------
    ## (B-E) Per-submodel updates of w_j, beta_j (j >= 2), d_j and omega_j
    ## -------------------------------------------------------------------
    for (j in seq_len(p)) {
      dj  <- dvec[j]; thj <- theta_v[j]; dl2j <- delta2_v[j]
      Vscj <- dj^2 * dl2j; dthj <- dj * thj
      ej   <- Y[, j] - XB[, j] - Bid[, j]
      
      ## w_j ~ GIG
      chi_j <- ej^2 / Vscj
      psiG_j <- thj^2 / dl2j + 2
      W[, j] <- rgig_vec(0.5, chi_j, psiG_j)
      
      ## beta_j (skipped for j = 1, which is updated in step A)
      if (j > 1) {
        wj   <- W[, j]
        Etaj <- Y[, j] - Bid[, j] - wj * dthj
        Binv <- crossprod(X / wj, X)
        Bmat <- chol2inv(chol(Binv))
        mu_v <- drop(crossprod(Etaj / wj, X) %*% Bmat)
        Qj   <- Binv / Vscj
        bv   <- beta[j, ]; om2 <- omega[j, ]^2
        for (l in seq_len(K)) {
          cm    <- mu_v[l] - sum(Qj[l, -l] * (bv[-l] - mu_v[-l])) / Qj[l, l]
          bv[l] <- rtnorm(cm, 1 / sqrt(Qj[l, l]), -om2[l], om2[l])
        }
        beta[j, ] <- bv
        XB <- X %*% t(beta)
      }
      
      ## d_j (d_1 fixed; d_j updated for j >= 2 when estimate_d = TRUE)
      if (j > 1 && estimate_d) {
        jd <- j - 1   # index in log_d
        lp_dj <- function(ld) {
          dv  <- exp(ld); Vd <- dv^2 * dl2j
          er  <- Y[, j] - XB[, j] - Bid[, j] - W[, j] * dv * thj
          -nobs * ld - 0.5 * sum(er^2 / (W[, j] * Vd)) -
            0.5 * ld^2 / hyper$s2_d
        }
        hs_dj <- function(ld) {
          eps <- 1e-4
          H   <- -(lp_dj(ld + eps) - 2*lp_dj(ld) + lp_dj(ld - eps)) / eps^2
          matrix(max(H, 1e-6), 1, 1)
        }
        r_ <- mh_hess(log_d[jd], lp_dj, hs_dj, sc_d[jd])
        if (r_$acc == 1) {
          log_d[jd] <- r_$x; dvec[j] <- exp(r_$x)
          acc_d[jd] <- acc_d[jd] + 1
        }
      }
      
      ## Auxiliary variables omega_j and lambda_j of the bridge prior
      omega[j, ] <- rexp(K, lam[j, ]) + sqrt(abs(beta[j, ]))
      lam[j, ]   <- rgamma(K, 3 + hyper$a0, hyper$b0 + omega[j, ])
    }
    
    ## -------------------------------------------------------------------
    ## (F) Scalar random effects b_ij (including the survival contribution)
    ## -------------------------------------------------------------------
    ## Full conditional:
    ##   log p(b_ij | .) = -b^2 / (2 s2b_j)
    ##                     - s_ij / (2 Vsc_j) b^2 + b G_ij                 [longitudinal]
    ##                     + Delta_i psi_j b
    ##                     - exp(etaz_i + sum_{k != j} psi_k b_ik + psi_j b) H0_i  [survival]
    
    Bid <- B[id, , drop = FALSE]
    for (j in seq_len(p)) {
      dj   <- dvec[j]; thj <- theta_v[j]; dl2j <- delta2_v[j]
      Vscj <- dj^2 * dl2j; dthj <- dj * thj
      wj   <- W[, j]; psij <- psi[j]
      
      Cresj <- Y[, j] - XB[, j] - wj * dthj
      s_ij  <- as.vector(rowsum(1 / wj, id))     # length N
      Gvj   <- as.vector(rowsum(Cresj / wj, id)) / Vscj  # length N
      
      ## The hazard depends on all random effects: eta_i = z_i'gamma + sum_k psi_k b_ik.
      ## When updating b_ij, the contribution sum_{k != j} psi_k b_ik must be kept inside
      ## exp(); it can be dropped from Delta * eb because it cancels in the MH ratio.
      eb_oth <- drop(B[, -j, drop = FALSE] %*% psi[-j])
      logtbj <- function(bv) {
        eb <- psij * bv
        -0.5 * bv^2 / s2b[j] -
          0.5 * s_ij * bv^2 / Vscj +
          bv * Gvj +
          dat$Delta * eb - exp(etaz + eb_oth + eb) * H0
      }
      prop_sdj <- sc_b[j] / sqrt(1/s2b[j] + s_ij/Vscj)
      bprop    <- B[, j] + prop_sdj * rnorm(N)
      ok       <- log(runif(N)) < (logtbj(bprop) - logtbj(B[, j]))
      B[ok, j] <- bprop[ok]
      acc_b[j] <- acc_b[j] + mean(ok)
    }
    Bid <- B[id, , drop = FALSE]
    
    ## -------------------------------------------------------------------
    ## (G) Scalar random-effect variances s2b_j (independent inverse-gamma updates)
    ## -------------------------------------------------------------------
    for (j in seq_len(p)) {
      s2b[j] <- 1 / rgamma(1,
                           shape = (hyper$mb + N) / 2,
                           rate  = (hyper$Phib + sum(B[, j]^2)) / 2)
    }
    
    ## -------------------------------------------------------------------
    ## (H) gamma: coefficients of the baseline survival covariates
    ## -------------------------------------------------------------------
    H0 <- H0_fun(zeta, sp)
    eb <- drop(B %*% psi)
    lp_g <- function(gg) -0.5 * sum(gg^2 * diag(Sig_g_inv)) +
      sum(dat$Delta * drop(dat$Z %*% gg)) -
      sum(exp(drop(dat$Z %*% gg) + eb) * H0)
    hs_g <- function(gg) Sig_g_inv +
      crossprod(dat$Z * (exp(drop(dat$Z %*% gg) + eb) * H0), dat$Z)
    r_ <- mh_hess(gamma, lp_g, hs_g, sc_gamma)
    gamma <- r_$x; acc_gamma <- acc_gamma + r_$acc
    etaz  <- drop(dat$Z %*% gamma)
    
    ## -------------------------------------------------------------------
    ## (I) psi: association between random effects and the hazard
    ## -------------------------------------------------------------------
    lp_p <- function(ps) -0.5 * sum(ps^2 * diag(Sig_p_inv)) +
      sum(dat$Delta * drop(B %*% ps)) -
      sum(exp(etaz + drop(B %*% ps)) * H0)
    hs_p <- function(ps) Sig_p_inv +
      crossprod(B * (exp(etaz + drop(B %*% ps)) * H0), B)
    r_ <- mh_hess(psi, lp_p, hs_p, sc_psi)
    psi <- r_$x; acc_psi <- acc_psi + r_$acc
    
    ## -------------------------------------------------------------------
    ## (J) zeta: B-spline coefficients of the log baseline hazard
    ## -------------------------------------------------------------------
    ee   <- exp(etaz + drop(B %*% psi))
    lp_z <- function(zt) -0.5 * sum(zt^2 * diag(Sig_z_inv)) +
      sum(dat$Delta * drop(sp$CT %*% zt)) - sum(ee * H0_fun(zt, sp))
    hs_z <- function(zt) {
      wv <- sp$wgl * exp(drop(sp$Cgl %*% zt)) * rep(ee, sp$K_GL)
      Sig_z_inv + crossprod(sp$Cgl * wv, sp$Cgl)
    }
    r_ <- mh_hess(zeta, lp_z, hs_z, sc_zeta)
    zeta <- r_$x; acc_zeta <- acc_zeta + r_$acc
    
    ## -------------------------------------------------------------------
    ## Adaptive step sizes (adjusted every 50 iterations during burn-in)
    ## -------------------------------------------------------------------
    nacc <- nacc + 1
    if (it %% 50 == 0) {
      if (it <= nburn) {
        adj <- function(sc, rate) {
          if      (rate < 0.20) sc * 0.70
          else if (rate < 0.25) sc * 0.85
          else if (rate > 0.50) sc * 1.30
          else if (rate > 0.45) sc * 1.15
          else sc
        }
        sc_alpha <- pmin(pmax(adj(sc_alpha, acc_alpha / nacc), 1e-4), 10)
        sc_beta1 <- pmin(pmax(adj(sc_beta1, acc_beta1 / nacc), 1e-4), 10)
        sc_gamma <- pmin(pmax(adj(sc_gamma, acc_gamma / nacc), 1e-4), 10)
        sc_psi   <- pmin(pmax(adj(sc_psi,   acc_psi   / nacc), 1e-4), 10)
        sc_zeta  <- pmin(pmax(adj(sc_zeta,  acc_zeta  / nacc), 1e-4), 10)
        for (j in seq_len(p)) {
          sc_b[j] <- pmin(pmax(adj(sc_b[j], acc_b[j] / nacc), 1e-4), 10)
        }
        for (jd in seq_len(p - 1)) {
          sc_d[jd] <- pmin(pmax(adj(sc_d[jd], acc_d[jd] / nacc), 1e-4), 10)
        }
      }
      ## Reset acceptance counters
      acc_alpha <- 0; acc_beta1 <- 0
      acc_b[]   <- 0; acc_d[]   <- 0
      acc_gamma <- 0; acc_psi   <- 0; acc_zeta <- 0
      nacc      <- 0
    }
    
    ## -------------------------------------------------------------------
    ## Store posterior draws
    ## -------------------------------------------------------------------
    if (it %in% keep) {
      g <- g + 1
      out$beta[g, , ]  <- beta
      out$alpha[g, ]   <- alpha
      out$gamma[g, ]   <- gamma
      out$psi[g, ]     <- psi
      out$zeta[g, ]    <- zeta
      out$s2b[g, ]     <- s2b
      out$dvec[g, ]    <- dvec
    }
    if (verbose > 0 && it %% verbose == 0)
      cat(sprintf("  [Sep-AL] iter %5d\n", it))
  }
  out
}


## =========================================================================
## 5. One replicate: fit both methods to the same simulated data set
## =========================================================================
## Simulate one data set (seed = seed0 + rep) and fit both methods to it.
run_one_compare <- function(rep, tau, cfg, Psi) {
  set.seed(cfg$seed0 + rep)
  dat <- gen_data(N = cfg$N, tau = tau, u0 = cfg$u0, d1 = cfg$d1, Psi = Psi)
  p <- dat$p; K <- dat$K; nq <- ncol(dat$Z)
  tb <- dat$truth
  
  ## Method 1: Joint-MAL
  fj <- fit_jointQR(dat, niter = cfg$niter, nburn = cfg$nburn,
                    thin = cfg$thin, verbose = 0)
  J <- list(
    beta     = apply(fj$beta,  c(2,3), mean),
    beta_lo  = apply(fj$beta,  c(2,3), quantile, 0.025),
    beta_hi  = apply(fj$beta,  c(2,3), quantile, 0.975),
    gamma    = colMeans(fj$gamma),
    gamma_lo = apply(fj$gamma, 2, quantile, 0.025),
    gamma_hi = apply(fj$gamma, 2, quantile, 0.975),
    psi      = colMeans(fj$psi),
    psi_lo   = apply(fj$psi,   2, quantile, 0.025),
    psi_hi   = apply(fj$psi,   2, quantile, 0.975))
  
  ## Method 2: Sep-AL (univariate AL submodels with a shared survival submodel)
  fs <- fit_sepAL(dat, niter = cfg$niter, nburn = cfg$nburn,
                  thin = cfg$thin, verbose = 0)
  S_ <- list(
    beta     = apply(fs$beta,  c(2,3), mean),
    beta_lo  = apply(fs$beta,  c(2,3), quantile, 0.025),
    beta_hi  = apply(fs$beta,  c(2,3), quantile, 0.975),
    gamma    = colMeans(fs$gamma),
    gamma_lo = apply(fs$gamma, 2, quantile, 0.025),
    gamma_hi = apply(fs$gamma, 2, quantile, 0.975),
    psi      = colMeans(fs$psi),
    psi_lo   = apply(fs$psi,   2, quantile, 0.025),
    psi_hi   = apply(fs$psi,   2, quantile, 0.975))
  
  ## Parameter-level results in long format (for Bias / RMSE / CP)
  mk <- function(M, est, lo, hi, nm, truth) data.frame(
    method = M, param = nm, true = truth, est = est, lo = lo, hi = hi,
    stringsAsFactors = FALSE)
  nm_beta  <- as.vector(outer(seq_len(p), seq_len(K),
                              function(a, b) sprintf("beta%d%d", a, b)))
  nm_gamma <- sprintf("gamma%d", seq_len(nq))
  nm_psi   <- sprintf("psi%d",   seq_len(p))
  
  long <- rbind(
    mk("Joint", as.vector(J$beta),  as.vector(J$beta_lo),  as.vector(J$beta_hi),
       nm_beta,  as.vector(tb$beta)),
    mk("Joint", J$gamma, J$gamma_lo, J$gamma_hi, nm_gamma, tb$gamma),
    mk("Joint", J$psi,   J$psi_lo,   J$psi_hi,   nm_psi,   tb$psi),
    mk("Sep",   as.vector(S_$beta), as.vector(S_$beta_lo), as.vector(S_$beta_hi),
       nm_beta,  as.vector(tb$beta)),
    mk("Sep",   S_$gamma, S_$gamma_lo, S_$gamma_hi, nm_gamma, tb$gamma),
    mk("Sep",   S_$psi,   S_$psi_lo,   S_$psi_hi,   nm_psi,   tb$psi))
  long$rep <- rep
  
  ## Contribution of this replicate to the APMSE
  sq <- function(a, b) sum((a - b)^2)
  apm <- data.frame(
    rep = rep,
    APMSE_Joint = sq(J$beta, tb$beta) + sq(J$gamma, tb$gamma) + sq(J$psi, tb$psi),
    APMSE_Sep   = sq(S_$beta, tb$beta) + sq(S_$gamma, tb$gamma) + sq(S_$psi, tb$psi),
    beta_Joint  = sq(J$beta,  tb$beta),  beta_Sep  = sq(S_$beta,  tb$beta),
    gamma_Joint = sq(J$gamma, tb$gamma), gamma_Sep = sq(S_$gamma, tb$gamma),
    psi_Joint   = sq(J$psi,   tb$psi),   psi_Sep   = sq(S_$psi,   tb$psi),
    cens = 1 - mean(dat$Delta))
  
  list(long = long, apm = apm)
}


## =========================================================================
## 6. Summary functions
## =========================================================================
## Bias, RMSE and coverage probability (CP) of each parameter for each method.
summarise_compare <- function(long) {
  key <- paste(long$method, long$param, sep = "|")
  sp  <- split(long, factor(key, levels = unique(key)))
  do.call(rbind, lapply(sp, function(d) data.frame(
    method = d$method[1], param = d$param[1], true = d$true[1],
    mean = mean(d$est),
    Bias = mean(d$est - d$true),
    RMSE = sqrt(mean((d$est - d$true)^2)),
    CP   = mean(d$lo <= d$true & d$true <= d$hi),
    row.names = NULL)))
}

## APMSE with Monte Carlo standard error, relative efficiency (RelEff) and win rate.
summarise_apmse <- function(apm) {
  mc <- function(x) sd(x) / sqrt(length(x))
  data.frame(
    APMSE_Joint = mean(apm$APMSE_Joint), MCSE_Joint = mc(apm$APMSE_Joint),
    APMSE_Sep   = mean(apm$APMSE_Sep),   MCSE_Sep   = mc(apm$APMSE_Sep),
    RelEff      = mean(apm$APMSE_Sep) / mean(apm$APMSE_Joint),
    WinRate     = mean(apm$APMSE_Joint < apm$APMSE_Sep),
    beta_Joint  = mean(apm$beta_Joint),  beta_Sep  = mean(apm$beta_Sep),
    gamma_Joint = mean(apm$gamma_Joint), gamma_Sep = mean(apm$gamma_Sep),
    psi_Joint   = mean(apm$psi_Joint),   psi_Sep   = mean(apm$psi_Sep),
    cens        = mean(apm$cens),
    row.names = NULL)
}


## =========================================================================
## 7. Main program: N = 500, low censoring, 2 correlation levels x 3 quantile combinations x 2 methods
## =========================================================================
if (sys.nframe() == 0L || interactive()) {
  
  cfg <- list(
    N = 500, S = 100,
    niter = 20000, nburn = 5000, thin = 10,
    u0 = 6.0,           # low censoring (about 20%)
    d1 = 1,
    seed0 = 20260921,
    ncores = max(1, parallel::detectCores() - 1)
  )
  
  taus <- list(c(0.50, 0.50, 0.50),
               c(0.25, 0.25, 0.75),
               c(0.75, 0.75, 0.25))
  psis <- list(
    low  = matrix(0.3, 3, 3) + diag(0.7, 3),
    high = matrix(0.6, 3, 3) + diag(0.4, 3)
  )
  psi_labels <- list(low = "Low", high = "High")
  
  ## Objects exported to the parallel workers (all defined in this file)
  exports <- c(
    ## Utility functions
    "rtnorm", "rinvgauss_vec", "rgig_vec", "riwish", "ldmvn",
    "gauss_legendre", "u2cholcorr", "logpost_Psi",
    "theta_fun", "delta2_fun",
    ## Model functions
    "gen_data", "prep_spline", "H0_fun", "surv_loglik", "mh_hess",
    "fit_jointQR",
    ## Comparison functions
    "fit_sepAL", "run_one_compare",
    ## Configuration
    "cfg"
  )
  
  cl <- makeCluster(cfg$ncores)
  clusterEvalQ(cl, library(splines))
  clusterExport(cl, exports)
  
  all_cmp <- list()
  for (pn in names(psis)) {
    Psi  <- psis[[pn]]; plab <- psi_labels[[pn]]
    clusterExport(cl, "Psi")
    for (tau in taus) {
      lab <- sprintf("%s_tau(%s)", plab, paste(tau, collapse = ","))
      cat("Running comparison:", lab, "...\n")
      res <- parLapply(cl, seq_len(cfg$S), function(s, tau, cfg, Psi)
        tryCatch(run_one_compare(s, tau, cfg, Psi), error = function(e) NULL),
        tau = tau, cfg = cfg, Psi = Psi)
      n_ok <- sum(!vapply(res, is.null, logical(1)))
      cat(sprintf("  成功完成 %d / %d 次重复\n", n_ok, cfg$S))
      res  <- res[!vapply(res, is.null, logical(1))]
      long <- do.call(rbind, lapply(res, `[[`, "long"))
      apm  <- do.call(rbind, lapply(res, `[[`, "apm"))
      
      tab_par <- summarise_compare(long)
      tab_apm <- summarise_apmse(apm)
      tab_apm$n_ok <- n_ok
      all_cmp[[lab]] <- list(long = long, apm = apm,
                             par = tab_par, apmse = tab_apm)
      
      cat("\n-- 逐参数对比 --\n")
      print(cbind(setting = lab, format(tab_par, digits = 3)), row.names = FALSE)
      cat("\n-- APMSE --\n")
      print(cbind(setting = lab, format(tab_apm, digits = 4)), row.names = FALSE)
      cat("\n")
    }
    saveRDS(all_cmp, sprintf("compare-%s-results.rds", plab))
  }
  stopCluster(cl)
  
  ## Export summary tables
  par_all <- do.call(rbind, lapply(names(all_cmp), function(nm)
    cbind(setting = nm, all_cmp[[nm]]$par)))
  apm_all <- do.call(rbind, lapply(names(all_cmp), function(nm)
    cbind(setting = nm, all_cmp[[nm]]$apmse)))
  write.csv(par_all, "compare_two_methods_param_summary.csv", row.names = FALSE)
  write.csv(apm_all, "compare_two_methods_APMSE.csv",         row.names = FALSE)
  saveRDS(all_cmp,   "compare_two_methods_results.rds")
  cat("对比模拟全部完成。\n")
  
  print(apm_all[, c("setting", "APMSE_Joint", "APMSE_Sep", "RelEff", "WinRate")])
}