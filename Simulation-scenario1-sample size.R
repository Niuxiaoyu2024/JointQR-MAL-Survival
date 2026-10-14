###############################################################################
## Main simulation study for the proposed joint quantile regression model
##
## Companion code for
##   "Bayesian Joint Quantile Regression of Longitudinal Ordinal-Continuous
##    Biomarkers and Time-to-Event Outcomes" (Niu and Zhao).
##
## This script contains the data-generating mechanism, the MH-within-Gibbs
## sampler and the simulation driver used in Section 4 of the paper. The main
## program below runs the sample-size dimension (N = 300 vs. 500); the other
## dimensions are obtained by changing the corresponding settings in the main
## program (e.g. u0, Psi_fixed) or the arguments of gen_data() (e.g. psi, alpha).
##
## Model (notation as in the paper):
##   y_it = beta_tau x_it + b_i + D theta w_it + sqrt(w_it) D Sigma^{1/2} m_it
##   Sigma = Lambda Psi Lambda,  w_it ~ Exp(1),  b_i ~ N_p(0, Sigma_b)
##   y*_1it = r  iff  alpha_{r-1} < y_1it <= alpha_r
##   h(t | b_i) = exp{ zeta' C(t) + gamma' z_i + psi' b_i }  (B-spline baseline hazard)
##
## The MH-within-Gibbs steps correspond to Section 3.2 of the paper:
##   (1) b_i: random-walk MH with proposal covariance proportional to the
##       inverse precision P_i^{-1} of the longitudinal part, which does not
##       depend on b_i, so the proposal is symmetric;
##   (2) gamma, psi, zeta and log d: MH with proposal covariance given by the
##       inverse negative Hessian at the current point, with the Hastings
##       correction q(cur | prop) / q(prop | cur);
##   (3) cut-points alpha and beta_{1.}: collapsed MH of Cowles (1996), which
##       integrates out the latent scores y_1 before updating alpha and beta_{1.};
##   (4) Psi: correlation matrix (unit diagonal) with an LKJ(eta) prior,
##       updated by random-walk MH on the atanh scale of the canonical partial
##       correlations (including the Jacobian).
## Step sizes are tuned adaptively during burn-in and held fixed afterwards.
##
## Dependencies: base R, splines, parallel; GIGrvg is required only when the
## number of longitudinal outcomes p is not 1 or 3; posterior is required only
## for the (optional) convergence diagnostics in Section 5.
###############################################################################

library(splines)
library(parallel)
# library(posterior)   # Vehtari et al. (2021): rank-normalized Rhat, bulk/tail ESS

## =========================================================================
## 0. Utility functions
## =========================================================================

## Draw from N(mu, sd^2) truncated to (lo, hi), vectorized (inverse-CDF method).
rtnorm <- function(mu, sd, lo, hi) {
  n <- length(mu)
  a <- (lo - mu) / sd
  b <- (hi - mu) / sd
  flip <- a > 0
  a2 <- ifelse(flip, -b, a)
  b2 <- ifelse(flip, -a, b)
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
  n  <- length(mu)
  y  <- rnorm(n)^2
  x  <- mu + mu^2 * y / (2 * shape) -
    mu / (2 * shape) * sqrt(4 * mu * shape * y + mu^2 * y^2)
  u  <- runif(n)
  out <- ifelse(u <= mu / (mu + x), x, mu^2 / x)
  pmax(out, 1e-12)
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
  J[cbind(i, i + 1)] <- bb
  J[cbind(i + 1, i)] <- bb
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
## 1. Data generation (MAL errors)
## =========================================================================
## Simulate one data set: longitudinal ordinal-continuous outcomes with MAL errors,
## survival times with baseline hazard exp{0.8 sin(0.5 t)}, uniform censoring on
## (0, u0), and irregular visit times (exponential gaps with rate visit_rate).
gen_data <- function(N = 300,
                     tau = c(0.5, 0.5, 0.5),
                     beta  = matrix(1, 3, 3),
                     Sigma_b = diag(3),
                     Psi = matrix(0.3, 3, 3) + diag(0.7, 3),
                     d1 = 1,
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
  D      <- diag(dvec)
  Sigma  <- Lambda %*% Psi %*% Lambda
  Veps   <- D %*% Sigma %*% D
  Leps   <- t(chol(Veps))
  
  b <- matrix(rnorm(N * p), N, p) %*% chol(Sigma_b)
  z <- matrix(rnorm(N * length(gamma)), N)
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
  
  ## MAL errors
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
  list(CT   = bfun(Tobs),
       Cgl  = bfun(as.vector(S)),
       wgl  = as.vector(Wt),
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
  S0 <- scale^2 * chol2inv(chol(hess(cur)))
  prop <- cur + drop(t(chol(S0)) %*% rnorm(length(cur)))
  H1 <- tryCatch(hess(prop), error = function(e) NULL)
  if (is.null(H1) || any(!is.finite(H1))) return(list(x = cur, acc = 0))
  S1 <- scale^2 * chol2inv(chol(H1))
  lr <- logpost(prop) - logpost(cur) + ldmvn(cur, prop, S1) - ldmvn(prop, cur, S0)
  if (is.finite(lr) && log(runif(1)) < lr) list(x = prop, acc = 1)
  else list(x = cur, acc = 0)
}

## =========================================================================
## 3. MH-within-Gibbs sampler for the joint model
## =========================================================================
## Fit the joint model by MH-within-Gibbs. Returns post-burn-in, thinned draws of
## beta, alpha, gamma, psi, zeta, Sigma_b, Psi and d, full trajectories for trace
## plots, acceptance rates and the posterior mean of the random effects.
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
  ## Identification: d_1 is fixed at dat$d1 and the ordinal submodel has no intercept;
  ## all cut-points alpha_1, ..., alpha_{R-1} are free within (alpha_min, alpha_max).
  ab <- c(alpha[1] - s_y, alpha, alpha[R - 1] + s_y)
  y1 <- (ab[dat$y1star] + ab[dat$y1star + 1]) / 2
  Y  <- cbind(y1, dat$Y2)
  
  beta <- t(qr.solve(X, Y))
  omega <- sqrt(abs(beta)) + 1
  lam   <- matrix(1, p, K)
  Psi   <- diag(p)
  uPsi   <- rep(0, p * (p - 1) / 2)
  sPsi   <- 1.2 / sqrt(nobs)
  acc_Psi <- 0
  log_d  <- log(dvec[-1])
  acc_d  <- 0
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
  
  keep     <- seq(nburn + thin, niter, by = thin)  # posterior draws (after burn-in)
  keep_all <- seq(thin, niter, by = thin)           # full trajectory (including burn-in)
  b_sum <- matrix(0, N, p)
  G <- length(keep); g <- 0
  G_all <- length(keep_all); g_all <- 0
  out <- list(beta = array(NA, c(G, p, K)), alpha = matrix(NA, G, R - 1),
              gamma = matrix(NA, G, nq), psi = matrix(NA, G, p),
              zeta = matrix(NA, G, L), Sigma_b = array(NA, c(G, p, p)),
              Psi = array(NA, c(G, p, p)),
              dvec = matrix(NA, G, p),
              ## Full trajectories (including burn-in) for trace plots
              all_iter   = keep_all,
              all_alpha  = matrix(NA, G_all, R - 1),
              all_beta1  = matrix(NA, G_all, K),
              all_beta2  = matrix(NA, G_all, K),
              all_beta3  = matrix(NA, G_all, K),
              all_gamma  = matrix(NA, G_all, nq),
              all_psi    = matrix(NA, G_all, p),
              all_Psi    = array(NA, c(G_all, p, p)),
              all_d      = matrix(NA, G_all, p))
  acc_rate <- matrix(NA, 0, 7,
                     dimnames = list(NULL, c("b", "gamma", "psi", "zeta", "alpha", "d", "beta1")))
  
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
      ## Collapsed MH (Cowles, 1996): update alpha and beta_{1.} with the latent
      ## scores y_1it integrated out. Target: prod_{it} pi*_{it,r}, where
      ##   pi*_{itr} = Phi((alpha_r - m_{1|2,it}) / sd1) - Phi((alpha_{r-1} - m_{1|2,it}) / sd1),
      ## m1 is the conditional mean of y_1 given the concurrent continuous responses and
      ## sd1 the conditional SD. Both proposals are symmetric random walks, so no Jacobian
      ## or Hastings correction is needed.
      ll_ord <- function(al, m1_use) {
        cc <- c(-Inf, al, Inf)
        pr <- pnorm((cc[dat$y1star + 1] - m1_use) / sd1) -
          pnorm((cc[dat$y1star]     - m1_use) / sd1)
        sum(log(pmax(pr, 1e-300)))
      }
      
      ## --- (i) Update alpha (beta_{1.} fixed) ---
      ## Symmetric random walk on the gaps g_r = alpha_r - alpha_{r-1} (alpha_0 = alpha_min);
      ## the linear reparameterization has unit Jacobian.
      g_cur  <- diff(c(amin, alpha))           # length R - 1
      g_prop <- g_cur + rnorm(R - 1) * sc$alpha
      if (all(g_prop > 0)) {
        ap <- amin + cumsum(g_prop)
        if (ap[R - 1] < amax) {
          lr <- ll_ord(ap, m1) - ll_ord(alpha, m1)
          if (is.finite(lr) && log(runif(1)) < lr) {
            alpha <- ap
            acc$alpha <- acc$alpha + 1
          }
        }
      }
      
      ## --- (ii) Update beta_{1.} (alpha fixed) ---
      ## Uses the category probabilities directly rather than the current y_1it;
      ## a change delta in beta_{1.} shifts the conditional mean by X %*% delta.
      beta1_cur  <- beta[1, ]
      beta1_prop <- beta1_cur + rnorm(K) * sc$beta1
      if (all(abs(beta1_prop) < omega[1, ]^2)) {        # bridge-prior constraint
        m1_prop <- m1 + drop(X %*% (beta1_prop - beta1_cur))
        lr <- ll_ord(alpha, m1_prop) - ll_ord(alpha, m1)
        if (is.finite(lr) && log(runif(1)) < lr) {
          beta[1, ] <- beta1_prop
          m1        <- m1_prop
          XB        <- X %*% t(beta)
          acc$beta1 <- acc$beta1 + 1
        }
      }
      
      cuts <- c(-Inf, alpha, Inf)
      Y[, 1] <- rtnorm(m1, sd1, cuts[dat$y1star], cuts[dat$y1star + 1])
    }
    
    
    
    E   <- Y - XB - Bid
    chi <- rowSums((E %*% Vinv) * E)
    psiG <- drop(t(theta) %*% Sinv %*% theta) + 2
    w   <- rgig_vec(1 - p / 2, chi, psiG)
    
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
    ## beta[1, ] has already been updated by the collapsed MH step above, so only
    ## beta[2:p, ] is updated here. as.vector() stacks by column: beta[j, k] has index
    ## (k - 1) * p + j, so beta[1, ] corresponds to seq(1, p * K, by = p).
    idx_skip <- if (alpha_update == "cowles") seq(1, p * K, by = p) else integer(0)
    for (l in setdiff(seq_along(bv), idx_skip)) {
      cm <- mu_v[l] - sum(Q[l, -l] * (bv[-l] - mu_v[-l])) / Q[l, l]
      bv[l] <- rtnorm(cm, 1 / sqrt(Q[l, l]), -om2[l], om2[l])
    }
    beta <- matrix(bv, p, K)
    XB   <- X %*% t(beta)
    
    Lres <- (Y - XB - Bid - outer(w, Dtheta)) / sqrt(w)
    U    <- sweep(Lres, 2, DL_inv, "*")
    S_r  <- crossprod(U)
    lp_cur <- logpost_Psi(uPsi, S_r, nobs, hyper$eta, p)
    for (rep_ in 1:5) {
      ucand <- uPsi + sPsi * rnorm(length(uPsi))
      lp_new <- logpost_Psi(ucand, S_r, nobs, hyper$eta, p)
      if (log(runif(1)) < lp_new - lp_cur) {
        uPsi <- ucand; lp_cur <- lp_new; acc_Psi <- acc_Psi + 1
      }
    }
    Lc  <- u2cholcorr(uPsi, p)$L
    Psi <- Lc %*% t(Lc)
    Sigma <- Lambda %*% Psi %*% Lambda   # refresh Sigma for the update of d below
    
    if (estimate_d && p > 1) {
      lp_d <- function(ld) {
        dv  <- c(d1fix, exp(ld)); Dv <- diag(dv)
        Vd  <- Dv %*% Sigma %*% Dv
        Vdi <- chol2inv(chol(Vd))
        Er  <- Y - XB - Bid - outer(w, drop(Dv %*% theta))
        ll  <- -nobs * sum(log(dv[-1])) -
          0.5 * sum(rowSums((Er / w) * (Er %*% Vdi)))
        ll - 0.5 * sum(ld^2) / hyper$s2_d
      }
      hs_d <- function(ld) {
        eps <- 1e-4
        H   <- numeric(p - 1)
        lp0 <- lp_d(ld)
        for (j in seq_len(p - 1)) {
          ldp <- ld; ldp[j] <- ld[j] + eps
          ldm <- ld; ldm[j] <- ld[j] - eps
          H[j] <- -(lp_d(ldp) - 2 * lp0 + lp_d(ldm)) / eps^2
        }
        diag(pmax(H, 1e-6))
      }
      r_ <- mh_hess(log_d, lp_d, hs_d, sc$d)
      if (r_$acc == 1) {
        log_d  <- r_$x; acc$d <- acc$d + 1
        dvec   <- c(d1fix, exp(log_d))
        D      <- diag(dvec)
        Dtheta <- drop(D %*% theta)
        DL_inv <- 1 / (dvec * sqrt(delta2_fun(tau)))
        V      <- D %*% Sigma %*% D
        Vinv   <- chol2inv(chol(V))
      }
    }
    
    omega <- matrix(rexp(p * K, as.vector(lam)), p, K) + sqrt(abs(beta))
    lam   <- matrix(rgamma(p * K, 3 + hyper$a0, hyper$b0 + as.vector(omega)), p, K)
    
    if (alpha_update == "gibbs") for (r in seq_len(R - 1)) {  # all cut-points free
      lo <- max(c(Y[dat$y1star == r, 1], if (r > 1) alpha[r - 1], hyper$alpha_min))
      hi <- min(c(Y[dat$y1star == r + 1, 1], if (r < R - 1) alpha[r + 1], hyper$alpha_max))
      if (hi > lo) alpha[r] <- runif(1, lo, hi)
    }
    
    Sigma_b <- riwish(N + hyper$mb, crossprod(B) + hyper$Phib)
    
    eb <- drop(B %*% psi)
    lp_g <- function(g) -0.5 * sum(g^2 * diag(Sig_g_inv)) +
      sum(dat$Delta * drop(dat$Z %*% g)) - sum(exp(drop(dat$Z %*% g) + eb) * H0)
    hs_g <- function(g) Sig_g_inv +
      crossprod(dat$Z * (exp(drop(dat$Z %*% g) + eb) * H0), dat$Z)
    r_ <- mh_hess(gamma, lp_g, hs_g, sc$gamma); gamma <- r_$x; acc$gamma <- acc$gamma + r_$acc
    etaz <- drop(dat$Z %*% gamma)
    
    lp_p <- function(ps) -0.5 * sum(ps^2 * diag(Sig_p_inv)) +
      sum(dat$Delta * drop(B %*% ps)) - sum(exp(etaz + drop(B %*% ps)) * H0)
    hs_p <- function(ps) Sig_p_inv +
      crossprod(B * (exp(etaz + drop(B %*% ps)) * H0), B)
    r_ <- mh_hess(psi, lp_p, hs_p, sc$psi); psi <- r_$x; acc$psi <- acc$psi + r_$acc
    
    ee <- exp(etaz + drop(B %*% psi))
    lp_z <- function(zt) -0.5 * sum(zt^2 * diag(Sig_z_inv)) +
      sum(dat$Delta * drop(sp$CT %*% zt)) - sum(ee * H0_fun(zt, sp))
    hs_z <- function(zt) {
      wv <- sp$wgl * exp(drop(sp$Cgl %*% zt)) * rep(ee, sp$K_GL)
      Sig_z_inv + crossprod(sp$Cgl * wv, sp$Cgl)
    }
    r_ <- mh_hess(zeta, lp_z, hs_z, sc$zeta); zeta <- r_$x; acc$zeta <- acc$zeta + r_$acc
    
    nacc <- nacc + 1
    if (it %% 50 == 0) {
      ar  <- unlist(acc) / nacc
      acc_rate <- rbind(acc_rate, ar)
      acc[] <- 0; nacc <- 0
      
      if (it <= nburn) {
        for (nm in names(sc)) {
          if (nm == "alpha" && alpha_update == "gibbs") next
          if (nm == "beta1" && alpha_update == "gibbs") next
          if (nm == "d"     && !estimate_d)             next
          if (ar[nm] < 0.20) sc[[nm]] <- sc[[nm]] * 0.70
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
        
        ## Step sizes are frozen after burn-in.
      }
    }
    
    if (it %in% keep_all) {
      g_all <- g_all + 1
      out$all_alpha[g_all, ]    <- alpha
      out$all_beta1[g_all, ]    <- beta[1, ]
      out$all_beta2[g_all, ]    <- beta[2, ]
      out$all_beta3[g_all, ]    <- beta[3, ]
      out$all_gamma[g_all, ]    <- gamma
      out$all_psi[g_all, ]      <- psi
      out$all_Psi[g_all, , ]    <- Psi
      out$all_d[g_all, ]        <- dvec
    }
    if (it %in% keep) {
      g <- g + 1
      out$beta[g, , ] <- beta; out$alpha[g, ] <- alpha
      out$gamma[g, ] <- gamma; out$psi[g, ] <- psi; out$zeta[g, ] <- zeta
      out$Sigma_b[g, , ] <- Sigma_b; out$Psi[g, , ] <- Psi
      out$dvec[g, ] <- dvec
      b_sum <- b_sum + B
    }
    if (verbose > 0 && it %% verbose == 0)
      cat(sprintf("iter %5d | acc b=%.2f g=%.2f p=%.2f z=%.2f a=%.2f d=%.2f b1=%.2f | beta1=%s\n",
                  it, acc_rate[nrow(acc_rate), 1], acc_rate[nrow(acc_rate), 2],
                  acc_rate[nrow(acc_rate), 3], acc_rate[nrow(acc_rate), 4],
                  acc_rate[nrow(acc_rate), 5], acc_rate[nrow(acc_rate), 6],
                  acc_rate[nrow(acc_rate), 7],
                  paste(round(beta[1, ], 3), collapse = ",")))
  }
  out$acc_rate <- acc_rate
  out$b_mean <- b_sum / G
  ## Acceptance rate of Psi over the post-burn-in iterations
  out$acc_Psi <- acc_Psi / (5 * (niter - nburn))
  out$spline <- sp
  out
}

## =========================================================================
## 4. Posterior summaries and one simulation replicate
## =========================================================================
## Matrix of posterior draws (columns = reported parameters).
param_draws <- function(fit) {
  G <- dim(fit$beta)[1]; p <- dim(fit$beta)[2]; K <- dim(fit$beta)[3]
  Rm1 <- ncol(fit$alpha)
  dr <- list()
  ## beta_{1k} on the original scale
  for (k in seq_len(K))
    dr[[sprintf("beta1%d", k)]] <- fit$beta[, 1, k]
  for (r in seq_len(Rm1)) dr[[sprintf("alpha%d", r)]] <- fit$alpha[, r]
  for (j in 2:p) for (k in seq_len(K)) dr[[sprintf("beta%d%d", j, k)]] <- fit$beta[, j, k]
  for (k in seq_len(ncol(fit$gamma))) dr[[sprintf("gamma%d", k)]] <- fit$gamma[, k]
  for (k in seq_len(ncol(fit$psi)))   dr[[sprintf("psi%d", k)]]   <- fit$psi[, k]
  ## psi_1 on the original scale
  dr[["psi1"]] <- fit$psi[, 1]
  for (j in 2:p) for (k in 1:(j - 1)) dr[[sprintf("Psi%d%d", j, k)]] <- fit$Psi[, j, k]
  if (!is.null(fit$dvec))
    for (j in 2:p) dr[[sprintf("d%d", j)]] <- fit$dvec[, j]
  do.call(cbind, dr)
}

## True parameter values in the same order as param_draws().
true_vector <- function(truth) {
  p <- nrow(truth$beta); K <- ncol(truth$beta); Rm1 <- length(truth$alpha)
  c(truth$beta[1, ],                                   # beta_{1k}
    truth$alpha,                                       # alpha_1, ..., alpha_{R-1}
    as.vector(t(truth$beta[-1, , drop = FALSE])),
    truth$gamma,
    truth$psi,
    t(truth$Psi)[upper.tri(truth$Psi)],
    truth$dvec[-1])
}

## Simulate one data set (seed = seed0 + rep), fit the model and return posterior
## means and 95% equal-tail credible intervals of all reported parameters.
run_one <- function(rep, tau, cfg, Psi) {
  set.seed(cfg$seed0 + rep)
  dat <- gen_data(N = cfg$N, tau = tau, u0 = cfg$u0, d1 = cfg$d1, Psi = Psi)
  fit <- fit_jointQR(dat, niter = cfg$niter, nburn = cfg$nburn, thin = cfg$thin,
                     verbose = 0)
  dr  <- param_draws(fit)
  data.frame(rep = rep, param = as.character(colnames(dr)),
             true = true_vector(dat$truth),
             est  = colMeans(dr),
             lo   = apply(dr, 2, quantile, 0.025),
             hi   = apply(dr, 2, quantile, 0.975),
             cens = 1 - mean(dat$Delta),
             nobs_per_subj = length(dat$id) / dat$N,
             row.names = NULL)
}

## Bias, RMSE and coverage probability (CP) of each parameter over the replicates.
summarise_sim <- function(res) {
  res$param <- as.character(res$param)
  sp <- split(res, factor(res$param, levels = unique(res$param)))
  do.call(rbind, lapply(sp, function(d) data.frame(
    param = d$param[1], true = d$true[1],
    mean  = mean(d$est),
    Bias  = mean(d$est - d$true),
    RMSE  = sqrt(mean((d$est - d$true)^2)),
    CP    = mean(d$lo <= d$true & d$true <= d$hi),
    CI_lo = mean(d$lo), CI_hi = mean(d$hi))))
}

## =========================================================================
## 6. Main program
## =========================================================================
if (sys.nframe() == 0L || interactive()) {
  
  cfg <- list(
    ## ---- Sample size and number of replicates ----
    N = 300,          # default value; overwritten by cfg_n$N (300 or 500) in the loop
    S = 100,        # number of replicates
    
    ## ---- MCMC settings ----
    #niter = 10000, nburn = 3000, thin = 10,
    niter = 20000, nburn = 5000, thin = 10,
    # thin = 10 retains (niter - nburn) / 10 posterior draws
    
    ## ---- Censoring rate (controlled by u0) ----
    ## Censoring times C_i ~ Uniform(0, u0); a larger u0 gives less censoring.
    u0 = 6.0,       # low censoring: about 20%
    #u0 = 1.5,       # high censoring: about 50%
    
    ## ---- Scale of the latent ordinal outcome ----
    ## d1 is the true value of d_1, which is also held fixed at this value in
    ## fit_jointQR for identification; d_2 and d_3 are estimated (estimate_d = TRUE).
    d1 = 1,          # d_1 = 1, as in the paper
    
    seed0 = 20260921,
    ncores = max(1, detectCores() - 1)
  )
  
  ## ---------- Quick check on a single data set ----------
  set.seed(1)
  dat <- gen_data(N = cfg$N, tau = c(0.5, 0.5, 0.5), u0 = cfg$u0, d1 = cfg$d1)
  cat(sprintf("删失率 = %.3f, 平均观测次数 = %.2f\n",
              1 - mean(dat$Delta), length(dat$id) / dat$N))
  t0 <- Sys.time()
  fit <- fit_jointQR(dat, niter = 20000, nburn = 5000, thin = 10, verbose = 2500)
  print(Sys.time() - t0)
  dr <- param_draws(fit)
  print(round(cbind(true = true_vector(dat$truth), est = colMeans(dr),
                    lo = apply(dr, 2, quantile, .025),
                    hi = apply(dr, 2, quantile, .975)), 3))
  
  ## ---------- Simulation: 2 sample sizes x 3 quantile combinations ----------
  taus <- list(c(0.50, 0.50, 0.50), c(0.25, 0.25, 0.75), c(0.75, 0.75, 0.25))
  
  ## Fixed: low longitudinal correlation (off-diagonal 0.3) and low censoring (u0 = 6.0)
  Psi_fixed <- matrix(0.3, 3, 3) + diag(0.7, 3)
  
  Ns        <- list(N300 = 300, N500 = 500)
  N_labels  <- list(N300 = "N=300", N500 = "N=500")
  
  cl <- makeCluster(cfg$ncores)
  clusterEvalQ(cl, library(splines))
  clusterExport(cl, c("rtnorm", "rinvgauss_vec", "rgig_vec", "riwish", "ldmvn",
                      "u2cholcorr", "logpost_Psi",
                      "gauss_legendre", "theta_fun", "delta2_fun", "gen_data",
                      "prep_spline", "H0_fun", "surv_loglik", "mh_hess",
                      "fit_jointQR", "param_draws", "true_vector", "run_one",
                      "cfg", "Psi_fixed"))
  
  all_res <- list()
  for (n_name in names(Ns)) {
    N_val   <- Ns[[n_name]]
    N_lab   <- N_labels[[n_name]]
    ## Local copy of cfg for the current sample size (only N changes)
    cfg_n   <- cfg
    cfg_n$N <- N_val
    clusterExport(cl, "cfg_n")                  ## update cfg_n on the workers
    
    for (tau in taus) {
      tau_str <- paste(tau, collapse = ",")
      lab     <- sprintf("MAL_%s_tau(%s)", N_lab, tau_str)
      cat("Running", lab, "...\n")
      res <- parLapply(cl, seq_len(cfg$S), function(s, tau, cfg_n, Psi_fixed)
        tryCatch(run_one(s, tau, cfg_n, Psi_fixed), error = function(e) NULL),
        tau = tau, cfg_n = cfg_n, Psi_fixed = Psi_fixed)
      res <- do.call(rbind, res)
      tab <- summarise_sim(res)
      all_res[[lab]] <- list(raw = res, table = tab)
      print(cbind(setting = lab, format(tab, digits = 3)), row.names = FALSE)
    }
    
    ## Save intermediate results after each sample size
    out_rds <- sprintf("低删失-低相关-%s-sim_results.rds", N_lab)
    saveRDS(all_res, out_rds)
    cat("Saved intermediate results to", out_rds, "\n")
  }
  stopCluster(cl)
  
  final <- do.call(rbind, lapply(names(all_res), function(nm)
    cbind(setting = nm, all_res[[nm]]$table)))
  write.csv(final, "低删失-低相关-两种N-sim_summary.csv", row.names = FALSE)
  saveRDS(all_res, "低删失-低相关-两种N-sim_results.rds")
  cat("全部完成，最终结果已保存。\n")
}