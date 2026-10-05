# Portable base-R checks: no real GWAS input or new test dependencies.
library(MESuSiE)
converged <- getFromNamespace(".mesusie_elbo_converged", "MESuSiE")
stopifnot(converged(0, .001), converged(.0001, .001),
          !converged(-.0001, .001), !converged(-1, .001),
          !converged(.001, .001), !converged(Inf, .001),
          !converged(NaN, .001))

# Stable Gaussian posterior moments agree with a direct reference, including
# singular priors (for which an inverse-prior formula cannot be used).
beta <- matrix(c(.05, -.04, .03, .01, .02, -.02), 2, 3, byrow=TRUE)
se2 <- matrix(rep(c(1e-8, 2e-8, 3e-8), each=2), 2, 3)
for (V in list(diag(c(.2, .3, .1)), matrix(.1, 3, 3), matrix(0, 3, 3))) {
  out <- mes_stable_mvlmm(beta, se2, V)
  for (j in seq_len(nrow(beta))) {
    S <- diag(se2[j, ])
    if (all(eigen(V, symmetric=TRUE)$values > 1e-10)) {
      C <- solve(solve(V) + solve(S))
      m <- C %*% (beta[j, ] / se2[j, ])
    } else {
      root <- matrix(sqrt(.1), 3, 1)
      if (all(V == 0)) root[] <- 0
      H <- solve(diag(ncol(root)) + t(root) %*% solve(S, root))
      C <- root %*% H %*% t(root)
      m <- C %*% (beta[j, ] / se2[j, ])
    }
    stopifnot(max(abs(out$post_mean[j, ] - m)) < 1e-9,
              max(abs(out$post_mean2[j, ] - (m^2 + diag(C)))) < 1e-9,
              all(out$post_mean2 >= 0), all(is.finite(out$lbf)))
  }
}

# Nonexchangeable SNP/configuration weights remain attached to identity.
# Six populations exercise the generic path beyond the optimized <=5 kernels.
set.seed(364)
for (n in c(2L, 5L, 6L)) {
  configs <- unlist(lapply(seq_len(n), function(k)
    combn(seq_len(n), k, simplify=FALSE)), recursive=FALSE)
  b <- matrix(rnorm(8*n, sd=.015), 8, n)
  s <- matrix(rep(seq(.0001, .0003, length.out=n), each=8), 8, n)
  V <- diag(seq(.001, .002, length.out=n)) + matrix(.0005, n, n)
  weights <- matrix(seq_len(8*length(configs)), 8, length(configs))
  weights <- weights / sum(weights)
  fit <- mes_covariance_em(b, s, V, weights, configs, 10L, 1e-9)
  order <- rev(seq_len(n))
  permuted <- mes_covariance_em(b[, order], s[, order], V[order, order],
    weights, lapply(configs, function(x) match(x, order)), 10L, 1e-9)
  stopifnot(max(abs(fit$V - permuted$V[order, order])) < 1e-10,
            abs(fit$objective - permuted$objective) < 1e-8,
            abs(fit$objective - mes_stable_loglik(b, s, fit$V, weights,
                                                 configs)) < 1e-6)
  limited <- mes_covariance_em(b, s, V, weights, configs, 1L, 1e-9)
  stopifnot(limited$status == 5L, limited$iterations == 1L)
}

quiet_fit <- function(...) {
  invisible(capture.output(fit <- meSuSie_core(...)))
  fit
}
sets <- function(fit) sort(vapply(fit$cs$cs, function(x)
  paste(sort(x), collapse=","), ""))
all_orders <- function(x) {
  if (length(x) == 1L) return(list(x))
  unlist(lapply(x, function(k) lapply(all_orders(setdiff(x, k)),
                                     function(rest) c(k, rest))), recursive=FALSE)
}
for (n in c(2L, 3L, 5L)) {
  p <- 16L
  ss <- lapply(seq_len(n), function(k) {
    N <- 8000 + 2000*k
    z <- c(8 + k/10, rep(.1, p-1))
    data.frame(SNP=paste0("s", seq_len(p)), Beta=z/sqrt(N),
               Se=rep(1/sqrt(N), p), Z=z, N=rep(N, p))
  })
  ld <- replicate(n, diag(p), simplify=FALSE)
  names(ss) <- names(ld) <- paste0("P", seq_len(n))
  orders <- if (n == 3L) all_orders(seq_len(n)) else
    list(seq_len(n), rev(seq_len(n)))
  fits <- lapply(orders, function(order) quiet_fit(ld[order], ss[order], L=2L,
    optim_method="em", max_iter=50L, cor_threshold=0,
    estimate_residual_variance=FALSE))
  for (fit in fits) {
    delta <- tail(diff(fit$ELBO[seq_len(fit$niter+1L)]), 1)
    stopifnot(isTRUE(fit$converged), converged(delta, .001),
              all(is.finite(fit$pip)), fit$em_calls == 2L*fit$niter,
              fit$em_maxiter <= fit$em_calls,
              max(abs(fit$pip - fits[[1]]$pip)) < 1e-6,
              identical(sets(fit), sets(fits[[1]])))
  }
  cat("EM_ORDER_PASS", n, "populations", length(orders), "orders\n")
  if (n == 2L) {
    default <- quiet_fit(ld, ss, L=1L, max_iter=1L, cor_threshold=0)
    native <- quiet_fit(ld, ss, L=1L, max_iter=1L, cor_threshold=0,
                        optim_method="optim")
    stopifnot(identical(default$pip, native$pip), !default$converged,
              is.null(default$em_calls), identical(default$V, native$V))
  }
}
cat("MESUSIE_NUMERICAL_REGRESSIONS_PASS\n")
