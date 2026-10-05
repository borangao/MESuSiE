# Optional covariance estimation; the author optimizer remains the default.
.mesusie_elbo_converged <- function(delta, tol) {
  is.finite(delta) && delta >= 0 && delta < tol
}

.mesusie_update_covariance_em <- function(beta, se2, obj, effect) {
  old <- obj$V[[effect]]
  have <- all(diag(old) > 0)
  initial <- if (have) old else
    diag(pmax(1e-10, apply(beta^2 - se2, 2, max)), ncol(beta))
  fit <- mes_covariance_em(beta, se2, initial, obj$pi,
                           obj$column_config, obj$em_max_iter, obj$em_tol)
  value <- mes_stable_loglik(beta, se2, fit$V, obj$pi, obj$column_config)
  if (!is.finite(value) || abs(value - fit$objective) >= 1e-6)
    stop("Covariance EM objective disagrees with the Gaussian likelihood")
  oldvalue <- if (have)
    mes_stable_loglik(beta, se2, old, obj$pi, obj$column_config) else Inf
  fallback <- is.finite(oldvalue) && oldvalue < value
  obj$em_calls <- obj$em_calls + 1L
  obj$em_maxiter <- obj$em_maxiter + as.integer(fit$status == 5L)
  obj$em_fallbacks <- obj$em_fallbacks + as.integer(fallback)
  obj$em_status[[effect]] <- list(status = fit$status,
                                 iterations = fit$iterations,
                                 objective = value, fallback = fallback)
  if (fallback) old else fit$V
}
