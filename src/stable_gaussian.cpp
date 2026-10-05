// Adapted from credtools 0.9.8; Copyright (c) 2024, Jianhua Wang.
// Original MIT license retained in inst/CREDTOOLS_LICENSE.
// Numerically stable evaluations of the unchanged Gaussian MESuSiE model.
#include <RcppArmadillo.h>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
void validate(const arma::mat& beta, const arma::mat& se2, const arma::mat& V) {
  if (beta.is_empty() || se2.n_rows != beta.n_rows ||
      se2.n_cols != beta.n_cols || V.n_rows != beta.n_cols ||
      V.n_cols != beta.n_cols || !beta.is_finite() || !se2.is_finite() ||
      !V.is_finite() || arma::any(arma::vectorise(se2) <= 0))
    Rcpp::stop("Invalid MESuSiE Gaussian dimensions or data");
}

arma::mat prior_root(const arma::mat& V) {
  arma::vec eig;
  arma::mat U;
  const double scale = arma::norm(V, "inf");
  const double roundoff = 64 * std::numeric_limits<double>::epsilon() *
                          V.n_rows * scale;
  if (arma::norm(V - V.t(), "inf") > roundoff ||
      !arma::eig_sym(eig, U, arma::symmatu(V)) || eig.min() < -roundoff)
    Rcpp::stop("MESuSiE prior covariance is not symmetric positive semidefinite");
  // Clip only eigensolver-scale negative roundoff, not substantive negatives.
  eig.transform([](double x) { return std::sqrt(std::max(0.0, x)); });
  return U * arma::diagmat(eig);
}

Rcpp::List gaussian(const arma::mat& beta, const arma::mat& se2,
                    const arma::mat& V, bool posterior) {
  validate(beta, se2, V);
  const arma::mat root = prior_root(V);
  const unsigned p = beta.n_rows, n = beta.n_cols;
  arma::vec lbf(p);
  arma::mat mean(p, n, arma::fill::zeros), second(p, n, arma::fill::zeros);
  for (unsigned j = 0; j < p; ++j) {
    const arma::vec b = beta.row(j).t(), s = se2.row(j).t();
    arma::mat chol;
    if (!arma::chol(chol, V + arma::diagmat(s), "lower"))
      Rcpp::stop("Non-positive MESuSiE sampling covariance");
    const arma::vec u = arma::solve(arma::trimatl(chol), b,
                                  arma::solve_opts::fast);
    lbf(j) = .5 * (arma::accu(arma::square(b) / s) - arma::dot(u, u) +
                      arma::accu(arma::log(s))) - arma::accu(arma::log(chol.diag()));
    if (posterior) {
      // Latent square-root conditioning avoids S - inv(W + W V W),
      // whose cancellation can make posterior variances negative.
      arma::mat weighted = root;
      weighted.each_col() /= arma::sqrt(s);
      arma::mat H;
      if (!arma::chol(H, arma::eye(n, n) + weighted.t() * weighted, "lower"))
        Rcpp::stop("Non-positive MESuSiE latent precision");
      const arma::mat X = arma::solve(arma::trimatl(H), root.t(),
                                      arma::solve_opts::fast);
      const arma::vec t = arma::solve(arma::trimatl(H), root.t() * (b / s),
                                      arma::solve_opts::fast);
      const arma::vec m = X.t() * t;
      mean.row(j) = m.t();
      second.row(j) = (arma::square(m) + arma::sum(arma::square(X), 0).t()).t();
    }
  }
  if (!lbf.is_finite() || !mean.is_finite() || !second.is_finite())
    Rcpp::stop("Nonfinite MESuSiE Gaussian evaluation");
  return Rcpp::List::create(Rcpp::Named("lbf") = lbf,
                            Rcpp::Named("post_mean") = mean,
                            Rcpp::Named("post_mean2") = second);
}
}

// [[Rcpp::export]]
Rcpp::List mes_stable_mvlmm(const arma::mat& beta, const arma::mat& se2,
                          const arma::mat& V) {
  return gaussian(beta, se2, V, true);
}

// [[Rcpp::export]]
double mes_stable_loglik(const arma::mat& beta, const arma::mat& se2,
                        const arma::mat& V, const arma::mat& prior,
                        Rcpp::List configs) {
  validate(beta, se2, V);
  if (configs.size() < 1 || prior.n_rows != beta.n_rows ||
      prior.n_cols != static_cast<unsigned>(configs.size()) ||
      !prior.is_finite() || arma::any(arma::vectorise(prior) < 0) ||
      arma::accu(prior) <= 0)
    Rcpp::stop("Invalid MESuSiE configuration weights");
  prior_root(V);
  double shift = -arma::datum::inf, mass = 0;
  for (int k = 0; k < configs.size(); ++k) {
    arma::uvec active = Rcpp::as<arma::uvec>(configs[k]);
    const arma::uvec unique = arma::unique(active);
    if (active.is_empty() || active.min() < 1 || active.max() > beta.n_cols ||
        unique.n_elem != active.n_elem)
      Rcpp::stop("Invalid MESuSiE ancestry configuration");
    active -= 1;
    const Rcpp::List result = gaussian(beta.cols(active), se2.cols(active),
                                        V.submat(active, active), false);
    const arma::vec lbf = Rcpp::as<arma::vec>(result["lbf"]);
    for (unsigned j = 0; j < beta.n_rows; ++j) {
      if (prior(j, k) == 0) continue;
      const double score = lbf(j) + std::log(prior(j, k));
      if (score > shift) {
        mass *= std::exp(shift - score);
        shift = score;
      }
      mass += std::exp(score - shift);
    }
  }
  const double value = -(shift + std::log(mass));
  if (!std::isfinite(value)) Rcpp::stop("Nonfinite MESuSiE likelihood");
  return value;
}
