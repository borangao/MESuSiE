// Adapted from credtools 0.9.8; Copyright (c) 2024, Jianhua Wang.
// Original MIT license retained in inst/CREDTOOLS_LICENSE.
// Exact-model covariance EM: cache invariant data and reuse one Cholesky solve.
// No SNP pruning, precision reduction, relaxed tolerance or changed EM budget.
#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {
struct ConfigData {
  arma::uvec active;
  std::vector<double> y, s, null_quad, null_logdet;
};

struct Accumulator {
  arma::mat H;
  double shift=-arma::datum::inf, mass=0;
  explicit Accumulator(int d):H(d,d,arma::fill::zeros) {}
  double weight(double lbf,double prior) {
    if(lbf>shift) {
      const double rescale=std::exp(shift-lbf);
      H*=rescale;mass*=rescale;shift=lbf;
    }
    const double w=std::exp(lbf-shift)*prior;
    mass+=w;
    return w;
  }
};
// Stream locally log-scaled sufficient moments; no SNP moment cube.
template<int D>
void moments(const ConfigData& data, const arma::mat& V, int p, int k,
             const arma::mat& prior, std::vector<Accumulator>& stats) {
  double VV[D*D];
  for(int a=0;a<D;++a) for(int b=0;b<D;++b)
    VV[a*D+b]=V(data.active[a],data.active[b]);
  for(int j=0;j<p;++j) {
    double chol[D*D]={}, u[D]={}, inverseL[D*D]={}, t[D]={};
    double logdet=0,quad=0;
    for(int a=0;a<D;++a) {
      for(int b=0;b<=a;++b) {
        double v=VV[a*D+b]+(a==b ? data.s[j*D+a] : 0.0);
        for(int c=0;c<b;++c) v-=chol[a*D+c]*chol[b*D+c];
        if(a==b) {
          if(!std::isfinite(v) || v<=0) Rcpp::stop("Non-positive EM sampling covariance");
          chol[a*D+b]=std::sqrt(v);
          logdet+=2*std::log(chol[a*D+b]);
        } else chol[a*D+b]=v/chol[b*D+b];
      }
      double v=data.y[j*D+a];
      for(int b=0;b<a;++b) v-=chol[a*D+b]*u[b];
      u[a]=v/chol[a*D+a];quad+=u[a]*u[a];
      for(int c=0;c<=a;++c) {
        v=(a==c ? 1.0 : 0.0);
        for(int b=c;b<a;++b) v-=chol[a*D+b]*inverseL[b*D+c];
        inverseL[a*D+c]=v/chol[a*D+a];
      }
    }
    for(int a=0;a<D;++a) for(int b=a;b<D;++b) t[a]+=inverseL[b*D+a]*u[b];
    const double score=.5*(data.null_quad[j]-quad+data.null_logdet[j]-logdet);
    const double weight=stats[k].weight(score,prior(j,k));
    double* dst=stats[k].H.memptr();
    for(int b=0;b<D;++b) for(int a=0;a<=b;++a) {
      double precision=0;
      for(int c=b;c<D;++c) precision+=inverseL[c*D+a]*inverseL[c*D+b];
      const double value=t[a]*t[b]-precision;
      dst[a+b*D]+=weight*value;
      if(a!=b)dst[b+a*D]=dst[a+b*D];
    }
  }
}

void dispatch(const ConfigData& data, const arma::mat& V, int p, int k,
              const arma::mat& prior, std::vector<Accumulator>& stats) {
  switch(data.active.n_elem) {
    case 1: moments<1>(data,V,p,k,prior,stats);break;
    case 2: moments<2>(data,V,p,k,prior,stats);break;
    case 3: moments<3>(data,V,p,k,prior,stats);break;
    case 4: moments<4>(data,V,p,k,prior,stats);break;
    case 5: moments<5>(data,V,p,k,prior,stats);break;
    default: Rcpp::stop("Unsupported fast-path configuration");
  }
}
}

// [[Rcpp::export]]
Rcpp::List mes_covariance_em(const arma::mat& beta, const arma::mat& se2,
                            arma::mat V, const arma::mat& prior,
                            Rcpp::List configs, int maxiter=100, double tol=1e-9) {
  const int p=beta.n_rows,n=beta.n_cols,nc=configs.size();
  if(p<1 || n<1 || nc<1 || maxiter<1 || !std::isfinite(tol) || tol<=0 ||
     se2.n_rows!=beta.n_rows || se2.n_cols!=beta.n_cols ||
     V.n_rows!=beta.n_cols || V.n_cols!=beta.n_cols ||
     prior.n_rows!=beta.n_rows || prior.n_cols!=(unsigned)nc ||
     !beta.is_finite() || !se2.is_finite() || !V.is_finite() || !prior.is_finite() ||
     arma::any(arma::vectorise(se2)<=0) || arma::any(arma::vectorise(prior)<0) || arma::accu(prior)<=0)
    Rcpp::stop("Invalid covariance-EM dimensions, data, weights or controls");
  std::vector<ConfigData> cached(nc);
  for(int k=0;k<nc;++k) {
    ConfigData& data=cached[k];
    data.active=Rcpp::as<arma::uvec>(configs[k]);
    if(data.active.is_empty() || data.active.min()<1 || data.active.max()>(unsigned)n)
      Rcpp::stop("Invalid ancestry configuration");
    data.active-=1;
    const int d=data.active.n_elem;
    data.y.resize(p*d);data.s.resize(p*d);
    data.null_quad.resize(p);data.null_logdet.resize(p);
    for(int j=0;j<p;++j) for(int a=0;a<d;++a) {
      const double y=beta(j,data.active[a]),s=se2(j,data.active[a]);
      data.y[j*d+a]=y;data.s[j*d+a]=s;
      data.null_quad[j]+=y*y/s;data.null_logdet[j]+=std::log(s);
    }
  }
  double previous=-arma::datum::inf,objective=previous;
  int iteration=0,status=5;
  std::vector<Accumulator> stats;
  for(const auto& data:cached) stats.emplace_back(data.active.n_elem);
  for(iteration=0;iteration<maxiter;++iteration) {
    Rcpp::checkUserInterrupt();
    for(auto& item:stats) {item.H.zeros();item.mass=0;item.shift=-arma::datum::inf;}
    for(int k=0;k<nc;++k) {
      const ConfigData& data=cached[k];
      if(n<=5) dispatch(data,V,p,k,prior,stats);
      else {
          // Use generic Armadillo solves beyond the specialized small dimensions.
          const arma::uvec& active=data.active;
          arma::mat VA=V.cols(active),VV=V.submat(active,active);
          for(int j=0;j<p;++j) {
            arma::vec y=beta.row(j).t(),s=se2.row(j).t();
            y=y.elem(active);s=s.elem(active);
            arma::mat T=VV+arma::diagmat(s),precision=arma::inv_sympd(T);
            arma::vec t=precision*y;
            double score=.5*(arma::accu(arma::square(y)/s)-arma::as_scalar(y.t()*precision*y)+
                          arma::accu(arma::log(s))-arma::log_det_sympd(T));
            const double weight=stats[k].weight(score,prior(j,k));
            stats[k].H+=weight*arma::symmatu(t*t.t()-precision);
          }
      }
    }
    double shift=-arma::datum::inf,mass=0;
    for(const auto& item:stats)shift=std::max(shift,item.shift);
    for(const auto& item:stats)mass+=std::exp(item.shift-shift)*item.mass;
    if(!std::isfinite(shift) || !std::isfinite(mass) || mass<=0)
      Rcpp::stop("Nonfinite covariance-EM objective");
    objective=shift+std::log(mass);
    if(iteration>0 && std::abs(objective-previous)<tol) {status=3;break;}
    if(iteration==maxiter-1)break;
    previous=objective;
    arma::mat next=V;
    for(int k=0;k<nc;++k) {
      const arma::uvec& active=cached[k].active;
      arma::mat H=stats[k].H*(std::exp(stats[k].shift-shift)/mass);
      arma::mat VA=V.cols(active);
      next+=VA*H*VA.t();
    }
    arma::vec eig;arma::mat U;arma::eig_sym(eig,U,arma::symmatu(next));
    eig.transform([](double x){return std::max(x,1e-12);});
    V=arma::symmatu(U*arma::diagmat(eig)*U.t());
  }
  return Rcpp::List::create(Rcpp::Named("V")=V,Rcpp::Named("objective")=-objective,
                           Rcpp::Named("iterations")=iteration+1,Rcpp::Named("status")=status);
}
