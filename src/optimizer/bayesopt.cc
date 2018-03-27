#include "bayesopt.h"
#include "common.h"

namespace MetaOpt {

template <typename Real>
BayesOpt<Real>::BayesOpt(const int               n_dim,
         const int               n_obj,
         const int               n_con,
         const Func<Real>        func,
         const std::vector<Real> lower,
         const std::vector<Real> upper,
         const int acc)
  :Optimizer<Real>(n_dim, n_obj, n_con, func, lower, upper){
  gp_.reset(new Gp<Real>());
  ga_.reset(new Ga<Real>(n_dim, n_obj, n_con, nullptr, lower, upper));
};

template <typename Real>
void BayesOpt<Real>::opt(){
  samples_ = Doe<Real, LHS>::gen(10, this->n_dim_, this->lower_, this->upper_, this->n_obj_, this->n_con_);
  for(auto& s: samples_)
    this->evaluate(s);
  std::vector<Real> t(0);
  for(int i=samples_.size();i!=100;++i){
    gp_->fit(samples_);
    ga_->set_func(std::bind(&Gp<Real>::evaluate, gp_.get(), std::placeholders::_1, t));
    ga_->opt();
    this->evaluate(ga_->best_);
    samples_.push_back(ga_->best_);
    LOG(INFO) << samples_.back();
  }
};


template class BayesOpt<double>;
}  // namespace MetaOpt