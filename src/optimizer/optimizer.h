#ifndef METAOPT_OPTIMIZER_H
#define METAOPT_OPTIMIZER_H

#include "common.h"
#include "sample.h"

namespace MetaOpt {

template <typename Real>
using Func = std::function<void(Real *x, Real *y, Real *c)>;

template <typename Real> class Optimizer {
public:
  Optimizer(const int n_dim = 0, const int n_obj = 0, const int n_con = 0,
            const Func<Real> func = nullptr,
            const std::vector<Real> upper = std::vector<Real>(),
            const std::vector<Real> lower = std::vector<Real>())
      : n_dim_(n_dim), n_obj_(n_obj), n_con_(n_con), func_(func), upper_(upper),
        lower_(lower), n_evaluation_(0) {
    best_ = std::move(sample());
    for (int i = 0; i != n_dim_; ++i) {
      if (upper_[i] < lower_[i]) {
        std::swap(upper_[i], lower_[i]);
      }
    }
  };
  ~Optimizer(){};
  Sample<Real> sample() {
    std::vector<Real> x(n_dim_, 0);
    for (int i = 0; i != n_dim_; ++i)
      x[i] = Random::get<Real>(lower_[i], upper_[i]);
    return Sample<Real>(
        std::move(x),
        std::vector<Real>(n_obj_, std::numeric_limits<Real>::infinity()),
        std::vector<Real>(n_con_, -1));
  };
  Sample<Real> get_best() { return best_; };
  void set_func(Func<Real> func) { func_ = func; };
  void evaluate(Sample<Real> &sample, int n_thread = 1) {
    func_(sample.x().data(), sample.obj().data(), sample.con().data());
    ++n_evaluation_;
  };
  void opt();
  int n_dim_;
  int n_obj_;
  int n_con_;
  Func<Real> func_;
  std::vector<Real> upper_;
  std::vector<Real> lower_;
  Sample<Real> best_;
  int n_evaluation_;
};
} // namespace MetaOpt

#endif // METAOPT_OPTIMIZER_H
