#ifndef METAOPT_OPTIMIZER_H
#define METAOPT_OPTIMIZER_H

#include "sample.h"

namespace MetaOpt {

template <typename Real>
using Func = std::function<void(Real* x, Real* y, Real* c)>;
template <typename Real>
using Func2 = std::function<void(Sample<Real>& sample)>;

template <typename Real>
class Optimizer
{
 public:
  Optimizer(const int                n_x = 0,
            const int                n_obj = 1,
            const int                n_con = 0,
            const Func<Real>         func = nullptr,
            const std::vector<Real>& lower = std::vector<Real>(),
            const std::vector<Real>& upper = std::vector<Real>())
      : n_x_(n_x),
        n_obj_(n_obj),
        n_con_(n_con),
        func_(func),
        lower_(lower),
        upper_(upper),
        n_evaluation_(0) {
    best_ = std::move(sample());
    if (!upper_.empty()) {
      for (int i = 0; i != n_x_; ++i) {
        if (upper_[i] < lower_[i]) {
          std::swap(upper_[i], lower_[i]);
        }
      }
    }
  };
  ~Optimizer(){};
  Sample<Real> sample() {
    Sample<Real> x(n_x_, n_obj_, n_con_);
    if (!upper_.empty()) {
      for (int i = 0; i != n_x_; ++i)
        x.x()[i] = exd::random::random<Real>(lower_[i], upper_[i]);
    }
    return x;
  };
  Sample<Real> get_best() { return best_; };
  void         set_func(Func<Real> func) { func_ = func; };
  void         set_func(Func2<Real> func) {
    func2_ = func;
    func_ = nullptr;
  };
  void evaluate(Sample<Real>& sample, int n_thread = 1) {
    if (func_)
      sample.evaluate(func_);
    else if (func2_)
      func2_(sample);
    ++n_evaluation_;
  };
  void              opt();
  int               n_x_;
  int               n_obj_;
  int               n_con_;
  Func<Real>        func_;
  Func2<Real>       func2_;
  std::vector<Real> upper_;
  std::vector<Real> lower_;
  Sample<Real>      best_;
  int               n_evaluation_;
};
}  // namespace MetaOpt

#endif  // METAOPT_OPTIMIZER_H
