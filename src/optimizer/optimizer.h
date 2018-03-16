#ifndef METAOPT_OPTIMIZER_H
#define METAOPT_OPTIMIZER_H

#include "common.h"

namespace MetaOpt {

template <typename Real>
using Func = std::function<void(Real* x, Real* y, Real* c)>;

template <typename Real>
using Sample =
    std::tuple<std::vector<Real>, std::vector<Real>, std::vector<Real>>;

// template <typename Real>
// class Sample : public SampleBase<Real>{
//   // std::vector<Real> x(){
//   //   return std::get<0>(*this);
//   // };
//   // std::vector<Real> obj(){
//   // };
//   // std::vector<Real> con(){
//   // };
// };

template <typename Real>
class Optimizer
{
 public:
  Optimizer(const int               n_dim = 0,
            const int               n_obj = 0,
            const int               n_con = 0,
            const Func<Real>        func = nullptr,
            const std::vector<Real> upper = std::vector<Real>(),
            const std::vector<Real> lower = std::vector<Real>())
      : n_dim_(n_dim),
        n_obj_(n_obj),
        n_con_(n_con),
        func_(func),
        upper_(upper),
        lower_(lower) {
    sample().swap(best_);
    for (int i = 0; i != n_dim_; ++i){
      if(upper_[i]< lower_[i]){
        std::swap(upper_[i], lower_[i]);
      }
    }
    LOG(INFO) << "Optimizer constructed, "
              << "n_dim = " << n_dim_;
    LOG(INFO) << best_;
  };
  ~Optimizer() { LOG(INFO) << "Optimizer destructed"; };
  Sample<Real> sample() {
    std::vector<Real> x(n_dim_, 0);
    for (int i = 0; i != n_dim_; ++i)
      x[i] = Random::get<Real>(lower_[i], upper_[i]);
    return Sample<Real>(
        std::move(x),
        std::vector<Real>(n_obj_, std::numeric_limits<Real>::infinity()),
        std::vector<Real>(n_con_, 0));
  };
  Sample<Real> get_best() { return best_; };
  void set_func(Func<Real> func) { func_ = func; };
  void evaluate(Sample<Real>& sample, int n_thread = 1) {
    func_(std::get<0>(sample).data(), std::get<1>(sample).data(),
          std::get<2>(sample).data());
  };
  void              opt();
  int               n_dim_;
  int               n_obj_;
  int               n_con_;
  Func<Real>        func_;
  std::vector<Real> upper_;
  std::vector<Real> lower_;
  Sample<Real>      best_;
  int               n_evaluation_;
};
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
  s.put('[');
  char comma[3] = {'\0', ' ', '\0'};
  for (const auto& e : v) {
    s << comma << e;
    comma[0] = ',';
  }
  return s << ']';
}

template <typename Real>
std::ostream& operator<<(std::ostream& s, const MetaOpt::Sample<Real>& sample) {
  s << "x = " << std::get<0>(sample) << ", obj = " << std::get<1>(sample)
    << ", con = " << std::get<2>(sample);
  return s;
}

#endif  // METAOPT_OPTIMIZER_H
