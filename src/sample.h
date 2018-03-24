#ifndef METAOPT_SAMPLE_H
#define METAOPT_SAMPLE_H

#include <iostream>
#include <vector>
#include <functional>

namespace MetaOpt {

template <typename Real>
using Func = std::function<void(Real* x, Real* y, Real* c)>;

template <typename Real>
class Sample
{
 public:
  Sample(){};
  Sample(int n_x, int n_obj, int n_con): x_(std::vector<Real>(n_x, 0)),
      obj_(std::vector<Real>(n_obj, std::numeric_limits<Real>::infinity())),
      con_(std::vector<Real>(n_con, -1)){};
  Sample(std::vector<Real> x, std::vector<Real> obj, std::vector<Real> con)
      : x_(x), obj_(obj), con_(con){};
  Sample(const Sample& s) : x_(s.x_), obj_(s.obj_), con_(s.con_){};
  Sample(Sample&& s) : x_(move(s.x_)), obj_(move(s.obj_)), con_(move(s.con_)){};
  Sample& operator=(const Sample& s) {
    x_ = s.x_;
    obj_ = s.obj_;
    con_ = s.con_;
    return *this;
  };
  Sample& operator=(Sample&& s) {
    x_ = move(s.x_);
    obj_ = move(s.obj_);
    con_ = move(s.con_);
    return *this;
  };
  void swap(Sample& s) {
    auto t(std::move(*this));
    *this = std::move(s);
    s = std::move(t);
  };
  friend std::ostream& operator<<(std::ostream& os, const Sample<Real>& s) {
    os << "x = " << s.x_ << ", obj = " << s.obj_ << ", con = " << s.con_;
    return os;
  };
  void evaluate(const Func<Real>& func) {
    func(x_.data(), obj_.data(), con_.data());
  };
  std::vector<Real>& x() { return x_; };
  std::vector<Real>& obj() { return obj_; };
  std::vector<Real>& con() { return con_; };

 private:
  std::vector<Real> x_, obj_, con_;
};

template <typename Real>
using Samples = std::vector<Sample<Real>>;
}  // namespace MetaOpt

#endif  // METAOPT_SAMPLE_H