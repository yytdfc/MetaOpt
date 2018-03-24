#ifndef METAOPT_MODEL_KERNEL_H
#define METAOPT_MODEL_KERNEL_H
#include <vector>
namespace MetaOpt {

enum KernelEnum {
  // GAUSS:
  GAUSS,     // 0 GAUSS function
  CSPL,      // 1 Cubic Spline function
  SPLINE,    // 2 Spline function
  MATERN32,  // 3 Spline function
  MATERN52,  // 4 Spline function
  // RBF:
  IMQ,  // inverse multi-quadratic RBFs
  MQ,   // multi-quadratic RBFs
  TPS,  // thin plate spline RBFs
  POW,  // power RBFs
  // Legacy:
  GEXP = 99,  // 99 Gaussian Exponential function
};

template <typename Real, int k>
class Kernel
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale);
};

template <typename Real>
Real R2(const std::vector<Real>& x1, const std::vector<Real>& x2) {
  Real r2 = 0.0;
  for (int i = 0; i != x1.size(); i++) {
    r2 += (x1[i] - x2[i]) * (x1[i] - x2[i]);
  }
  return r2;
}

template <typename Real>
class Kernel<Real, GAUSS>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real r2 = R2<Real>(x1, x2);
    return exp(-r2 / (2 * length_scale * length_scale));
  };
};

template <typename Real>
class Kernel<Real, MATERN32>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real r = sqrt(3 * R2<Real>(x1, x2)) / length_scale;
    return (1 + r) * exp(-r);
  };
};

template <typename Real>
class Kernel<Real, MATERN52>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real r = sqrt(5 * R2<Real>(x1, x2)) / length_scale;
    return (1 + r + r * r / 3.) * exp(-r);
  };
};

template <typename Real>
class Kernel<Real, CSPL>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real temp = 1.0;
    for (int i = 0; i != x1.size(); ++i) {
      Real r = fabs(x1[i] - x2[i]);
      Real kesi = r / length_scale;
      if (kesi >= 0.0 && kesi <= 0.2)
        temp *= 1.0 - 15.0 * kesi * kesi * (1.0 - 2.0 * kesi);
      else if (kesi > 0.2 && kesi < 1)
        temp *= 1.25 * (1.0 - kesi) * (1.0 - kesi) * (1.0 - kesi);
      else
        temp *= 0.0;
    }
    return (temp);
  };
};

template <typename Real>
class Kernel<Real, SPLINE>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real temp = 1.0;
    for (int i = 0; i != x1.size(); ++i) {
      Real r = fabs(x1[i] - x2[i]);
      Real kesi = r / length_scale;
      if (kesi >= 0.0 && kesi <= 0.4)
        temp *= 1.0 + kesi * kesi * (-15.0 + kesi * (35.0 + kesi * (-24.375)));
      else if (kesi > 0.4 && kesi < 1)
        temp *=
            (5.0 / 3.0) +
            kesi *
                ((-20.0 / 3.0) +
                 kesi * (10.0 + kesi * ((-20.0 / 3.0) + kesi * (5.0 / 3.0))));
      else
        temp *= 0.0;
    }
    return (temp);
  };
};

template <typename Real>
class Kernel<Real, IMQ>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    return 1.0 / sqrt(1.0 + R2<Real>(x1, x2));
  };
};

template <typename Real>
class Kernel<Real, MQ>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    return sqrt(1.0 + R2<Real>(x1, x2));
  };
};

template <typename Real>
class Kernel<Real, TPS>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real r2 = R2<Real>(x1, x2);
    r2 *= 10.0;
    if (r2 <= 1E-6)
      return 0;
    else
      return r2 * log(r2);
  };
};

template <typename Real>
class Kernel<Real, POW>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    return pow(R2<Real>(x1, x2), 0.75);
  };
};

template <typename Real>
class Kernel<Real, GEXP>
{
 public:
  static Real correlation(const std::vector<Real>& x1,
                          const std::vector<Real>& x2,
                          const Real               length_scale) {
    Real r2 = R2<Real>(x1, x2);
    return exp(-pow(r2, 0.5) / (2 * length_scale * length_scale));
  };
};
}  // namespace MetaOpt
#endif  // METAOPT_MODEL_KERNEL_H
