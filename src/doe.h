#ifndef METAOPT_DOE_H
#define METAOPT_DOE_H

#include "random.h"
#include "sample.h"

// using Random = effolkronium::random_static;

namespace MetaOpt {

enum DoeEnum {
  ZERO,  // 0 zero to lower bound
  MC,    // 1 Monte-Carlo
  LHS,   // 2 Latin Hypercubic Sampling
};

template <typename Real, int k>
class Doe
{
 public:
  static Samples<Real> gen(const int                n_sample,
                           const int                n_x,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0);
};

template <typename Real>
class Doe<Real, ZERO>
{
 public:
  static Samples<Real> gen(const int                n_sample,
                           const int                n_x,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_sample);
    if (upper.empty()) {
      for (auto& s : samples) {
        s = Sample<Real>(n_x, n_obj, n_con);
      }
    } else {
      for (auto& s : samples) {
        s = Sample<Real>(n_x, n_obj, n_con);
        for (int i = 0; i != n_x; ++i)
          s.x()[i] = lower[i];
      }
    }
    return samples;
  }
};

template <typename Real>
class Doe<Real, MC>
{
 public:
  static Samples<Real> gen(const int                n_sample,
                           const int                n_x,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_sample);
    for (auto& s : samples) {
      s = Sample<Real>(n_x, n_obj, n_con);
    }
    if (upper.empty()) {
      for (int d = 0; d != n_x; ++d) {
        for (int i = 0; i != n_sample; ++i) {
          samples[i].x()[d] = Random.gen<Real>(0, 1);
        }
      }

    } else {
      for (int d = 0; d != n_x; ++d) {
        for (int i = 0; i != n_sample; ++i) {
          samples[i].x()[d] = Random.gen<Real>(lower[d], upper[d]);
        }
      }
    }
    return samples;
  }
};

template <typename Real>
class Doe<Real, LHS>
{
 public:
  static Samples<Real> gen(const int                n_sample,
                           const int                n_x,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_sample);
    for (auto& s : samples) {
      s = Sample<Real>(n_x, n_obj, n_con);
    }
    std::vector<int> list(n_sample);
    int              i = 0;
    for (auto& t : list)
      t = i++;
    if (upper.empty()) {
      for (int d = 0; d != n_x; ++d) {
        Random.shuffle(list);
        for (int i = 0; i != n_sample; ++i) {
          samples[i].x()[d] = (Random.gen<Real>(0, 1) + list[i]) / n_sample;
        }
      }

    } else {
      for (int d = 0; d != n_x; ++d) {
        Random.shuffle(list);
        Real range = upper[d] - lower[d];
        for (int i = 0; i != n_sample; ++i) {
          samples[i].x()[d] =
              lower[d] + range * (Random.gen<Real>(0, 1) + list[i]) / n_sample;
        }
      }
    }
    return samples;
  }
};

}  // namespace MetaOpt
#endif  // METAOPT_DOE_H