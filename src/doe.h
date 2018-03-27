#ifndef METAOPT_DOE_H
#define METAOPT_DOE_H

#include "random.hpp"
#include "sample.h"

using Random = effolkronium::random_static;

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
  static Samples<Real> gen(const int                n_x,
                           const int                n_dim,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0);
};

template <typename Real>
class Doe<Real, ZERO>
{
 public:
  static Samples<Real> gen(const int                n_x,
                           const int                n_dim,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_x);
    if (upper.empty()) {
      for (auto& s : samples) {
        s = Sample<Real>(n_dim, n_obj, n_con);
      }
    } else {
      for (auto& s : samples) {
        s = Sample<Real>(
            std::vector<Real>(n_dim, 0),
            std::vector<Real>(n_obj, std::numeric_limits<Real>::infinity()),
            std::vector<Real>(n_con, -1));
      }
    }
    return samples;
  }
};

template <typename Real>
class Doe<Real, MC>
{
 public:
  static Samples<Real> gen(const int                n_x,
                           const int                n_dim,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_x);
    for (auto& s : samples) {
      s = Sample<Real>(n_dim, n_obj, n_con);
    }
    if (upper.empty()) {
      for (int d = 0; d != n_dim; ++d) {
        for (int i = 0; i != n_x; ++i) {
          samples[i].x()[d] = Random::get<Real>(0, 1);
        }
      }

    } else {
      for (int d = 0; d != n_dim; ++d) {
        for (int i = 0; i != n_x; ++i) {
          samples[i].x()[d] = Random::get<Real>(lower[d], upper[d]);
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
  static Samples<Real> gen(const int                n_x,
                           const int                n_dim,
                           const std::vector<Real>& lower = std::vector<Real>(),
                           const std::vector<Real>& upper = std::vector<Real>(),
                           const int                n_obj = 1,
                           const int                n_con = 0) {
    Samples<Real> samples(n_x);
    for (auto& s : samples) {
      s = Sample<Real>(n_dim, n_obj, n_con);
    }
    std::vector<int> list(n_x);
    int              i = 0;
    for (auto& t : list)
      t = i++;
    if (upper.empty()) {
      for (int d = 0; d != n_dim; ++d) {
        Random::shuffle(list);
        for (int i = 0; i != n_x; ++i) {
          samples[i].x()[d] = (Random::get<Real>(0, 1) + list[i]) / n_x;
        }
      }

    } else {
      for (int d = 0; d != n_dim; ++d) {
        Random::shuffle(list);
        Real range = upper[d] - lower[d];
        for (int i = 0; i != n_x; ++i) {
          samples[i].x()[d] =
              lower[d] + range * (Random::get<Real>(0, 1) + list[i]) / n_x;
        }
      }
    }
    return samples;
  }
};

}  // namespace MetaOpt
#endif  // METAOPT_DOE_H