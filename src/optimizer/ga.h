#ifndef METAOPT_OPTIMIZER_GP_H
#define METAOPT_OPTIMIZER_GP_H
#include "exd/random.hpp"
#include "optimizer.h"

namespace MetaOpt {
template <typename Real>
class Ga : public Optimizer<Real>
{
 public:
  Ga(const int                n_x = 0,
     const int                n_obj = 1,
     const int                n_con = 0,
     const Func<Real>         func = nullptr,
     const std::vector<Real>& lower = std::vector<Real>(),
     const std::vector<Real>& upper = std::vector<Real>(),
     const int                n_population = 0,
     const int                n_generation = 0,
     const Real               p_crossover = 0.9,
     const Real               p_mutation = 0.05);
  ~Ga();
  void          opt();
  void          tour_select(int k);
  void          simubinary_crossover(int k);
  void          liner_crossover(int k);
  void          PBMutation(int k);
  Sample<Real>& select(Sample<Real>& s1, Sample<Real>& s2);
  void          crossover(Sample<Real>& s1, Sample<Real>& s2);
  void          mutation(Sample<Real>& s);
  void          statistics();

  int              n_population_;
  int              n_generation_;
  int              i_generation_;
  Samples<Real>    population_;
  Samples<Real>    newpop_;
  std::vector<int> tourlist_;
  Real             p_crossover_;
  Real             p_mutation_;
  Real             i_mutation_;
  Real             i_crossover_;
  Real             nm = 1, nc = 1;
};

}  // namespace MetaOpt
#endif  // METAOPT_OPTIMIZER_GP_H
