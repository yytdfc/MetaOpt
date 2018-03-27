#ifndef METAOPT_OPTIMIZER_BAYESOPT_H
#define METAOPT_OPTIMIZER_BAYESOPT_H
#include "doe.h"
#include "model/gp.h"
#include "optimizer.h"
#include "optimizer/ga.h"

namespace MetaOpt {
template <typename Real>
class BayesOpt : public Optimizer<Real>
{
 public:
  BayesOpt(const int               n_dim = 0,
           const int               n_obj = 1,
           const int               n_con = 0,
           const Func<Real>        func = nullptr,
           const std::vector<Real> lower = std::vector<Real>(),
           const std::vector<Real> upper = std::vector<Real>(),
           const int               acc = 12);
  ~BayesOpt(){};
  void                      opt();
  std::unique_ptr<Gp<Real>> gp_;
  std::unique_ptr<Ga<Real>> ga_;
  Samples<Real>             samples_;
};

}  // namespace MetaOpt

#endif  // METAOPT_OPTIMIZER_BAYESOPT_H