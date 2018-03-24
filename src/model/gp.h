#ifndef METAOPT_MODEL_GP_H
#define METAOPT_MODEL_GP_H

#include <eigen3/Eigen/Core>
#include "model.h"

namespace MetaOpt {

template <typename Real>
using Matrixx = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Real>
using KernelFun = std::function<Real(const std::vector<Real>&,
                                     const std::vector<Real>&,
                                     const Real length_scale)>;

template <typename Real>
class Gp : public Model<Real>
{
 public:
  Gp(){};
  Gp(const Samples<Real>& samples) : Model<Real>(samples){};
  ~Gp(){};
  void init(const Samples<Real>& samples);
  void evaluate(Sample<Real>&      sample,
                std::vector<Real>& mse = std::vector<Real>(0));
  void evaluate(Samples<Real>& samples);
  void fit(const Samples<Real>& samples = Samples<Real>(0));
  void obj_mle(Real* hyper_l, Real* obj, Real* con = nullptr);

  // paras
  int             i_kernel_ = 3;
  int             i_order_ = 1;
  KernelFun<Real> kernel_fun_;
  Real            mu_ = 0;
  Real            sigma_ = 0;
  int             f_dim_ = 1;
  // model data
  Matrixx<Real>                                           k_;
  Real                                                    length_scale_;
  std::vector<Real>                                       pk_;
  Matrixx<Real>                                           alpha_;
  std::shared_ptr<Eigen::LDLT<Eigen::Ref<Matrixx<Real>>>> ldlt_;
};
}  // namespace MetaOpt

#endif  // METAOPT_MODEL_GP_H
