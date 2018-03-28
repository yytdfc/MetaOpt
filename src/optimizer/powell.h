#ifndef METAOPT_OPTIMIZER_POWELL_H
#define METAOPT_OPTIMIZER_POWELL_H
#include "optimizer.h"

namespace MetaOpt {
template <typename Real>
class Powell : public Optimizer<Real>
{
 public:
  Powell(const int n_x = 0, const Func<Real> func = nullptr)
      : Optimizer<Real>(n_x, 1, 0, func){};
  ~Powell(){};
  void SetFx(std::function<Real(Real* x)>);
  void SetBracktRange(Real range);
  void InitDirec(Real** direc);
  void InitReverseDirec(Real** direc);
  void InitRandomDirec(Real** direc);
  void LineSearch(Sample<Real>&, Real xi[]);
  Real f1dim2(Real alpha, Sample<Real>& x, Real* p, Sample<Real>& temp);
  Real Brent(Real                        xa,
             Real                        xb,
             Real                        xc,
             Real                        tol,
             Real&                       xmin,
             std::function<Real(Real x)> func1d);
  void erase(Real pbar[], Real prr[], Real pr[]);
  int  evolve(Sample<Real>& x0,
              Real**        direc,
              int           maxiter,
              Real          ftol,
              int           dispIter = 0,
              Real          terminalLine = -10000);
  void opt(Sample<Real>& p,
           Real          ftol = 1e-8,
           int           maxiter = 1000,
           int           dispIter = 0,
           Real          terminalLine = -10000);
  int  amoeba(Real* x, Real ftol, int dispIter);
  void Mnbrak(Real&                       ax,
              Real&                       bx,
              Real&                       cx,
              Real&                       fa,
              Real&                       fb,
              Real&                       fc,
              std::function<Real(Real x)> func1d);
  void Bracket(Real&                       ax,
               Real&                       bx,
               Real&                       cx,
               Real&                       fa,
               Real&                       fb,
               Real&                       fc,
               std::function<Real(Real x)> func1d);

 private:
  Real brackeRange;

  inline Real sgn(Real x) { return x < 0 ? -1 : (x == 0 ? 0 : 1); }
};
}  // namespace MetaOpt
#endif  // METAOPT_OPTIMIZER_POWELL_H