#ifndef METAOPT_POWELL_H
#define METAOPT_POWELL_H
#include <string>
#include <vector>
#include <functional>

template <typename Real>
class Powell
{
 public:
  Powell();
  Powell(int _nv);
  Powell(int _nv, std::function<Real(Real* x)> f);
  ~Powell();
  void SetFx(std::function<Real(Real* x)>);
  void SetBracktRange(Real range);
  void InitDirec(Real** direc);
  void InitReverseDirec(Real** direc);
  void InitRandomDirec(Real** direc);
  void LineSearch(Real p[], Real xi[]);
  Real f1dim2(Real alpha, Real* x, Real* p, Real* temp);
  Real Brent(Real  xa,
             Real  xb,
             Real  xc,
             Real  tol,
             Real& xmin,
             std::function<Real(Real x)> func1d);
  void erase(Real pbar[], Real prr[], Real pr[]);
  int evolve(Real*  x0,
             Real** direc,
             int    maxiter,
             Real   ftol,
             int    dispIter = 0,
             Real   terminalLine = -10000);
  int Optimize(Real p[],
               Real ftol,
               int  maxiter = 1000,
               int  dispIter = 0,
               Real terminalLine = -10000);
  int amoeba(Real* x, Real ftol, int dispIter);
  void Mnbrak(Real& ax,
              Real& bx,
              Real& cx,
              Real& fa,
              Real& fb,
              Real& fc,
              std::function<Real(Real x)> func1d);
  void Bracket(Real& ax,
               Real& bx,
               Real& cx,
               Real& fa,
               Real& fb,
               Real& fc,
               std::function<Real(Real x)> func1d);

 private:
  std::function<Real(Real* x)> func;
  int                          nv;
  Real                         brackeRange;

  inline Real sgn(Real x) { return x < 0 ? -1 : (x == 0 ? 0 : 1); }
};
#endif // METAOPT_POWELL_H