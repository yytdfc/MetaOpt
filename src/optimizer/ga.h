#ifndef METAOPT_GA_H
#define METAOPT_GA_H
#include <string>
#include <functional>

template <typename Real>
class GA
{
 public:
  ~GA();
  void GAinit(int   nv,
              int   size,
              int   nGen,
              int   nCons,
              Real  pCr,
              Real  pMu,
              Real* up,
              Real* low);
  void  setFx(std::function<Real(Real*)>);
  void  evolve();
  Real* best;
  void  test();

  void tourSelect(int k);
  void SBXover(int k);
  void linerXover(int k);
  void PBMutation(int k);
  void rdMutation(int k);
  void initialpop();
  void statistics();
  void tourney();

 private:
  std::function<Real(Real*)> fx;
  int*                       tourlist;

  int    nPop;  /* population size */
  int    nGens; /* max. number of generations */
  int    nVar;  /* no. of problem variables */
  int    nCons;
  int    generation;
  Real*  consP;
  Real   pCrossover; /* probability of SBXover */
  Real   pMutation;  /* probability of PBMutation */
  Real** pop;
  Real** newpop;
  Real*  upper; /* GT\'s variables upper bound */
  Real*  lower; /* GT\'s variables lower bound */

  Real nMutation;
  Real nCrossover;
  Real nm, nc;

};
#endif // METAOPT_GA_H