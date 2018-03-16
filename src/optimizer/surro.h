#ifndef METAOPT_SURRO_H
#define METAOPT_SURRO_H
#include "ga.h"
#include "doe/doe.h"
#include "model/gkrig.h"
#include <string>
#include <vector>
#include <ctime>

template <typename Real>
class Surro
{
 public:
  void initialize(int   nv,
                  int   nS,
                  int   init,
                  int   nCon,
                  Real* up,
                  Real* low,
                  std::function<Real(Real* x)>);
  ~Surro();
  int readInput(std::string infile);
  void initSample();
  void getBest(const int&);
  void add(int);
  void opt();
  void dealAdd();
  void dealRestart();
  void allocate();
  void setFx(Real(Real* x));
  void setGKrigFx(Real(Real* x));
  void backupResult();
  int  n_dim_;
  int  nCons;
  int  nInit;
  int  nSample;
  int  nNow;
  int  vecLengh;

  int                   nProb;
  bool                  isInit;
  int                   isRestart;
  string                restartFrom;
  int                   addN;
  int*                  addMethod;
  int                   nDoE;
  Real**                sample;
  Real*                 response;
  Real**                cons;
  Real*                 upper;
  Real*                 lower;
  Real*                 convergence;
  Real*                 best;
  Real                  bestfit, *bestcons;
  function<Real(Real*)> fx;

  // Real predictorMP(Real *, Real *);
  // Real predictorEI(Real *, Real *);
  // Real predictorME(Real *, Real *);

 private:
  GKrig<Real> krig;
  GA<Real>    ga;
  clock_t     timer;
  // GA parameters
  int  ga_Pops;
  int  ga_Gens;
  Real ga_Pcr;
  Real ga_Pmu;

  // Kriging parameters
  int krig_corr;
  int krig_const_theta;
  int krig_porder;
  int krig_norm;
  int krig_dcmp;
  int krig_ParaOpt;
  int krig_regular;
  int krig_out_points;
  int krig_EIcons;
};

#endif  // METAOPT_SURRO_H