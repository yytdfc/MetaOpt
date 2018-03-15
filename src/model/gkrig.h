#ifndef METAOPT_GKRIG_H
#define METAOPT_GKRIG_H
#include "kriging.h"
#include <string>

template <typename Real>
class GKrig : public Kriging<Real>
{
 public:
  void GKinitialize(int    ncorr,
                    int    nconst_theta,
                    int    nporder,
                    int    nnorm,
                    int    ndcmp,
                    int    nParaOpt,
                    int    nregular,
                    int    ndim,
                    int    points,
                    int    nout_points,
                    int    nny,
                    Real** xx,
                    Real*  up,
                    Real*  low);
  void GKsetPredict(int k);
  void GKtraining();
  void GKoutputRSM();
  void GKprediction();
  void setFx(Real f(Real*));
  Real fxIndex(Real* x, int);
  Real (*fx)(Real*);

  Real GKpredictorMP(Real* x);
  Real GKpredictorEI(Real* x);

 private:
  int  yIndex;
  Real GKpredictor(Real*);
  Real GKMSE(Real*);
  int  t; /* Number of levels of data                        */
};
#endif  // METAOPT_GKRIG_H