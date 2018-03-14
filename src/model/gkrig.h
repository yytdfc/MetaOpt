#include "kriging.h"
#include <string>

class GKrig:public Kriging{
public:
	void GKinitialize(int ncorr, int nconst_theta, int nporder,
		int nnorm, int ndcmp, int nParaOpt, int nregular, int ndim,
		int points, int nout_points, int nny, 
		double **xx, double *up, double *low);
	void GKsetPredict(int k);
	void GKtraining();
	void GKoutputRSM();
	void GKprediction();
	void setFx(double f(double*));
	double fxIndex(double *x, int);
	double(*fx)(double *);

	double GKpredictorMP(double *x);
	double GKpredictorEI(double *x);

private:

	int yIndex;
	double GKpredictor(double *);
	double GKMSE(double *);
	int t; /* Number of levels of data                        */


};