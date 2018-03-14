#include "ga.h"
#include "doe/doe.h"
#include "model/gkrig.h"
#include <string>
#include <random>
#include <vector>
#include <ctime>
typedef function< double(double *) > FUNC;
class Surro {
public:
	void initialize(int nv, int nS, int init, int nCon, double *up, double *low, double(*f)(double *));
	~Surro();
	int readInput(std::string infile);
	void initSample();
	void getBest(const int &);
	void add(int);
	void opt();
	void dealAdd();
	void dealRestart();
	void allocate();
	void setFx(double(*f)(double *));
	void setGKrigFx(double(*f)(double *));
	void backupResult();
	int nVar;
	int nCons;
	int nInit;
	int nSample;
	int nNow;
	int vecLengh;

	int nProb;
	bool isInit;
	int isRestart;
	string restartFrom;
	int addN;
	int *addMethod;
	int nDoE;
	double **sample;
	double *response;
	double **cons;
	double *upper;
	double *lower;
	double *convergence;
	double *best;
	double bestfit, *bestcons;
	FUNC fx;

	//double predictorMP(double *, double *);
	//double predictorEI(double *, double *);
	//double predictorME(double *, double *);

private:
	GKrig krig;
	GA ga;
	clock_t timer;
	// GA parameters
	int ga_Pops;
	int ga_Gens;
	double ga_Pcr;
	double ga_Pmu;

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

	std::random_device rd;
	std::mt19937_64 gen;
	/*return a random double in [0,1]*/
	inline double randomDouble(){
		return (double)gen() / gen.max();
	}
	/*return a random double in [a,b], numbers can be gened*/
	inline double randomDouble(double a, double b){
		return (b - a) * gen() / gen.max() + a;
	}
	/*return a random int in [a,b], (b-a+1) numbers can be gened*/
	inline int randomInt(int a, int b){
		return gen() % (b - a + 1) + a;
	}
	inline double max(double a, double b){
		return (a > b) ? a : b;
	}
	inline double min(double a, double b){
		return (a < b) ? a : b;
	}
	inline void swap(double &a, double &b){
		double temp = a;
		a = b;
		b = temp;
	}
};

