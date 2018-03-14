#pragma once
#include <string>
#include <random>
#include <vector>
#include <functional>
typedef std::function< double (double *x) > FUNC;
typedef std::function< double (double x) > FUNC1D;
class Powell {
public:
	Powell();
	Powell(int _nv);
	Powell(int _nv, FUNC f);
	~Powell();
	void SetFx(FUNC);
	void SetBracktRange(double range);
	void InitDirec(double **direc);
	void InitReverseDirec(double **direc);
	void InitRandomDirec(double **direc);
	void LineSearch(double p[], double xi[]);
	double f1dim2(double alpha, double *x, double *p, double *temp);
	double Brent(double xa,double xb,double xc,double tol,double& xmin,FUNC1D func1d);
	void erase(double pbar[], double prr[], double pr[]);
	int evolve(double *x0,double **direc,int maxiter, double ftol, int dispIter=0, double terminalLine=-10000);
	int Optimize(double p[], double ftol, int maxiter=1000, int dispIter=0, double terminalLine=-10000);
	int amoeba(double *x, double ftol, int dispIter);
	void Mnbrak(double& ax, double& bx, double& cx, double& fa,
			double& fb, double& fc,FUNC1D func1d);
	void Bracket(double& ax, double& bx, double& cx, double& fa,
			double& fb, double& fc,FUNC1D func1d);
private:
	FUNC func;
	int nv;
	double brackeRange;
//	std::random_device rd;  //�����������
	std::mt19937_64 gen;
	///*return a random double in [0,1]*/
	inline double randomDouble(){
		return (double)gen() / gen.max();
	}
	inline double sgn(double x){
		return x<0 ? -1:(x==0 ? 0:1);
	}
};