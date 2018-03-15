#include <string>
#include <random>
#include <functional>
typedef std::function<double(double*)> FUNC;
class GA
{
 public:
  ~GA();
  void GAinit(int     nv,
              int     size,
              int     nGen,
              int     nCons,
              double  pCr,
              double  pMu,
              double* up,
              double* low);
  void    setFx(FUNC);
  void    evolve();
  double* best;
  void    test();

  void tourSelect(int k);
  void SBXover(int k);
  void linerXover(int k);
  void PBMutation(int k);
  void rdMutation(int k);
  void initialpop();
  void statistics();
  void tourney();

 private:
  FUNC fx;
  int* tourlist;

  int      nPop;  /* population size */
  int      nGens; /* max. number of generations */
  int      nVar;  /* no. of problem variables */
  int      nCons;
  int      generation;
  double*  consP;
  double   pCrossover; /* probability of SBXover */
  double   pMutation;  /* probability of PBMutation */
  double** pop;
  double** newpop;
  double*  upper; /* GT\'s variables upper bound */
  double*  lower; /* GT\'s variables lower bound */

  double nMutation;
  double nCrossover;
  double nm, nc;

  std::random_device rd;
  std::mt19937_64    gen;
  /*return a random double in [0,1]*/
  inline double randomDouble() { return (double)gen() / gen.max(); }
  /*return a random double in [a,b], numbers can be gened*/
  inline double randomDouble(double a, double b) {
    return (b - a) * gen() / gen.max() + a;
  }
  /*return a random int in [a,b], (b-a+1) numbers can be gened*/
  inline int randomInt(int a, int b) { return gen() % (b - a + 1) + a; }
  inline double max(double a, double b) { return (a > b) ? a : b; }
  inline double min(double a, double b) { return (a < b) ? a : b; }
  inline void swap(double& a, double& b) {
    double temp = a;
    a = b;
    b = temp;
  }
};
