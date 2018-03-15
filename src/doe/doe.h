#include <iostream>
#include <random>
#include <string>
using namespace std;
class Sample
{
 public:
  Sample(int n, int num);
  ~Sample();
  void genLHS();
  void genMC();
  void readSample(string str);
  void addLHS(int num);
  void printSample(string str);
  void printAdd(string str);
  double** sample;

 private:
  int                addNum;
  int                ndim;
  int                nSample;
  double**           addSample;
  bool               ifadd;
  bool               hasSample;
  std::random_device rd;
  std::mt19937_64    gen;
  /*return a random double in [0,1]*/
  inline double randomDouble() { return ((double)gen() / gen.max()); }
  /*return a random int in [a,b], (b-a+1) numbers can be gened*/
  inline int randomInt(int a, int b) { return (gen() % (b - a + 1)) + a; }
};
