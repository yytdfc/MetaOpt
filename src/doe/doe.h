#ifndef METAOPT_DOE_H
#define METAOPT_DOE_H

#include <string>
#include <vector>

template <typename Real>
class Doe
{
 public:
  Doe(int n, int num);
  ~Doe();
  void gen();
  void genMC();
  void readSample(std::string str);
  void addLHS(int num);
  void printSample(std::string str);
  void printAdd(std::string str);
  Real** sample;

 private:
  int    addNum;
  int    ndim;
  int    nSample;
  Real** addSample;
  bool   ifadd;
  bool   hasSample;
};

#endif  // METAOPT_DOE_H