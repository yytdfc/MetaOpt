#include "optimizer/bayesopt.h"
#include "some_func.h"
using namespace MetaOpt;
using namespace std;

int main() {
  // (n_dim, n_obj, n_con, func, upper, lower)
  const int dim = 2;
  MetaOpt::BayesOpt<double> a(dim, 1, 0, rosenbrock<double, dim>, vector<double>(dim, 0),
      vector<double>(dim, 2), 12);
  a.opt();
  // auto s = a.sample();
  // std::cout << s << std::endl;
  std::cout << a.n_evaluation_ << std::endl;

  // cout << Doe<double, LHS>::gen(10, dim, vector<double>(dim, 0),
      // vector<double>(dim, 2), 1, 0) << endl;

  return 0;
}
