#include "optimizer/ga.h"
#include "some_func.h"
// #include "optimizer/optimizer.h"
using namespace MetaOpt;
using namespace std;

int main() {
  // (n_dim, n_obj, n_con, func, upper, lower)
  const int dim = 2;
  MetaOpt::Ga<double> a(dim, 1, 0, rosenbrock<double, dim>, vector<double>(dim,0),
  vector<double>(dim, 2), 10 * dim, 50 * dim);
  // MetaOpt::Ga<double> a(7, 1, 4, G9<double>, vector<double>(7, 0),
                        // vector<double>(7, 10), 100, 500);
  a.opt();

}
