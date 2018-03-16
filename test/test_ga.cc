#include "optimizer/ga.h"
// #include "optimizer/optimizer.h"
// using namespace MetaOpt;

void rosenbrock(double* x, double* y, double* c) {
  y[0] = (100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
          (1 - x[0]) * (1 - x[0]));
}

using namespace std;
int main(){
  // (n_dim, n_obj, n_con, func, upper, lower)
  MetaOpt::Ga<double> a(2, 1, 0, rosenbrock, vector<double>{-2, -2},
                        vector<double>{2, 2}, 60);
  a.opt();
  // auto s = a.sample();
  // std::cout << s << std::endl;
  // a.evaluate(s);
  // std::cout << s << std::endl;
  // MetaOpt::Optimizer<double> a(2, 1, 0, rosenbrock, vector<double>{0, 0},
  // vector<double>{2, 0});
  // MetaOpt::<float> a(2);
  return 0;
}
