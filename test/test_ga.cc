#include "optimizer/ga.h"
// #include "optimizer/optimizer.h"
using namespace MetaOpt;

void rosenbrock(double *x, double *y, double *c) {
  y[0] = (100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
          (1 - x[0]) * (1 - x[0]));
}

void G9(double *x, double *y, double *c) {
  c[0] = -(2 * pow(x[0], 2) + 3 * pow(x[1], 4) + x[2] + 4 * pow(x[3], 2) +
           5 * x[4] - 127);
  c[1] = -(7 * x[0] + 3 * x[1] + 10 * pow(x[2], 2) + x[3] - x[4] - 282);
  c[2] = -(23 * x[0] + pow(x[1], 2) + 6 * pow(x[5], 2) - 8 * x[6] - 196);
  c[3] = -(4 * pow(x[0], 2) + pow(x[1], 2) - 3 * x[0] * x[1] +
           2 * pow(x[2], 2) + 5 * x[5] - 11 * x[6]);
  y[0] = pow((x[0] - 10), 2) + 5 * pow((x[1] - 12), 2) + pow(x[2], 4) +
         3 * pow((x[3] - 11), 2) + 10 * pow(x[4], 6) + 7 * pow(x[5], 2) +
         pow(x[6], 4) - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6];
}

namespace MetaOpt {
template class Ga<double>;
}

using namespace std;

int main() {
  // (n_dim, n_obj, n_con, func, upper, lower)
  // MetaOpt::Ga<double> a(2, 1, 0, rosenbrock, vector<double>{-2, -2},
  // vector<double>{2, 2}, 80, 100);
  MetaOpt::Ga<double> a(7, 1, 4, G9, vector<double>(7, 0),
                        vector<double>(7, 10), 100, 500);
  a.opt();
  // auto s = a.sample();
  // std::cout << s << std::endl;
  // a.evaluate(s);
  // std::cout << s << std::endl;
  // MetaOpt::Ga<double> a(2, 1, 0, rosenbrock, vector<double>{0, 0},
  // vector<double>{2, 0}, 10,0.1,0.1);
  // MetaOpt::<float> a(2);

  // MetaOpt::Optimizer<double> a(2, 1, 0, rosenbrock, vector<double>{-2, -2},
  //                         vector<double>{2, 2});

  // MetaOpt::Sample<double> a(vector<double>{2,
  // 2},vector<double>(2,0),vector<double>());
  // // MetaOpt::Sample<double> aa(vector<double>{2,
  // 2},vector<double>(2,0),vector<double>()); MetaOpt::Sample<double>
  // aa(move(vector<double>{2,
  // 2}),move(vector<double>(2,0)),move(vector<double>()));
  // // MetaOpt::Sample<double> ab(vector<double>{2,
  // 2},vector<double>(2,0),vector<double>()); LOG(INFO) << a; LOG(INFO) << aa;
  // MetaOpt::Sample<double> b(a);
  // MetaOpt::Sample<double> d(move(aa));
  // LOG(INFO) << a;
  // LOG(INFO) << aa;
  // LOG(INFO) << a.x_.size();
  // LOG(INFO) << aa.x_.size();
  // a.x_[0] = 3;
  // // aa.x_[0] = 3;
  // LOG(INFO) << a;
  // LOG(INFO) << aa;
  // LOG(INFO) << "(move(Sample()))";
  // MetaOpt::Sample<double> dss(move(MetaOpt::Sample<double>(
  //     vector<double>{2, 2}, vector<double>(2, 0), vector<double>())));
  // LOG(INFO) << "(Sample())";
  // MetaOpt::Sample<double> cc(vector<double>{2, 3}, vector<double>(2, 0),
  //                            vector<double>());
  // MetaOpt::Sample<double> dsss(cc);
  // cc.x_[0] = 3;
  // LOG(INFO) << dss;
  // LOG(INFO) << cc;
  // LOG(INFO) << dsss;
  // cc = dss;
  // LOG(INFO) << dss;
  // LOG(INFO) << cc;
  // LOG(INFO) << dsss;
  // dss.x_[0] = 9;
  // LOG(INFO) << dss;
  // LOG(INFO) << cc;
  // LOG(INFO) << dsss;
  // cc = move(dss);
  // LOG(INFO) << dss;
  // LOG(INFO) << cc;
  // LOG(INFO) << dsss;
  // return 0;
}
