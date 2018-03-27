#include "optimizer/powell.h"
#include "some_func.h"
using namespace MetaOpt;
using namespace std;

int main() {
  // (n_dim, n_obj, n_con, func, upper, lower)
  const int dim = 2;
  MetaOpt::Powell<double> a(dim, rosenbrock<double, dim>);
  auto s = a.sample();
  a.opt(s);
  // auto s = a.sample();
  std::cout << s << std::endl;
  std::cout << a.n_evaluation_ << std::endl;
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
