#include <gtest/gtest.h>
#include "optimizer/powell.h"
#include "doe.h"
#include "some_func.h"
using namespace MetaOpt;
using namespace std;

TEST(metaopt, doe) {
  const int               dim = 2;
  MetaOpt::Powell<double> a(dim, rosenbrock<double, dim>);
  auto r = Doe<double, MC>::gen(1, dim, vector<double>(dim, 0),
                                vector<double>(dim, 2))[0];
  std::cout << r << std::endl;
  a.opt(r);
  std::cout << r << std::endl;
  std::cout << a.n_evaluation_ << std::endl;
  EXPECT_LE(fabs(r[dim]), 1e-6);
  for (auto i : exd::range::range(dim)) {
    EXPECT_LE(fabs(r[i] - 1), 1e-3);
  }
}
