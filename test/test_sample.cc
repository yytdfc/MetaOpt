#include <gtest/gtest.h>
#include "some_func.h"
#include "doe.h"
using namespace MetaOpt;

TEST(metaopt, doe) {
  const std::size_t n = 10;
  auto              lower = exd::range::range<double>(n).to_vector();
  auto              upper = exd::range::range<double>(1, n + 1).to_vector();
  {
    auto ss = Doe<double, LHS>::gen(n, n, lower, upper);
    EXPECT_EQ(ss.size(), n);
    for (auto& s : ss) {
      for (auto i : exd::range::range(n)) {
        EXPECT_GE(s[i], lower[i]);
        EXPECT_LE(s[i], upper[i]);
      }
    }
  }
  {
    auto ss = Doe<double, MC>::gen(n, n, lower, upper);
    EXPECT_EQ(ss.size(), n);
    for (auto& s : ss) {
      for (auto i : exd::range::range(n)) {
        EXPECT_GE(s[i], lower[i]);
        EXPECT_LE(s[i], upper[i]);
      }
    }
  }
  {
    auto ss = Doe<double, ZERO>::gen(n, n, lower, upper);
    EXPECT_EQ(ss.size(), n);
    for (auto& s : ss) {
      for (auto i : exd::range::range(n)) {
        EXPECT_GE(s[i], lower[i]);
        EXPECT_LE(s[i], upper[i]);
      }
    }
  }
}
