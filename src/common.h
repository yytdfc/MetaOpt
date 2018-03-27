#ifndef METAOPT_COMMON_H
#define METAOPT_COMMON_H
// std
#include <functional>
#include <iostream>
#include <string>
#include <vector>
// glog
#include <glog/logging.h>
// random
#include "random.hpp"
using Random = effolkronium::random_static;

// ostream
template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
  s.put('[');
  char comma[3] = {'\0', ' ', '\0'};
  for (const auto& e : v) {
    s << comma << e;
    comma[0] = ',';
  }
  return s << ']';
}

#endif  // METAOPT_COMMON_H
