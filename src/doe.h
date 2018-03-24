#ifndef METAOPT_DOE_H
#define METAOPT_DOE_H

#include "sample.h"

#include "random.hpp"
using Random = effolkronium::random_static;

namespace MetaOpt {
namespace Doe{
  template <typename Real>
Samples<Real> gen(int n_x, int n_dim, int n_obj, int n_con, int method=0) {
  Samples<Real> samples(n_x);
  for(auto& s: samples){
    s = Sample<Real>(n_dim, n_obj, n_con);
  }
  switch(method){
    case 0:
      std::vector<int> list(n_x);
      int i = 0;
      for (auto& t : list)
        t = i++;
      for (int d = 0; d < n_dim; ++d) {
        Random::shuffle(list);
        for(int i=0;i!=n_x;++i){
          samples[i].x()[d] = (Random::get<Real>(0, 1) + list[i] )/ n_x;

        }
      }

      // int* list = new int[n_sample_];
      // int  rand;
      // 
      //   for (j = 0; j < n_sample_; ++j) {
      //     list[j] = j;
      //   }
      //   
      // }
      // hasSample = 1;
      // delete[] list;
      break;
  }
  return samples;
}
}


}
#endif  // METAOPT_DOE_H