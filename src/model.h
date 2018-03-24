#ifndef METAOPT_MODEL_H
#define METAOPT_MODEL_H

#include "common.h"
#include "sample.h"

namespace MetaOpt {

template <typename Real>
using Func = std::function<void(Real* x, Real* y, Real* c)>;

template <typename Real>
class Model
{
 public:
  Model(){};
  Model(const Samples<Real>& samples){
    init(samples);
  };
  ~Model(){};
  void init(const Samples<Real>& samples){
    if(samples.empty())
      return;
    samples_ = samples;
    n_sample_ = samples_.size();
    n_dim_ = samples_[0].x().size();
    n_y_ = samples_[0].obj().size() + samples_[0].con().size();
    // norm
    mean_ = std::vector<Real> (n_dim_, 0);
    st_ = std::vector<Real> (n_dim_, 1);
    // for(auto& s :samples_){
    //   for(int i = 0; i!=n_dim_;++i){
    //     mean_[i]+=s.x()[i];
    //   }
    // }
    // for(int i = 0; i!=n_dim_;++i){
    //   mean_[i] /= n_sample_;
    // }
    // for(auto& s :samples_){
    //   for(int i = 0; i!=n_dim_;++i){
    //     st_[i]+=pow(s.x()[i] - mean_[i], 2);
    //   }
    // }
    // for(int i = 0; i!=n_dim_;++i){
    //   st_[i] = sqrt(st_[i]) / n_sample_;
    //   if(st_[i] <= std::numeric_limits<Real>::min())
    //     st_[i] = std::numeric_limits<Real>::min();
    // }
    for(auto& s :samples_){
      for(int i = 0; i!=n_dim_;++i){
        s.x()[i] = (s.x()[i] - mean_[i]) / st_[i];
      }
    }

    
    // for(int i = 0; i!=n_dim_;++i){
    //   upper_[i] = mean_[i] + 0.5 * st_[i];
    //   lower_[i] = mean_[i] - 0.5 * st_[i];
    // }

    // upper_ = std::vector<Real>(n_dim_, - std::numeric_limits<Real>::max());
    // lower_ = std::vector<Real>(n_dim_,   std::numeric_limits<Real>::max());
    // for(auto& s :samples_){
    //   for(int i = 0; i!=n_dim_;++i){
    //     upper_[i]=std::max(upper_[i], s.x()[i]);
    //     lower_[i]=std::min(lower_[i], s.x()[i]);
    //   }
    // }
    // for(int i = 0; i!=n_dim_;++i){
    //   if(upper_[i]==lower_[i]){
    //     upper_[i]=lower_[i] + std::numeric_limits<Real>::min();
    //     lower_[i]=lower_[i] - std::numeric_limits<Real>::min();
    //   }
    // }
    // if(i_norm_){
      // for(auto& s :samples_){
      //   for(int i = 0; i!=n_dim_;++i){
      //     s.x()[i] = (s.x()[i] - lower_[i]) / (upper_[i] - lower_[i]);
      //   }
      // }
    // }
  };
  void              evaluate(Sample<Real>& sample);
  void              evaluate(Samples<Real>& samples, int n_thread = 1);
  void              fit(const Samples<Real>& samples = nullptr);
  int               n_dim_;
  int               n_y_;
  int               n_sample_;
  std::vector<Real> upper_;
  std::vector<Real> lower_;
  std::vector<Real> mean_;
  std::vector<Real> st_;
  Samples<Real>     samples_;
  int               i_norm_ = true;
};
}  // namespace MetaOpt

#endif // METAOPT_MODEL_H
