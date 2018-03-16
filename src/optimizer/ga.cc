#include "ga.h"
using namespace std;
namespace MetaOpt {

template class Ga<double>;

template <typename Real>
Ga<Real>::Ga(const int               n_dim,
             const int               n_obj,
             const int               n_con,
             const Func<Real>        func,
             const std::vector<Real> upper,
             const std::vector<Real> lower,
             const int               n_population,
             const Real              p_crossover,
             const Real              p_mutation)
    : Optimizer<Real>(n_dim, n_obj, n_con, func, upper, lower) {
  n_population_ = n_population > 0 ? n_population : 10 * n_dim;
  n_generation_ = 10 * n_dim;
  vector<int>(n_population_).swap(tourlist_);
  int i = 0;
  for (auto& t : tourlist_)
    t = i++;
  vector<Sample<Real>>(n_population_).swap(population_);
  for (auto& p : population_)
    this->sample().swap(p);
  vector<Sample<Real>>(n_population_).swap(newpop_);
  for (auto& p : newpop_)
    this->sample().swap(p);
  if (n_con > 0) vector<Real>(n_dim).swap(cons_p_);
  p_crossover_ = 1 - pow(1 - p_crossover, 1.0 / n_dim);
  p_mutation_ = 1 - pow(1 - p_mutation, 1.0 / n_dim);
  LOG(INFO) << "Ga constructed.";
}
template <typename Real>
Ga<Real>::~Ga() {
  LOG(INFO) << "Ga destructed";
}
template <typename Real>
void Ga<Real>::opt() {
  for (i_generation_ = 0; i_generation_ < n_generation_; ++i_generation_) {
    statistics();

    Random::shuffle(tourlist_.begin(), tourlist_.end());
    {
      LOG(INFO) << "stat";
      for (int i = 0; i < n_population_; ++i)
        LOG(INFO) << "population[" << i << "]:  " << population_[i];
    }
    for (int k = 0; k != n_population_; ++k) {
      select(population_[k], population_[tourlist_[k]], newpop_[k]);
      // if (k % 2 == 1) crossover(newpop_[k], newpop_[k - 1]);
    }
    {
      LOG(INFO) << "select";
      for (int i = 0; i < n_population_; ++i)
        LOG(INFO) << "population[" << i << "]:  " << newpop_[i];
    }
    for (int k = 0; k != n_population_; k += 2) {
      crossover(newpop_[k], newpop_[k + 1]);
    }
    {
      LOG(INFO) << "crossover";
      for (int i = 0; i < n_population_; ++i)
        LOG(INFO) << "population[" << i << "]:  " << newpop_[i];
    }
    for (auto& p : newpop_) {
      mutation(p);
    }
    {
      LOG(INFO) << "mutation:";
      for (int i = 0; i < n_population_; ++i)
        LOG(INFO) << "population[" << i << "]:  " << newpop_[i];
    }
    for (int k = 0; k < n_population_; ++k) {
      population_[k].swap(newpop_[k]);
    }
    LOG(INFO) << "xxxx:";
    for (int i = 0; i < n_population_; ++i)
      LOG(INFO) << "population[" << i << "]:  " << population_[i];
    if (i_generation_ % 1 == 0) {
      LOG(INFO) << "Generation " << i_generation_ + 1
                << ", best: " << this->best_;
    }
  }
  LOG(INFO) << "GA evolve finished.";
  LOG(INFO) << "best:  " << this->best_;

  LOG(INFO) << "nMutation = " << i_mutation_;
  LOG(INFO) << "nCrossover = " << i_crossover_;
}

template <typename Real>
void Ga<Real>::statistics() {
  bool flag = 1;
  int  bestn = -1;
  auto best = numeric_limits<Real>::infinity();
  for (int i = 0; i < n_population_; ++i) {
    this->evaluate(population_[i]);

    // if (n_con_ > 0) {
    //   flag = 1;
    //   cons_p_[i] = 0;
    //   for (int j = 0; j < n_con_; ++j) {
    //     if (population_[i][n_dim_ + j + 1] < 0) {
    //       cons_p_[i] -= population_[i][n_dim_ + j + 1];
    //       flag = 0;
    //     }
    //   }
    // }
    if (best > std::get<1>(population_[i])[0] && flag){
      bestn = i;
      best = std::get<1>(population_[i])[0];
    }
  }
  if (bestn != -1) {
    this->best_ = population_[bestn];
  }
}

template <typename Real>
void Ga<Real>::select(Sample<Real>& s1, Sample<Real>& s2, Sample<Real>& s) {
  // LOG(INFO) << "select:" << std::get<1>(s1)[0] << ":" << std::get<1>(s2)[0];
  // s.swap(std::get<1>(s1)[0] < std::get<1>(s2)[0] ? s1 : s2);
  auto& ss = std::get<1>(s1)[0] < std::get<1>(s2)[0] ? s1 : s2;
  for (int i = 0; i != this->n_dim_; ++i) {
    std::get<0>(s)[i] = std::get<0>(ss)[i];
  }
  std::get<1>(s)[0] = std::get<1>(ss)[0];
  // LOG(INFO) << "   get:" << std::get<1>(s)[0];
}

template <typename Real>
void Ga<Real>::crossover(Sample<Real>& s1, Sample<Real>& s2) {
  switch (1) {
    case 0:
      for (int i = 0; i != this->n_dim_; ++i) {
        if (Random::get<Real>(0, 1) < p_crossover_) {
          Real dif = std::get<0>(s2)[i] - std::get<0>(s1)[i];
          std::get<0>(s2)[i] =
              std::get<0>(s1)[i] + Random::get<Real>(0, 1) * dif;
          std::get<0>(s1)[i] =
              std::get<0>(s1)[i] + Random::get<Real>(0, 1) * dif;
          ++i_crossover_;
        }
      }
      break;
    case 1:
      for (int i = 0; i != this->n_dim_; ++i) {
        if (Random::get<Real>(0, 1) < p_crossover_) {
          auto& x1i = std::get<0>(s1)[i];
          auto& x2i = std::get<0>(s2)[i];
          if (x1i > x2i) swap(x1i, x2i);
          Real mid = (x2i + x1i) / 2;
          Real dif = x2i - x1i;
          if (dif < 1e-6) {
            dif = 1e-6;
          }
          Real beta =
              1 + 2 / dif * min(x1i - this->lower_[i], this->upper_[i] - x2i);
          Real alpha = 2 - pow(beta, -(nc + 1));
          if (alpha > 2 - 1e-6) alpha = 2 - 1e-6;
          alpha *= Random::get<Real>(0, 1);
          if (alpha <= 1)
            beta = pow(alpha, 1.0 / (nc + 1));
          else
            beta = 1.0 / pow(2 - alpha, 1.0 / (nc + 1));

          x1i = mid - 0.5 * beta * dif;
          x2i = mid + 0.5 * beta * dif;
          if (isnan(x1i) || isinf(x1i))
            x1i = Random::get<Real>(this->lower_[i], this->upper_[i]);
          if (x1i > this->upper_[i]) x1i = this->upper_[i];
          if (x1i < this->lower_[i]) x1i = this->lower_[i];
          if (isnan(x2i) || isinf(x2i))
            x2i = Random::get<Real>(this->lower_[i], this->upper_[i]);
          if (x2i > this->upper_[i]) x2i = this->upper_[i];
          if (x2i < this->lower_[i]) x2i = this->lower_[i];
          ++i_crossover_;
        }
      }
      break;
  }
  
}

template <typename Real>
void Ga<Real>::mutation(Sample<Real>& s) {
  switch (0) {
    case 0:
      for (int i = 0; i != this->n_dim_; ++i) {
        if (Random::get<Real>(0, 1) < p_mutation_) {
          std::get<0>(s)[i] =
              Random::get<Real>(this->lower_[i], this->upper_[i]);
          ++i_mutation_;
        }
      }
    case 1:
      for (int i = 0; i != this->n_dim_; ++i) {
        if (Random::get<Real>(0, 1) < p_mutation_) {
          auto& xi = std::get<0>(s)[i];
          Real  u = Random::get<Real>(0, 1);
          Real  deltamax = this->upper_[i] - this->lower_[i];
          Real  delta =
              min(xi - this->lower_[i], this->upper_[i] - xi) / deltamax;
          if (u > 0.5) {
            u = 1 - u;
            delta = 1 - pow((2 * u + (1 - 2 * u) * pow(1 - delta, nm + 1)),
                            1 / (nm + 1));
          } else {
            delta = pow((2 * u + (1 - 2 * u) * pow(1 - delta, nm + 1)),
                        1 / (nm + 1)) -
                    1;
          }
          xi += delta * deltamax;
          if (isnan(xi) || isinf(xi))
            xi = Random::get<Real>(this->lower_[i], this->upper_[i]);
          if (xi > this->upper_[i])
            xi = this->upper_[i];
          else if (xi < this->lower_[i])
            xi = this->lower_[i];
          ++i_mutation_;
        }
      }
      break;
  }
}
}