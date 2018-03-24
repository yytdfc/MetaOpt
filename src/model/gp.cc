#include "gp.h"

#include <eigen3/Eigen/Cholesky>
#include <eigen3/Eigen/LU>
#include "common.h"
#include "kernel.h"
namespace MetaOpt {
using namespace Eigen;
using namespace std;

template <typename Real>
void Gp<Real>::init(const Samples<Real>& samples) {
  Model<Real>::init(samples);

  // length_scale_
  length_scale_ = 1.0;
  // mu_
  mu_ = 1e-10;
  // k = | C F |
  //     | F 0 |
  switch (i_order_) {
    case 0:
      f_dim_ = 0;
      break;
    case 1:
      f_dim_ = 1;
      break;
    case 2:
      f_dim_ = 1 + this->n_dim_;
      break;
  }
  k_ = Matrixx<Real>(this->n_sample_ + f_dim_, this->n_sample_ + f_dim_);
  switch (i_kernel_) {
    case GAUSS:
      kernel_fun_ = Kernel<Real, GAUSS>::correlation;
      break;
    case MATERN32:
      kernel_fun_ = Kernel<Real, MATERN32>::correlation;
      break;
    case MATERN52:
      kernel_fun_ = Kernel<Real, MATERN52>::correlation;
      break;
  }
};

template <typename Real>
void Gp<Real>::fit(const Samples<Real>& samples) {
  if (samples.empty()) {
    if (this->samples_.empty()) {
      LOG(ERROR) << "No input samples!";
    }
  } else {
    init(samples);
  }
  Real sigma;
  obj_mle(&length_scale_, &sigma);
  LOG(INFO) << "MLE: " << sigma;
};

template <typename Real>
void Gp<Real>::obj_mle(Real* length_scale, Real* obj, Real* con) {
  // correlation
  for (int i = 0; i != this->n_sample_; ++i) {
    k_(i, i) = 1 + mu_;
    for (int j = i + 1; j != this->n_sample_; ++j) {
      k_(j, i) = k_(i, j) = kernel_fun_(this->samples_[i].x(),
                                        this->samples_[j].x(), *length_scale);
    }
  }
  // init F
  // i_order_ = 1
  if (i_order_ >= 1) {
    for (int i = 0; i != this->n_sample_; ++i) {
      k_(i, this->n_sample_) = k_(this->n_sample_, i) = 1;
    }
    k_(this->n_sample_, this->n_sample_) = mu_;
  }
  // // i_order_ = 2
  if (i_order_ == 2) {
    for (int i = 0; i != this->n_sample_; ++i) {
      for (int j = 0; j < this->n_dim_; j++) {
        k_(i, this->n_sample_ + j + 1) = k_(this->n_sample_ + j + 1, i) =
            this->samples_[i].x()[j];
        k_(this->n_sample_ + j + 1, this->n_sample_ + j + 1) = mu_;
      }
    }
  }
  LOG(INFO) << "K\n" << k_;
  Matrixx<Real> y(this->n_sample_ + f_dim_, 1);
  for (int i = 0; i != this->n_sample_; ++i) {
    y(i, 0) = this->samples_[i].obj()[0];
  }
  ldlt_.reset(new LDLT<Ref<Matrixx<Real>>>(k_));
  alpha_ = ldlt_->solve(y);
  LOG(INFO) << "K\n" << k_;
  LOG(INFO) << "alpha_\n" << alpha_;
  Matrixx<Real> diag = ldlt_->vectorD();
  Real log_diag = 0.0;
  for (int i = 0; i != this->n_sample_ + f_dim_; ++i) {
    log_diag += log(fabs(diag(i, 0)));
  }
  sigma_ = (y.transpose() * alpha_)(0, 0) / this->n_sample_;
  obj[0] = sigma_ + log_diag;
  LOG(INFO) << "sigma: " << sigma_;
  return;
}

template <typename Real>
void Gp<Real>::evaluate(Sample<Real>& sample, vector<Real>& mse) {
  if (this->i_norm_) {
    for (int i = 0; i != this->n_dim_; ++i) {
      sample.x()[i] = (sample.x()[i] - this->mean_[i]) / this->st_[i];
    }
  }
  Matrixx<Real> k(this->n_sample_ + f_dim_, 1);
  for (int i = 0; i != this->n_sample_; ++i) {
    k(i, 0) = kernel_fun_(sample.x(), this->samples_[i].x(), length_scale_);
  }
  if (i_order_ >= 1) {
    k(this->n_sample_, 0) = 1;
  }
  if (i_order_ == 2) {
    for (int i = 0; i < this->n_dim_; i++) {
      k(this->n_sample_ + i + 1, 0) = sample.x()[i];
    }
  }
  sample.obj()[0] = (k.transpose() * alpha_)(0, 0);
  if (!mse.empty()) {
    mse[0] = sigma_ * sqrt(fabs((1 - (k.transpose() * ldlt_->solve(k))(0, 0))));
  }
};

template class Gp<double>;
}  // namespace MetaOpt
