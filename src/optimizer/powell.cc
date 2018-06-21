#include "powell.h"
#include <cmath>
#include <ctime>
#include "common.h"

using namespace std;

namespace MetaOpt {

template <typename Real>
void Powell<Real>::InitDirec(Real** direc) {
  for (unsigned int i = 0; i != this->n_x_; ++i) {
    for (unsigned int j = 0; j != this->n_x_; ++j) {
      direc[i][j] = 0;
    }
    direc[i][i] = 1.0;
  }
}
template <typename Real>
void Powell<Real>::InitReverseDirec(Real** direc) {
  for (unsigned int i = 0; i != this->n_x_; ++i) {
    for (unsigned int j = 0; j != this->n_x_; ++j) {
      direc[i][j] = 0;
    }
    direc[i][i] = -1.0;
  }
}
template <typename Real>
void Powell<Real>::InitRandomDirec(Real** direc) {
  for (unsigned int i = 0; i != this->n_x_; ++i) {
    for (unsigned int j = 0; j != this->n_x_; ++j) {
      direc[i][j] = exd::random::random<Real>(0, 1);
    }
  }
}

template <typename Real>
Real Powell<Real>::f1dim2(Real          alpha,
                          Sample<Real>& x,
                          Real*         p,
                          Sample<Real>& temp) {
  for (int j = 0; j != this->n_x_; ++j) {
    temp.x()[j] = x.x()[j] + alpha * p[j];
  }
  this->evaluate(temp);
  return temp.obj()[0];
}

template <typename Real>
Real Powell<Real>::Brent(Real                        xa,
                         Real                        xb,
                         Real                        xc,
                         Real                        tol,
                         Real&                       xmin,
                         std::function<Real(Real x)> func1d) {
  int        done, maxiter = 100;
  const Real mintol = 1.0e-11, cgold = 0.381966;
  Real       rat, fu, r, q, p, xmid, tol1, tol2, a, b;
  Real       u, etemp, dum, v, w, x, deltax, fx, fv, fw;
  x = w = v = xb;
  fw = fv = fx = func1d(x);
  if (xa < xc) {
    a = xa;
    b = xc;
  } else {
    a = xc;
    b = xa;
  }
  deltax = 0.0;
  for (int iter = 1; iter <= maxiter; iter++) {
    tol1 = tol * fabs(x) + mintol;
    tol2 = 2.0 * tol1;
    xmid = 0.5 * (a + b);
    if (fabs(x - xmid) <= tol2 - 0.5 * (b - a)) {
      break;
    }
    done = -1;
    if (fabs(deltax) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) {
        p = -p;
      }
      q = fabs(q);
      etemp = deltax;
      deltax = rat;
      dum = fabs(0.5 * q * etemp);
      if (fabs(p) < dum && p > q * (a - x) && p < q * (b - x)) {
        rat = p / q;
        u = x + rat;
        if (u - a < tol2 || b - u < tol2) {
          rat = fabs(tol1) * sgn(xmid - x);
        }
        done = 0;
      }
    }
    if (done) {
      if (x >= xmid) {
        deltax = a - x;
      } else {
        deltax = b - x;
      }
      rat = cgold * deltax;
    }
    if (fabs(rat) >= tol1) {
      u = x + rat;
    } else {
      u = x + fabs(tol1) * sgn(rat);
    }
    fu = func1d(u);
    if (fu <= fx) {
      if (u >= x) {
        a = x;
      } else {
        b = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    } else {
      if (u < x) {
        a = u;
      } else {
        b = u;
      }
      if (fu <= fw || w == x) {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      } else {
        if (fu <= fv || v == x || v == w) {
          v = u;
          fv = fu;
        }
      }
    }
  }
  xmin = x;
  return fx;
}

template <typename Real>
void Powell<Real>::LineSearch(Sample<Real>& x, Real* direc) {
  int                         j;
  Real                        tol = 1e-4;
  Real                        fa, f, fb, xb, xa = 0.0;
  Real                        xmin, xx = 1.0;
  Sample<Real>                tempx = this->sample();
  std::function<Real(Real x)> func1d =
      std::bind(&Powell<Real>::f1dim2, this, placeholders::_1, x, direc, tempx);
  Mnbrak(xa, xx, xb, fa, f, fb, func1d);
  x.obj()[0] = Brent(xa, xx, xb, tol, xmin, func1d);
  for (j = 0; j != this->n_x_; ++j) {
    direc[j] = xmin * direc[j];
    x.x()[j] = x.x()[j] + direc[j];
  }
}

template <typename Real>
int Powell<Real>::evolve(Sample<Real>& x0,
                         Real**        direc,
                         int           maxiter,
                         Real          ftol,
                         int           iter,
                         Real          terminalLine) {
  this->evaluate(x0);
  int          ret = 0;
  Sample<Real> x1 = this->sample();
  Sample<Real> x2 = this->sample();
  Real*        direc1 = new Real[this->n_x_];
  int          bigind;
  Real         t, temp, fx, delta, fx2;
  Real&        fval = x0.obj()[0];
  while (1) {
    x1 = x0;
    fx = fval;
    bigind = 0;
    delta = 0.0;
    for (unsigned int i = 0; i < this->n_x_; i++) {
      fx2 = fval;
      LineSearch(x0, direc[i]);
      if (fval < terminalLine) {
        break;
      }
      if (fabs(fx2 - fval) > delta) {
        delta = fabs(fx2 - fval);
        bigind = i;
      }
      // cout<<" "<<i<<" vec search"<<fval<<endl;
    }
    if (fval < terminalLine) {
      ret = iter;
      break;
    }
    // Construct the extrapolated point
    for (unsigned int j = 0; j != this->n_x_; ++j) {
      direc1[j] = x0.x()[j] - x1.x()[j];
      x2.x()[j] = x0.x()[j] + direc1[j];
    }
    this->evaluate(x2);
    fx2 = x2.obj()[0];
    if (fx > fx2) {
      t = 2.0 * (fx + fx2 - 2.0 * fval);
      temp = fx - fval - delta;
      t *= temp * temp;
      temp = fx - fx2;
      t -= delta * temp * temp;
      if (t < 0.0) {
        LineSearch(x0, direc1);
        if (fval < terminalLine) {
          ret = iter;
          break;
        }
        for (unsigned int j = 0; j != this->n_x_; ++j) {
          direc[bigind][j] = direc[this->n_x_ - 1][j];
          direc[this->n_x_ - 1][j] = direc1[j];
        }
      }
    }
    iter = iter + 1;
    // LOG(INFO) << "Iter = " << iter << ", func = " << fval << ", n_evaluation
    // = " << this->n_evaluation_;
    if (fabs(fx - fval) < ftol) {
      ret = iter;
      break;
    }
    // if (2.0*fabs(fx-fval)<=ftol*(fabs(fx)+fabs(fval))+1e-20){
    //	break;
    //}
    if (iter >= maxiter) {
      // cout<<"powell exceeding maximum iterations"<<"  "<<ret<<endl;
      ret = iter;
      break;
    }
  }
  delete[] direc1;
  return ret;
}

template <typename Real>
void Powell<Real>::opt(Sample<Real>& x0,
                       Real          ftol,
                       int           maxiter,
                       int           iter,
                       Real          terminalLine) {
  Real** direc = new Real*[this->n_x_];
  direc[0] = new Real[this->n_x_ * this->n_x_];
  for (unsigned int i = 1; i != this->n_x_; ++i)
    direc[i] = direc[i - 1] + this->n_x_;
  Real fbest;
  bool flip = true;
  for (int i = 0; i < 100; ++i) {
    if (flip)
      InitDirec(direc);
    else
      InitReverseDirec(direc);
    flip = !flip;
    iter = evolve(x0, direc, maxiter, ftol, iter, terminalLine);
    LOG(INFO) << "Iter: " << iter << ", best: " << x0
              << ", n_evaluation: " << this->n_evaluation_;
    if (fabs(x0.obj()[0] - fbest) < ftol || iter >= maxiter) {
      break;
    }
    fbest = x0.obj()[0];
  }

  delete[] direc[0];
  delete[] direc;
  return;
}

template <typename Real>
void Powell<Real>::Mnbrak(Real&                       xa,
                          Real&                       xb,
                          Real&                       xc,
                          Real&                       fa,
                          Real&                       fb,
                          Real&                       fc,
                          std::function<Real(Real x)> func1d) {
  Real r, q, temp, gold = 1.618034;
  Real maxiter = 100;
  int  glimit = 100;
  Real u, ulim, fu, tiny = 1e-30;
  fa = func1d(xa);
  fb = func1d(xb);
  if (fb > fa) {
    temp = xa;
    xa = xb;
    xb = temp;
    temp = fb;
    fb = fa;
    fa = temp;
  }
  xc = xb + gold * (xb - xa);
  fc = func1d(xc);
  int iter = 0;
  while (fb >= fc) {
    iter++;
    if (iter < maxiter) break;

    r = (xb - xa) * (fb - fc);
    q = (xb - xc) * (fb - fa);
    temp = q - r;
    // cout<<"xa "<<xa<<", xb "<<xb<<", xc "<<xb<<endl;
    // cout<<"fa "<<fa<<", fb "<<fb<<", fc "<<fb<<endl;
    // if(xa>1e10)
    //	system("pause");
    // cout<<"r "<<r<<", q "<<q<<", temp "<<temp<<endl;
    if (fabs(temp) < tiny) {
      temp = tiny;
    }
    u = xb - ((xb - xc) * q - (xb - xa) * r) / (2 * temp);
    ulim = xb + glimit * (xc - xb);
    if ((xb - u) * (u - xc) > 0) {
      fu = func1d(u);
      if (fu < fc) {
        xa = xb;
        fa = fb;
        xb = u;
        fb = fu;
        return;
      } else {
        if (fu > fb) {
          xc = u;
          fc = fu;
          return;
        }
      }
      u = xc + gold * (xc - xb);
      fu = func1d(u);
    } else {
      if ((xc - u) * (u - ulim) > 0) {
        fu = func1d(u);
        if (fu < fc) {
          xb = xc;
          xc = u;
          u = xc + gold * (xc - xb);
          fb = fc;
          fc = fu;
          fu = func1d(u);
        }
      } else {
        if ((u - ulim) * (ulim - xc) >= 0) {
          u = ulim;
          fu = func1d(u);
        } else {
          u = xc + gold * (xc - xb);
          fu = func1d(u);
        }
      }
    }
    xa = xb;
    xb = xc;
    xc = u;
    fa = fb;
    fb = fc;
    fc = fu;
  }
}
template class Powell<double>;
}  // namespace MetaOpt
