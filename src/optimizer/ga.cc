#include "ga.h"
#include "common.h"
#include <fstream>
using namespace std;

template class GA<double>;

template <typename Real>
void GA<Real>::setFx(std::function<Real(Real*)> f) {
  fx = f;
}
template <typename Real>
void GA<Real>::GAinit(int   n,
                      int   Pop,
                      int   Gens,
                      int   Cons,
                      Real  pCr,
                      Real  pMu,
                      Real* up,
                      Real* low) {
  int i;
  nPop = Pop;
  nGens = Gens;
  nVar = n;
  nCons = Cons;
  tourlist = new int[2 * nPop];
  pop = new Real*[nPop];
  pop[0] = new Real[nPop * (nVar + nCons + 1)];
  for (i = 1; i < nPop; ++i) {
    pop[i] = pop[i - 1] + nVar + nCons + 1;
  }
  newpop = new Real*[nPop];
  newpop[0] = new Real[nPop * nVar];
  for (i = 1; i < nPop; ++i) {
    newpop[i] = newpop[i - 1] + nVar;
  }
  upper = new Real[nVar];
  lower = new Real[nVar];
  for (int i = 0; i < nVar; ++i) {
    upper[i] = up[i];
    lower[i] = low[i];
  }
  best = new Real[nVar + nCons + 1];
  if (nCons > 0) {
    consP = new Real[nPop];
  }
  pCrossover = 1 - pow(1 - pCr, 1.0 / nVar);
  pMutation = 1 - pow(1 - pMu, 1.0 / nVar);
  nm = 1, nc = 1;
  cout << "GA constructed." << endl;
}
template <typename Real>
GA<Real>::~GA() {
  delete[] tourlist;
  delete[] pop[0];
  delete[] pop;
  delete[] newpop[0];
  delete[] newpop;
  delete[] upper;
  delete[] lower;
  delete[] best;

  if (nCons > 0) {
    delete[] consP;
  }
  // cout << "GA destructed." << endl;
}

template <typename Real>
void GA<Real>::evolve() {
  int i, j, k;
  initialpop();
  for (generation = 0; generation < nGens; ++generation) {
    tourney();
    for (k = 0; k < nPop; k += 2) {
      tourSelect(k);
      SBXover(k);
      PBMutation(k);
      PBMutation(k + 1);
    }
    for (i = 0; i < nPop; ++i) {
      for (j = 0; j < nVar; ++j) {
        pop[i][j] = newpop[i][j];
      }
    }
    statistics();
    // if (generation % 20 == 0){
    //	cout << "Generation. " << generation + 1 <<
    //		"\t best fit = \t" << best[nVar] << endl;
    //}
  }
  // cout << "GA evolve finished. best fit = \t"
  //	<< best[nVar] << endl;
  // cout << "best x =" << endl;
  // for (j = 0; j < nVar; ++j){
  //	cout << best[j] << endl;
  //}
  // cout << "nMutation = \t"
  //	<< nMutation << endl;
  // cout << "nCrossover = \t"
  //	<< nCrossover << endl;
}

template <typename Real>
void GA<Real>::tourney() {
  int i, rand, temp;
  for (i = 0; i < nPop; ++i) {
    tourlist[i] = i;
    tourlist[i + nPop] = i;
  }
  for (i = 0; i < nPop; ++i) {
    rand = Random::get<int>(0, nPop - 1);
    temp = tourlist[rand];
    tourlist[rand] = tourlist[i];
    tourlist[i] = temp;
  }
  for (i = nPop; i < 2 * nPop; ++i) {
    rand = Random::get<int>(nPop, 2 * nPop - 1);
    temp = tourlist[rand];
    tourlist[rand] = tourlist[i];
    tourlist[i] = temp;
  }
}

template <typename Real>
void GA<Real>::tourSelect(int k) {
  int i;
  int s1, s2;
  int kk = 2 * k;
  s1 = tourlist[kk];
  s2 = tourlist[kk + 1];
  if (nCons > 0) {
    if (consP[s1] > consP[s2])
      s1 = s2;
    else if (consP[s1] == 0 && consP[s2] == 0) {
      if (pop[s1][nVar] > pop[s2][nVar]) s1 = s2;
    }
  } else {
    if (pop[s1][nVar] > pop[s2][nVar]) s1 = s2;
  }
  for (i = 0; i < nVar; ++i)
    newpop[k][i] = pop[s1][i];

  s1 = tourlist[kk + 2];
  s2 = tourlist[kk + 3];
  if (nCons > 0) {
    if (consP[s1] > consP[s2])
      s1 = s2;
    else if (consP[s1] == 0 && consP[s2] == 0) {
      if (pop[s1][nVar] > pop[s2][nVar]) s1 = s2;
    }
  } else {
    if (pop[s1][nVar] > pop[s2][nVar]) s1 = s2;
  }
  for (i = 0; i < nVar; ++i)
    newpop[k + 1][i] = pop[s1][i];
}

template <typename Real>
void GA<Real>::SBXover(int k) {
  // if (Random::get<Real>(0, 1) > pCrossover)
  //	return;
  int  i;
  Real alpha, beta, mid, dif;
  for (i = 0; i < nVar; ++i) {
    if (Random::get<Real>(0, 1) < pCrossover) {
      if (newpop[k][i] > newpop[k + 1][i]) swap(newpop[k][i], newpop[k + 1][i]);
      mid = (newpop[k + 1][i] + newpop[k][i]) / 2;
      dif = newpop[k + 1][i] - newpop[k][i];
      if (dif < 1e-6) {
        dif = 1e-6;
      }

      beta =
          1 +
          2 / dif * min(newpop[k][i] - lower[i], upper[i] - newpop[k + 1][i]);
      alpha = 2 - pow(beta, -(nc + 1));
      if (alpha > 2 - 1e-6) alpha = 2 - 1e-6;
      alpha *= Random::get<Real>(0, 1);
      if (alpha <= 1)
        beta = pow(alpha, 1.0 / (nc + 1));
      else
        beta = 1.0 / pow(2 - alpha, 1.0 / (nc + 1));

      newpop[k][i] = mid - 0.5 * beta * dif;
      newpop[k + 1][i] = mid + 0.5 * beta * dif;
      if (isnan(newpop[k][i]) || isinf(newpop[k][i]))
        newpop[k][i] = Random::get<Real>(lower[i], upper[i]);
      if (newpop[k][i] > upper[i]) newpop[k][i] = upper[i];
      if (newpop[k][i] < lower[i]) newpop[k][i] = lower[i];
      if (isnan(newpop[k + 1][i]) || isinf(newpop[k + 1][i]))
        newpop[k][i] = Random::get<Real>(lower[i], upper[i]);
      if (newpop[k + 1][i] > upper[i]) newpop[k + 1][i] = upper[i];
      if (newpop[k + 1][i] < lower[i]) newpop[k + 1][i] = lower[i];
      ++nCrossover;
    }
  }
}

template <typename Real>
void GA<Real>::linerXover(int k) {
  // if (Random::get<Real>(0, 1) > pCrossover)
  //	return;
  int  i;
  Real dif;
  for (i = 0; i < nVar; ++i) {
    if (Random::get<Real>(0, 1) < pCrossover) {
      if (newpop[k][i] > newpop[k + 1][i]) swap(newpop[k][i], newpop[k + 1][i]);
      dif = newpop[k + 1][i] - newpop[k][i];
      newpop[k + 1][i] = newpop[k][i] + Random::get<Real>(0, 1) * dif;
      newpop[k][i] = newpop[k][i] + Random::get<Real>(0, 1) * dif;
      ++nCrossover;
    }
  }
}

template <typename Real>
void GA<Real>::PBMutation(int k) {
  // if (Random::get<Real>(0, 1) > pMutation)
  //	return;
  int i;
  for (i = 0; i < nVar; ++i) {
    if (Random::get<Real>(0, 1) < pMutation) {
      Real u = Random::get<Real>(0, 1);
      Real delta, deltamax;
      deltamax = upper[i] - lower[i];
      delta = min(newpop[k][i] - lower[i], upper[i] - newpop[k][i]) / deltamax;
      if (u > 0.5) {
        u = 1 - u;
        delta = 1 - pow((2 * u + (1 - 2 * u) * pow(1 - delta, nm + 1)),
                        1 / (nm + 1));
      } else {
        delta =
            pow((2 * u + (1 - 2 * u) * pow(1 - delta, nm + 1)), 1 / (nm + 1)) -
            1;
      }
      newpop[k][i] += delta * deltamax;
      if (isnan(newpop[k][i]) || isinf(newpop[k][i]))
        newpop[k][i] = Random::get<Real>(lower[i], upper[i]);
      if (newpop[k][i] > upper[i]) newpop[k][i] = upper[i];
      if (newpop[k][i] < lower[i]) newpop[k][i] = lower[i];
      ++nMutation;
    }
  }
}

template <typename Real>
void GA<Real>::rdMutation(int k) {
  // if (Random::get<Real>(0, 1) > pMutation)
  //	return;
  int i;
  for (i = 0; i < nVar; ++i) {
    if (Random::get<Real>(0, 1) < pMutation) {
      newpop[k][i] = Random::get<Real>(lower[i], upper[i]);
      ++nMutation;
    }
  }
}

template <typename Real>
void GA<Real>::initialpop() {
  int i, j;
  for (i = 0; i < nPop; ++i) {
    for (j = 0; j < nVar; ++j) {
      pop[i][j] = Random::get<Real>(lower[j], upper[j]);
    }
  }
  // cout << "GA initialized." << endl;
  best[nVar] = 1e9;
  for (i = 0; i < nVar; ++i)
    best[i] = pop[0][i];
  statistics();
}

template <typename Real>
void GA<Real>::statistics() {
  int  i, j;
  bool flag = 1;
  int  bestn = -1;
  for (i = 0; i < nPop; ++i) {
    fx(pop[i]);
    if (nCons > 0) {
      flag = 1;
      consP[i] = 0;
      for (j = 0; j < nCons; ++j) {
        if (pop[i][nVar + j + 1] < 0) {
          consP[i] -= pop[i][nVar + j + 1];
          flag = 0;
        }
      }
    }
    if (best[nVar] > pop[i][nVar] && flag) bestn = i;
  }
  if (bestn != -1) {
    for (i = 0; i < nCons + nVar + 1; ++i)
      best[i] = pop[bestn][i];
  }
}
