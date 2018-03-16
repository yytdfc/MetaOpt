#include "doe.h"
#include "common.h"
#include <fstream>
using namespace std;

template class Doe<double>;

template <typename Real>
Doe<Real>::Doe(int n_dim, int n_sample) : n_dim_(n_dim), n_sample_(n_sample) {
  vector<vector<Real>>(n_dim_, vector<Real>(n_sample_)).swap(sample_);
  ifadd = 0;
  hasSample = 0;
}

template <typename Real>
Doe<Real>::~Doe() {}

template <typename Real>
void Doe<Real>::printSample(string str) {
  if (!hasSample) {
    cerr << "doesn't have sample_!" << endl;
    return;
  }
  int      i, j;
  ofstream fout(str);
  fout.setf(ios::scientific);
  for (i = 0; i < n_sample_; ++i) {
    for (j = 0; j < n_dim_; ++j) {
      fout << sample_[j][i] << "\t";
    }
    fout << endl;
  }
  fout.close();
}

template <typename Real>
void Doe<Real>::printAdd(string str) {
  if (!ifadd) {
    cerr << "doesn't have sample_!" << endl;
    return;
  }
  int      i, j;
  ofstream fout(str);
  fout.setf(ios::scientific);
  for (i = 0; i < addNum; ++i) {
    for (j = 0; j < n_dim_; ++j) {
      fout << addSample[j][i] << "\t";
    }
    fout << endl;
  }
  fout.close();
}

template <typename Real>
void Doe<Real>::gen() {
  int  i, j;
  int* list = new int[n_sample_];
  int  rand;
  for (i = 0; i < n_dim_; ++i) {
    for (j = 0; j < n_sample_; ++j) {
      list[j] = j;
    }
    for (j = 0; j < n_sample_; ++j) {
      rand = Random::get<int>(j, n_sample_ - 1);
      sample_[i][j] = list[rand];
      sample_[i][j] = (sample_[i][j] + Random::get<Real>(0, 1)) / n_sample_;
      list[rand] = list[j];
    }
  }
  hasSample = 1;
  delete[] list;
}

template <typename Real>
void Doe<Real>::genMC() {
  int i, j;
  for (i = 0; i < n_dim_; ++i) {
    for (j = 0; j < n_sample_; ++j) {
      sample_[i][j] = Random::get<Real>(0, 1);
    }
  }
  hasSample = 1;
}

template <typename Real>
void Doe<Real>::readSample(string str) {
  int      i, j;
  ifstream fin(str);
  for (i = 0; i < n_sample_; ++i) {
    for (j = 0; j < n_dim_; ++j) {
      fin >> sample_[j][i];
    }
  }
  fin.close();
  hasSample = 1;
}

template <typename Real>
void Doe<Real>::addLHS(int num) {
  if (!hasSample) {
    cerr << "doesn't have sample_!" << endl;
    return;
  }
  int i, j, k;
  int position, rand;
  addNum = num;
  addSample = new Real*[n_dim_];
  for (i = 0; i < n_dim_; ++i) {
    addSample[i] = new Real[addNum];
  }
  int* map = new int[n_sample_ + addNum];
  int* list = new int[addNum];
  for (i = 0; i < n_dim_; ++i) {
    // initial map
    for (j = 0; j < n_sample_ + addNum; ++j) {
      map[j] = 0;
    }
    // mark exist points
    for (j = 0; j < n_sample_; ++j) {
      position = (int)(sample_[i][j] * (n_sample_ + addNum));
      for (k = 0; k < n_sample_ + addNum; ++k) {
        if (map[position + k] == 0 && (position + k) < (n_sample_ + addNum)) {
          map[position + k] = 1;
          break;
        } else if (map[position - k] == 0 && (position - k) >= 0) {
          map[position - k] = 1;
          break;
        }
      }
    }
    // get the blank points sequence
    for (j = 0, k = 0; j < n_sample_ + addNum; ++j) {
      if (map[j] == 0) {
        list[k] = j;
        ++k;
      }
    }
    // get LHS for add samples
    for (j = 0; j < addNum; ++j) {
      rand = Random::get<int>(j, addNum - 1);
      addSample[i][j] = list[rand];
      addSample[i][j] =
          (addSample[i][j] + Random::get<Real>(0, 1)) / (n_sample_ + addNum);
      list[rand] = list[j];
    }
  }
  ifadd = 1;
  delete[] map;
  delete[] list;
}
