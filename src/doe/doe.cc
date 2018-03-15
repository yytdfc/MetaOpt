#include "doe.h"
#include "common.h"
#include <fstream>
using namespace std;



template class Doe<double>;

template <typename Real>
Doe<Real>::Doe(int n, int num) {
  ndim = n;
  nSample = num;
  sample = new Real*[ndim];
  for (int i = 0; i < ndim; ++i) {
    sample[i] = new Real[nSample];
  }
  ifadd = 0;
  hasSample = 0;
}

template <typename Real>
Doe<Real>::~Doe() {
  for (int i = 0; i < ndim; ++i) {
    delete[] sample[i];
  }
  delete[] sample;
  if (ifadd) {
    for (int i = 0; i < ndim; ++i) {
      delete[] addSample[i];
    }
    delete[] addSample;
  }
}

template <typename Real>
void Doe<Real>::printSample(string str) {
  if (!hasSample) {
    cerr << "doesn't have sample!" << endl;
    return;
  }
  int      i, j;
  ofstream fout(str);
  fout.setf(ios::scientific);
  for (i = 0; i < nSample; ++i) {
    for (j = 0; j < ndim; ++j) {
      fout << sample[j][i] << "\t";
    }
    fout << endl;
  }
  fout.close();
}

template <typename Real>
void Doe<Real>::printAdd(string str) {
  if (!ifadd) {
    cerr << "doesn't have sample!" << endl;
    return;
  }
  int      i, j;
  ofstream fout(str);
  fout.setf(ios::scientific);
  for (i = 0; i < addNum; ++i) {
    for (j = 0; j < ndim; ++j) {
      fout << addSample[j][i] << "\t";
    }
    fout << endl;
  }
  fout.close();
}

template <typename Real>
void Doe<Real>::gen() {
  int  i, j;
  int* list = new int[nSample];
  int  rand;
  for (i = 0; i < ndim; ++i) {
    for (j = 0; j < nSample; ++j) {
      list[j] = j;
    }
    for (j = 0; j < nSample; ++j) {
      rand = Random::get<int>(j, nSample - 1);
      sample[i][j] = list[rand];
      sample[i][j] = (sample[i][j] + Random::get<Real>(0, 1)) / nSample;
      list[rand] = list[j];
    }
  }
  hasSample = 1;
  delete[] list;
}

template <typename Real>
void Doe<Real>::genMC() {
  int i, j;
  for (i = 0; i < ndim; ++i) {
    for (j = 0; j < nSample; ++j) {
      sample[i][j] = Random::get<Real>(0, 1);
    }
  }
  hasSample = 1;
}

template <typename Real>
void Doe<Real>::readSample(string str) {
  int      i, j;
  ifstream fin(str);
  for (i = 0; i < nSample; ++i) {
    for (j = 0; j < ndim; ++j) {
      fin >> sample[j][i];
    }
  }
  fin.close();
  hasSample = 1;
}

template <typename Real>
void Doe<Real>::addLHS(int num) {
  if (!hasSample) {
    cerr << "doesn't have sample!" << endl;
    return;
  }
  int i, j, k;
  int position, rand;
  addNum = num;
  addSample = new Real*[ndim];
  for (i = 0; i < ndim; ++i) {
    addSample[i] = new Real[addNum];
  }
  int* map = new int[nSample + addNum];
  int* list = new int[addNum];
  for (i = 0; i < ndim; ++i) {
    // initial map
    for (j = 0; j < nSample + addNum; ++j) {
      map[j] = 0;
    }
    // mark exist points
    for (j = 0; j < nSample; ++j) {
      position = (int)(sample[i][j] * (nSample + addNum));
      for (k = 0; k < nSample + addNum; ++k) {
        if (map[position + k] == 0 && (position + k) < (nSample + addNum)) {
          map[position + k] = 1;
          break;
        } else if (map[position - k] == 0 && (position - k) >= 0) {
          map[position - k] = 1;
          break;
        }
      }
    }
    // get the blank points sequence
    for (j = 0, k = 0; j < nSample + addNum; ++j) {
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
          (addSample[i][j] + Random::get<Real>(0, 1)) / (nSample + addNum);
      list[rand] = list[j];
    }
  }
  ifadd = 1;
  delete[] map;
  delete[] list;
}
