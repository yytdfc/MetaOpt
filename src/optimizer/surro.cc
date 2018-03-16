#include "surro.h"
#include "common.h"
#include <fstream>
#include <sstream>

using namespace std;

template class Surro<double>;

template <typename Real>
void Surro<Real>::setFx(Real f(Real* x)) {
  fx = f;
  LOG(INFO) << "Opt Function setted." << endl;
}
template <typename Real>
void Surro<Real>::setGKrigFx(Real f(Real* x)) {
  krig.setFx(f);
}
template <typename Real>
int Surro<Real>::readInput(string infile) {
  ifstream fin(infile);
  string   line;
  if (!fin) {
    cerr << "Input file doesn't exist!" << endl;
    return 0;
  }
  int i;
  // input Parameters for Opt
  getline(fin, line);
  getline(fin, line);
  fin >> nProb >> isRestart >> restartFrom;
  getline(fin, line);
  getline(fin, line);
  fin >> n_dim_ >> nCons >> nSample >> addN;
  getline(fin, line);
  getline(fin, line);
  allocate();
  dealAdd();
  for (i = 0; i < n_dim_; ++i)
    fin >> lower[i];
  getline(fin, line);
  for (i = 0; i < n_dim_; ++i) {
    fin >> upper[i];
    if (lower[i] > upper[i]) {
      swap(lower[i], upper[i]);  //	? unknow condition
    }
  }
  getline(fin, line);
  // input Parameters for DoE
  getline(fin, line);
  getline(fin, line);
  fin >> nInit >> nDoE;
  getline(fin, line);
  // input Parameters for Kriging model
  getline(fin, line);
  getline(fin, line);
  fin >> krig_corr >> krig_const_theta >> krig_porder;
  getline(fin, line);
  getline(fin, line);
  fin >> krig_dcmp >> krig_ParaOpt >> krig_EIcons;
  getline(fin, line);
  // input Parameters for GA
  getline(fin, line);
  getline(fin, line);
  fin >> ga_Pops >> ga_Gens >> ga_Pcr >> ga_Pmu;
  fin.close();
  krig_norm = 1;
  krig_regular = 1;
  krig_out_points = 2;
  krig.EIcons = krig_EIcons;
  Real* GAup = new Real[n_dim_];
  Real* GAlow = new Real[n_dim_];
  for (i = 0; i < n_dim_; ++i) {
    GAup[i] = 1;
    GAlow[i] = 0;
  }
  ga.GAinit(n_dim_, ga_Pops, ga_Gens, nCons, ga_Pcr, ga_Pmu, GAup, GAlow);
  nNow = 0;
  delete[] GAup;
  delete[] GAlow;
  return nProb;
}

template <typename Real>
void Surro<Real>::dealAdd() {
  if (addN > 10000) {
    addMethod[0] = addN / 10000;
    addN = addN % 10000;
    addMethod[1] = addN / 1000;
    addN = addN % 1000;
    addMethod[2] = addN / 100;
    addN = addN % 100;
    addMethod[3] = addN / 10;
    addMethod[4] = addN % 10;
    addN = 5;
  } else if (addN > 1000) {
    addMethod[0] = addN / 1000;
    addN = addN % 1000;
    addMethod[1] = addN / 100;
    addN = addN % 100;
    addMethod[2] = addN / 10;
    addMethod[3] = addN % 10;
    addN = 4;
  } else if (addN > 100) {
    addMethod[0] = addN / 100;
    addN = addN % 100;
    addMethod[1] = addN / 10;
    addMethod[2] = addN % 10;
    addN = 3;
  } else if (addN > 10) {
    addMethod[0] = addN / 10;
    addMethod[1] = addN % 10;
    addN = 2;
  } else {
    addMethod[0] = addN;
    addN = 1;
  }
  for (int i = 0; i < addN; ++i) {
    if (addMethod[i] > 5 || addMethod[i] == 0) {
      cerr << "WRONG add method!" << endl;
      exit(0);
    }
  }
}

template <typename Real>
void Surro<Real>::backupResult() {
  stringstream ss;
  string       filename = "result.dat";
  string       cmd;
  ifstream     fin(filename);
  if (fin.is_open()) {
    fin.close();
    for (int i = 1; true; i++) {
      ss << "result_back";
      ss << i;
      ss << ".dat";
      ss >> filename;
      LOG(INFO) << filename << endl;
      ss.clear();
      fin.open(filename);
      if (fin.is_open()) {
        fin.close();
        continue;
      } else {
        break;
      }
    }
    fin.close();
    cmd = "cp result.dat " + filename;
    LOG(INFO) << cmd << endl;
    system(cmd.c_str());
  }
  fin.close();
}

template <typename Real>
void Surro<Real>::dealRestart() {
  int      i, j;
  ifstream fin(restartFrom);
  if (fin.is_open()) {
    int    num, check = 1;
    string s;
    getline(fin, s);
    getline(fin, s);
    while (fin >> num) {
      LOG(INFO) << num << endl;
      if (check != num) {
        cerr << "Error in restarting!" << endl;
        exit(0);
      }
      for (i = 0; i < vecLengh; ++i) {
        fin >> sample[num - 1][i];
      }
      fin >> convergence[num - 1];
      ++check;
    }
    nNow = num;
    LOG(INFO) << "Restarting from " << nNow << " samples." << endl;
  } else {
    cerr << "Can't open " << restartFrom << endl;
  }
  fin.close();
  best[n_dim_] = 1e9;
  ofstream fout("result.dat");
  if (fout) {
    fout << "VARIABLES=num,x1";
    for (i = 1; i < n_dim_; ++i)
      fout << ",x" << i + 1;
    fout << ",y";
    for (i = 0; i < nCons; ++i)
      fout << ",con" << i + 1;
    fout << ",best" << endl;
    fout << "ZONE T=\"opt result\"" << endl;
    for (i = 0; i < nNow; ++i) {
      getBest(i);
      fout.setf(ios::scientific);
      fout << i + 1 << "\t";
      for (j = 0; j < vecLengh; ++j) {
        fout << sample[i][j] << "\t";
      }
      fout << convergence[i] << "\t";
      fout << endl;
    }
  } else {
    cerr << "Can't write result!" << endl;
  }
  fout.close();
}

template <typename Real>
void Surro<Real>::allocate() {
  int i;
  addMethod = new int[5];
  vecLengh = n_dim_ + nCons + 1;
  sample = new Real*[nSample];
  sample[0] = new Real[nSample * vecLengh];
  for (i = 1; i < nSample; ++i) {
    sample[i] = sample[i - 1] + vecLengh;
  }
  best = new Real[vecLengh];
  convergence = new Real[nSample];
  // convergence = { 0 };
  if (nCons > 0) {
    cons = new Real*[nSample];
    for (i = 0; i < nSample; ++i) {
      cons[i] = sample[i] + n_dim_ + 1;
    }
  }
  upper = new Real[n_dim_];
  lower = new Real[n_dim_];
  isInit = true;
}
template <typename Real>
void Surro<Real>::initialize(int   nv,
                             int   nS,
                             int   init,
                             int   nCon,
                             Real* up,
                             Real* low,
                             std::function<Real(Real* x)> f) {
  int i;
  n_dim_ = nv;
  nInit = init;
  nSample = nS;
  nCons = nCon;
  fx = f;
  for (i = 0; i < n_dim_; ++i) {
    upper[i] = up[i];
    lower[i] = low[i];
  }
  Real* GAup = new Real[n_dim_];
  Real* GAlow = new Real[n_dim_];
  for (i = 0; i < n_dim_; ++i) {
    GAup[i] = 1;
    GAlow[i] = 0;
  }
  ga.GAinit(n_dim_, ga_Pops, ga_Gens, nCons, ga_Pcr, ga_Pmu, GAup, GAlow);
  delete[] GAup;
  delete[] GAlow;
}

template <typename Real>
Surro<Real>::~Surro() {
  if (isInit) {
    delete[] addMethod;
    delete[] sample[0];
    delete[] sample;
    delete[] best;
    delete[] convergence;
    if (nCons > 0) {
      delete[] cons;
    }
    delete[] upper;
    delete[] lower;
  }
}
template <typename Real>
void Surro<Real>::getBest(const int& k) {
  bool flag = true;
  int  i;
  if (sample[k][n_dim_] < best[n_dim_]) {
    for (i = 0; i < nCons; ++i) {
      if (sample[k][n_dim_ + 1 + i] < 0) {
        flag = false;
      }
    }
    if (flag)
      for (i = 0; i < vecLengh; ++i) {
        best[i] = sample[k][i];
      }
  }
  convergence[k] = best[n_dim_];
}
template <typename Real>
void Surro<Real>::initSample() {
  int i, j;
  best[n_dim_] = 1e9;
  Doe<Real> initSamp(n_dim_, nInit);
  initSamp.gen();
  LOG(INFO) << "LHS samples generated." << endl;

  for (i = 0; i < nInit; ++i) {
    for (j = 0; j < n_dim_; ++j) {
      sample[i][j] = lower[j] + initSamp.sample_[j][i] * (upper[j] - lower[j]);
    }
    fx(sample[i]);
    getBest(i);
    ofstream fout("result.dat", ofstream::app);
    if (fout) {
      fout.setf(ios::scientific);
      fout << i + 1 << "\t";
      for (j = 0; j < vecLengh; ++j) {
        fout << sample[i][j] << "\t";
      }
      fout << convergence[i] << "\t";
      fout << endl;
    } else {
      cerr << "Can't write result!" << endl;
    }
    fout.close();
  }

  nNow = nInit;
  LOG(INFO) << "Initial samples calculated." << endl;
}

template <typename Real>
void Surro<Real>::add(const int nadd) {
  int i;
  krig.initialize(krig_corr, krig_const_theta, krig_porder, krig_norm,
                  krig_dcmp, krig_ParaOpt, krig_regular, n_dim_, nNow,
                  krig_out_points, 1 + nCons, sample, upper, lower);
  krig.GKtraining();
  if (nadd == 1) {
    krig.EI = best[n_dim_];
    ga.setFx(bind(&GKrig<Real>::GKpredictorEI, &krig, placeholders::_1));
    LOG(INFO) << "adding EI." << endl;
  } else if (nadd == 2) {
    ga.setFx(bind(&GKrig<Real>::GKpredictorMP, &krig, placeholders::_1));
    LOG(INFO) << "adding MP." << endl;
  } else if (nadd == 3) {
    ga.setFx(bind(&GKrig<Real>::predictorME, &krig, placeholders::_1));
    LOG(INFO) << "adding ME." << endl;
  } else if (nadd == 4) {
    ga.setFx(bind(&GKrig<Real>::predictorPI, &krig, placeholders::_1));
    LOG(INFO) << "adding PI." << endl;
  } else if (nadd == 5) {
    ga.setFx(bind(&GKrig<Real>::predictorLCB, &krig, placeholders::_1));
    LOG(INFO) << "adding LCB." << endl;
  }
  ga.evolve();
  krig.destructor();
  for (i = 0; i < n_dim_; ++i) {
    sample[nNow][i] = lower[i] + ga.best[i] * (upper[i] - lower[i]);
  }
  fx(sample[nNow]);
  getBest(nNow);
  LOG(INFO) << nNow << " samples best = " << best[n_dim_] << endl;
  ofstream fout("result.dat", ofstream::app);
  if (fout) {
    fout.setf(ios::scientific);
    fout << nNow + 1 << "\t";
    for (i = 0; i < vecLengh; ++i) {
      fout << sample[nNow][i] << "\t";
    }
    fout << convergence[nNow] << "\t";
    fout << endl;
  } else {
    cerr << "Can't write result!" << endl;
  }
  fout.close();
  ++nNow;
}

template <typename Real>
void Surro<Real>::opt() {
  timer = clock();
  backupResult();
  int i, j;
  if (isRestart) dealRestart();
  if (nNow == 0) {
    ofstream fout("result.dat");
    if (fout) {
      fout << "VARIABLES=num,x1";
      for (i = 1; i < n_dim_; i++)
        fout << ",x" << i + 1;
      fout << ",y";
      for (i = 0; i < nCons; i++)
        fout << ",con" << i + 1;
      fout << ",best" << endl;
      fout << "ZONE T=\"opt result\"" << endl;
    } else {
      cerr << "Can't write result!" << endl;
    }
    fout.close();
    initSample();
  }
  while (1) {
    for (i = 0; i < addN; ++i) {
      if (nNow < nSample)
        add(addMethod[i]);
      else
        break;
    }
    if (nNow >= nSample) break;
  }
  LOG(INFO) << endl << "best x = " << endl;
  for (i = 0; i < n_dim_; ++i) {
    LOG(INFO) << "\t" << best[i] << endl;
  }
  LOG(INFO) << "using " << (Real)(clock() - timer) / (Real)CLOCKS_PER_SEC
            << " secs" << endl;
}
