#include "kriging.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <stddef.h>
#include <string>
#define TINYa 1.0e-20

double Kriging::predictorEI(double* x) {
  double s, dy;
  setPredict(0);
  dy = EI - predictor(x);
  s = MSE(x);
  x[n_dim] = -dy * NormSDist(dy / s) - s * Normal(dy / s);
  if (EIcons) {
    for (int i = 1; i < ny; i++) {
      setPredict(i);
      x[n_dim] = x[n_dim] * NormSDist(predictor(x) / MSE(x));
      x[n_dim + i] = 1;
    }
  } else {
    for (int i = 1; i < ny; i++) {
      setPredict(i);
      x[n_dim + i] = predictor(x);
    }
  }
  return x[n_dim];
}
double Kriging::predictorMP(double* x) {
  for (int i = 0; i < ny; i++) {
    setPredict(i);
    x[n_dim + i] = predictor(x);
  }
  return x[n_dim];
}
double Kriging::predictorME(double* x) {
  setPredict(0);
  x[n_dim] = -MSE(x);
  for (int i = 1; i < ny; i++) {
    x[n_dim + i] = 1;
  }
  return x[n_dim];
}
double Kriging::predictorPI(double* x) {
  double s, dy;
  setPredict(0);
  dy = EI - predictor(x);
  s = MSE(x);
  x[n_dim] = -NormSDist(dy / s);
  if (EIcons) {
    for (int i = 1; i < ny; i++) {
      setPredict(i);
      x[n_dim] = x[n_dim] * NormSDist(predictor(x) / MSE(x));
      x[n_dim + i] = 1;
    }
  } else {
    for (int i = 1; i < ny; i++) {
      setPredict(i);
      x[n_dim + i] = predictor(x);
    }
  }
  return x[n_dim];
}
double Kriging::predictorLCB(double* x) {
  setPredict(0);
  x[n_dim] = predictor(x) - 4 * MSE(x);
  for (int i = 1; i < ny; i++) {
    setPredict(i);
    x[n_dim + i] = predictor(x);
  }
  return x[n_dim];
}

void Kriging::initialize(int      ncorr,
                         int      nconst_theta,
                         int      nporder,
                         int      nnorm,
                         int      ndcmp,
                         int      nParaOpt,
                         int      nregular,
                         int      ndim,
                         int      npoints,
                         int      nout_points,
                         int      nny,
                         double** xxx,
                         double*  up,
                         double*  low) {
  int i, j, k;
  rho_flag = 0;
  cokrig = 0;
  corr = ncorr;
  const_theta = nconst_theta;
  porder = nporder;
  norm = nnorm;
  dcmp = ndcmp;
  ParaOpt = nParaOpt;
  regular = nregular;
  n_dim = ndim;
  points = npoints;
  out_points = nout_points;
  ny = nny;
  if (points == 0 || n_dim == 0 || ny == 0) exit(0);
  /*----------------------------------------------------------------------------
  | Robustness treatment of input parameters
  ----------------------------------------------------------------------------*/
  if (corr >= 10) {
    dcmp = 0;
    const_theta = 1;
    ParaOpt = 1;
  }
  if (corr == 0 || corr == 2 || const_theta == 1) {
    ParaOpt = 1;
  }

  allocate();

  for (i = 0; i < points; ++i) {
    for (j = 0; j < n_dim; ++j) {
      xx[i][j] = xxx[i][j];
    }
    for (j = 0; j < ny; ++j) {
      nyy[j][i] = xxx[i][n_dim + j];
    }
    flag[i] = 0;
    // here wait for GEK
    // fin >> flag[i];
  }
  for (i = points; i < allpoints; ++i)
    flag[i] = 0;

  //---------------------get the xbound-------------------------------

  for (i = 0; i < n_dim; ++i) {
    xbound[i][1] = up[i];
    xbound[i][0] = low[i];
  }
  /*-------------------------------------------------------------
  | normalization of sampling data
  -------------------------------------------------------------*/
  if (norm == 1) {
    for (i = 0; i < points; i++) {
      for (j = 0; j < n_dim; j++) {
        xx[i][j] = (xx[i][j] - xbound[j][0]) / (xbound[j][1] - xbound[j][0]);
      }
      if (flag[i] > 0) {
        k = flag[i] - 1;
        for (j = 0; j < ny; j++)
          nyy[j][i] *= xbound[k][1] - xbound[k][0];
      }
    }
  }
}

void Kriging::allocate() {
  int i, j;
  /*----------------------------------------------------------------------------
  | set parameter for kriging model
  ----------------------------------------------------------------------------*/
  /*
  int _points, int _n_dim,
  int _cokrig, int _corr, int _porder,
  int _const_theta, int _norm, int _dcmp,
  int _ParaOpt, int _regular, int _ny

  n_dim = _n_dim;
  ny = _ny;
  points = _points;

  cokrig = _cokrig;
  corr = _corr;
  porder = _porder;
  const_theta = _const_theta;

  norm = _norm;
  dcmp = _dcmp;
  ParaOpt = _ParaOpt;
  regular = _regular;*/

  h_dim = n_dim;
  if (corr == 0) h_dim += n_dim;
  if (regular == 2) h_dim += 1;

  /*----------------------------------------------------------------------------
  | rho is multiplicative scalling factor between low- and hi-fidelity models
  ----------------------------------------------------------------------------*/
  if (rho_flag != 1.0) rho_flag = 0.0;

  /*----------------------------------------------------------------------------
  | NULL pointer to ylf and yhf, they are only for hybrid bridge function
  ----------------------------------------------------------------------------*/
  ylf = NULL;
  yhf = NULL;

  /*---------------------------------------------------------------------------
  | mu is added to v for regulariation
  |  mu  =  (10.0 +  points)*2.22E-16;
  ----------------------------------------------------------------------------*/
  if (regular == 2 || regular == 1)
    mu = (1000.0 + points) * 2.22E-16;
  else
    mu = 0.0;

  ncal_mle = 0;

  rho = 0.0;
  init_F = 0;
  init_phi = 0;
  init_phibar = 0;

  /*---------------------------------------------------------------------------
  | _porder <= 0 for simple kriging or RBFs
  ----------------------------------------------------------------------------*/
  if (porder < 0) {
    np = 0;
  }
  /*---------------------------------------------------------------------------
  | _porder == 0 for ordinary kriging or RBFs
  ----------------------------------------------------------------------------*/
  else if (porder == 0) {
    if (cokrig == 0)
      np = porder + 1;
    else
      np = 2 * (porder + 1);
  }
  /*---------------------------------------------------------------------------
  | _porder == 1 for universal kriging or RBFs
  ----------------------------------------------------------------------------*/
  else if (porder == 1) {
    if (cokrig == 0)
      np = porder + n_dim;
    else
      np = 2 * (porder + n_dim);
  } else {
    printf("\n Currently no such regression functionitiy,_porder = %d \n",
           porder);
    porder = 0;
    if (cokrig == 0)
      np = porder + 1;
    else
      np = 2 * (porder + 1);
  }

  allpoints = points + np;

  xx = (double**)malloc((points) * sizeof(double*));
  for (i = 0; i < points; i++)
    xx[i] = (double*)malloc((n_dim) * sizeof(double));
  xbound = (double**)malloc((n_dim) * sizeof(double*));
  for (i = 0; i < n_dim; i++)
    xbound[i] = (double*)malloc((2) * sizeof(double));

  flag = (int*)malloc((allpoints) * sizeof(int));

  nyy = (double**)malloc((ny) * sizeof(double*));
  for (i = 0; i < ny; i++)
    nyy[i] = (double*)malloc((allpoints) * sizeof(double));

  nyvi = (double**)malloc((ny) * sizeof(double*));
  for (i = 0; i < ny; i++)
    nyvi[i] = (double*)malloc((allpoints) * sizeof(double));
  nsigma_sq = (double*)malloc((ny) * sizeof(double));
  yvi = nyvi[0];

  nv = (double***)malloc((ny) * sizeof(double**));
  for (i = 0; i < ny; i++) {
    nv[i] = (double**)malloc((allpoints) * sizeof(double*));
    for (j = 0; j < allpoints; j++)
      nv[i][j] = (double*)malloc((allpoints) * sizeof(double));
  }
  v = nv[0];

  rvi = (double*)malloc((allpoints) * sizeof(double));

  ndiag = (double**)malloc((ny) * sizeof(double*));
  for (i = 0; i < ny; i++)
    ndiag[i] = (double*)malloc((allpoints) * sizeof(double));
  diag = ndiag[0];

  nindx = (int**)malloc((ny) * sizeof(int*));
  for (i = 0; i < ny; i++)
    nindx[i] = (int*)malloc((allpoints) * sizeof(int));
  indx = nindx[0];

  ntheta = (double**)malloc((ny) * sizeof(double*));
  for (i = 0; i < ny; i++)
    ntheta[i] = (double*)malloc((n_dim) * sizeof(double));
  theta = ntheta[0];

  npk = (double**)malloc((ny) * sizeof(double));
  for (i = 0; i < ny; i++)
    npk[i] = (double*)malloc((n_dim) * sizeof(double));
  pk = npk[0];

  nF = (double***)malloc((ny) * sizeof(double**));
  for (i = 0; i < ny; i++) {
    nF[i] = (double**)malloc((allpoints) * sizeof(double*));
    for (j = 0; j < allpoints; j++)
      nF[i][j] = (double*)malloc((allpoints) * sizeof(double));
  }

  F = nF[0];

  phi = (double*)malloc(np * sizeof(double));
  phibar = (double*)malloc(np * sizeof(double));

  isreadXinput = 0;
}

void Kriging::readInput(string infile) {
  int    i, j, k, restart;
  string line;
  cokrig = 0;
  ifstream fin(infile);
  if (!fin.is_open()) {
    cout << infile << " doesn't exist!" << endl;
    exit(0);
  }
  getline(fin, line);
  fin >> restart >> corr >> const_theta >> porder;
  getline(fin, line);
  getline(fin, line);
  fin >> norm >> dcmp >> ParaOpt >> regular;
  getline(fin, line);
  getline(fin, line);
  fin >> n_dim >> points >> out_points >> ny;
  if (points == 0 || n_dim == 0 || ny == 0) exit(0);
  /*----------------------------------------------------------------------------
  | Robustness treatment of input parameters
  ----------------------------------------------------------------------------*/
  if (corr >= 10) {
    dcmp = 0;
    const_theta = 1;
    ParaOpt = 1;
  }
  if (corr == 0 || corr == 2 || const_theta == 1) {
    ParaOpt = 1;
  }

  allocate();

  getline(fin, line);
  getline(fin, line);
  for (i = 0; i < n_dim; i++) {
    fin >> xbound[i][0] >> xbound[i][1];
    if (xbound[i][0] > xbound[i][1]) swap(xbound[i][0], xbound[i][1]);
  }

  getline(fin, line);
  getline(fin, line);
  for (i = 0; i < points; i++) {
    for (j = 0; j < n_dim; j++) {
      fin >> xx[i][j];
    }
    for (j = 0; j < ny; j++)
      fin >> nyy[j][i];
    fin >> flag[i];
  }
  //---------------------get the xbound-------------------------------
  // for (i = 0; i < n_dim; i++){
  //	xbound[i][0] = xx[0][i];
  //	xbound[i][1] = xx[0][i];
  //	for (j = 1; j<points; j++){
  //		if (xbound[i][0]> xx[j][i])
  //			xbound[i][0] = xx[j][i];
  //		if (xbound[i][1] < xx[j][i])
  //			xbound[i][1] = xx[j][i];
  //	}
  //	if (xbound[i][1] == xbound[i][0]){
  //		xbound[i][1] += 1;
  //		xbound[i][0] -= 1;
  //	}
  //}

  /*-------------------------------------------------------------
  | normalization of sampling data
  -------------------------------------------------------------*/
  if (norm == 1) {
    for (i = 0; i < points; i++) {
      for (j = 0; j < n_dim; j++) {
        xx[i][j] = (xx[i][j] - xbound[j][0]) / (xbound[j][1] - xbound[j][0]);
      }
      if (flag[i] > 0) {
        k = flag[i] - 1;
        for (j = 0; j < ny; j++)
          nyy[j][i] *= xbound[k][1] - xbound[k][0];
      }
    }
  }
  fin.close();

  cout << restart << "\t" << corr << "\t" << const_theta << "\t" << porder
       << endl;
  cout << norm << "\t" << dcmp << "\t" << ParaOpt << "\t" << regular << endl;
  cout << n_dim << "\t" << points << "\t" << out_points << "\t" << ny << endl;
  for (i = 0; i < n_dim; i++) {
    cout << xbound[i][0] << "\t" << xbound[i][1] << endl;
  }
}

/*---------------------------------------------------------------------------
| Step 2 : Kriging/GEK model fitting based on sampling data
---------------------------------------------------------------------------*/
void Kriging::training() {
  if (ParaOpt == 1) {
    fitting();
  } else if (ParaOpt == 2) {
    double* x = (double*)malloc((n_dim + 2) * sizeof(double));
    mle_initialization(x);
    mle_optimization(x);
    free(x);
  } else {
    printf("No such choice for ParaOpt = %d\n ", ParaOpt);
  }
}

/*----------------------------------------------------------------------------
| Step 3: Implement interpolation by kriging/GEK predictor and
|         Output the interpolated value and its MSE as tecplot file
----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
|  read parameters and sampled data from input file
----------------------------------------------------------------------------*/
void Kriging::readXinput() {
  int   i, j;
  FILE* fpin;
  if ((fpin = fopen("Xinput.txt", "r")) == NULL) {
    printf("Failed to open file Xinput.txt \n");
    exit(0);
  }
  fscanf(fpin, "%d ", &(ninterp));
  printf("ninterp = %d points\n", ninterp);
  Xinput = (double**)malloc((ninterp) * sizeof(double*));
  for (i = 0; i < ninterp; i++)
    Xinput[i] = (double*)malloc((n_dim) * sizeof(double));

  for (i = 0; i < ninterp; i++) {
    for (j = 0; j < n_dim; j++) {
      fscanf(fpin, "%lf ", &(Xinput[i][j]));
    }
  }

  fclose(fpin);
  isreadXinput = 1;
}

void Kriging::predictionDB() {
  int i, j, k;
  readXinput();
  double* xinterp;
  xinterp = (double*)malloc((n_dim) * sizeof(double));
  FILE* fpout4;
  if ((fpout4 = fopen("Youtput.dat", "w")) == NULL) {
    printf("Failed to open file Youtput.dat \n");
    exit(0);
  }
  for (k = 0; k < ny; k++) {
    setPredict(k);
    for (i = 0; i < ninterp; i++) {
      for (j = 0; j < n_dim; j++) {
        if (norm == 1)
          xinterp[j] =
              (Xinput[i][j] - xbound[j][0]) / (xbound[j][1] - xbound[j][0]);
        else
          xinterp[j] = Xinput[i][j];
      }
      fprintf(fpout4, "%lf\n", predictor(xinterp));
    }
    // fsprintf(fpout4,"\n");
  }
  fclose(fpout4);
  free(xinterp);

  printf("\n Please check interpolated results from \"output_rsm.dat\" \n");
  printf("\n Please check interpolated results from \"Youtput.dat\" \n");
  printf("\n");
}

void Kriging::prediction() {
  int i, j, k;

  /*----------------------------------------------------------------------------
  |  read parameters and sampled data from input file
  ----------------------------------------------------------------------------*/
  readXinput();
  FILE* fpout4;
  if ((fpout4 = fopen("Youtput.dat", "w")) == NULL) {
    printf("Failed to open file Youtput.dat \n");
    exit(0);
  }
  double **Youtput, **Yrmse, *xinterp;
  Youtput = (double**)malloc((ninterp) * sizeof(double*));
  for (i = 0; i < ninterp; i++)
    Youtput[i] = (double*)malloc((ny) * sizeof(double));
  Yrmse = (double**)malloc((ninterp) * sizeof(double*));
  for (i = 0; i < ninterp; i++)
    Yrmse[i] = (double*)malloc((ny) * sizeof(double));
  xinterp = (double*)malloc((n_dim) * sizeof(double));
  for (k = 0; k < ny; k++) {
    setPredict(k);
    for (i = 0; i < ninterp; i++) {
      for (j = 0; j < n_dim; j++) {
        if (norm == 1)
          xinterp[j] =
              (Xinput[i][j] - xbound[j][0]) / (xbound[j][1] - xbound[j][0]);
        else
          xinterp[j] = Xinput[i][j];
      }
      Youtput[i][k] = predictor(xinterp);
      Yrmse[i][k] = MSE(xinterp);
    }
  }

  fprintf(fpout4, "VARIABLES=x1");
  for (i = 1; i < n_dim; i++)
    fprintf(fpout4, ",x%d", i + 1);
  for (i = 0; i < ny; i++)
    fprintf(fpout4, ",y%d", i + 1);
  fprintf(fpout4, "\nZONE T=\"Y_output\",I=%d\n", ninterp);

  for (i = 0; i < ninterp; i++) {
    for (j = 0; j < n_dim; j++)
      fprintf(fpout4, "%le\t", Xinput[i][j]);
    for (j = 0; j < ny; j++)
      fprintf(fpout4, "%le\t", Youtput[i][j]);
    fprintf(fpout4, "\n");
  }
  fprintf(fpout4, "\n\n");

  fprintf(fpout4, "VARIABLES=x1");
  for (i = 1; i < n_dim; i++)
    fprintf(fpout4, ",x%d", i + 1);
  for (i = 0; i < ny; i++)
    fprintf(fpout4, ",y%d", i + 1);
  fprintf(fpout4, "\nZONE T=\"Y_MSE\",I=%d\n", ninterp);

  for (i = 0; i < ninterp; i++) {
    for (j = 0; j < n_dim; j++)
      fprintf(fpout4, "%le\t", Xinput[i][j]);
    for (j = 0; j < ny; j++)
      fprintf(fpout4, "%le\t", Yrmse[i][j]);
    fprintf(fpout4, "\n");
  }

  fclose(fpout4);
  for (i = 0; i < ninterp; i++) {
    free(Youtput[i]);
    free(Yrmse[i]);
  }
  free(Youtput);
  free(Yrmse);
  free(xinterp);
}

void Kriging::outputRSM() {
  int i, j, k;
  //----------------write "output_rsm.dat"-------------------

  FILE* fpout;
  FILE* fpout2;
  if ((fpout = fopen("output_rsm.dat", "w")) == NULL) {
    printf(" Failed to open file output_rsm.dat\n");
    exit(0);
  }
  if ((fpout2 = fopen("output_rsm_mse.dat", "w")) == NULL) {
    printf(" Failed to open file output_rsm_mse.dat\n");
    exit(0);
  }
  if (n_dim == 1) {
    fprintf(fpout, "VARIABLES=x1");
    for (i = 0; i < ny; i++)
      fprintf(fpout, ",y%d", i + 1);
    fprintf(fpout, "\n");
    fprintf(fpout, "ZONE T = \"Interpolated\", I=%d\n", out_points);
    fprintf(fpout2, "VARIABLES=x1");
    for (i = 0; i < ny; i++)
      fprintf(fpout2, ",y%d", i + 1);
    fprintf(fpout2, "\n");
    fprintf(fpout2, "ZONE T = \"Interpolated_MSE\", I=%d\n", out_points);
  } else if (n_dim == 2) {
    fprintf(fpout, "VARIABLES=x1,x2");
    for (i = 0; i < ny; i++)
      fprintf(fpout, ",y%d", i + 1);
    fprintf(fpout, "\n");
    fprintf(fpout, "ZONE T = \"Interpolated\",I=%d,J=%d\n", out_points,
            out_points);
    fprintf(fpout2, "VARIABLES=x1,x2");
    for (i = 0; i < ny; i++)
      fprintf(fpout2, ",y%d", i + 1);
    fprintf(fpout2, "\n");
    fprintf(fpout2, "ZONE T = \"Interpolated_MSE\",I=%d,J=%d\n", out_points,
            out_points);
  } else {
    ; /** to be added... **/
  }
  double  xout, xout1, yout, ymse;
  double* dx;
  dx = (double*)malloc((n_dim) * sizeof(double));
  for (i = 0; i < n_dim; i++) {
    if (out_points == 1)
      dx[i] = 0.0;
    else
      dx[i] = (xbound[i][1] - xbound[i][0]) / (out_points - 1.0);
  }

  double* xstar;
  xstar = (double*)malloc((n_dim) * sizeof(double));
  if (n_dim == 1) {
    for (i = 0; i < out_points; i++) {
      xout = xbound[0][0] + dx[0] * i;
      /*-------------------------------------------------------------
      | normalization of sampling data
      -------------------------------------------------------------*/
      if (norm == 1)
        xstar[0] = (xout - xbound[0][0]) / (xbound[0][1] - xbound[0][0]);
      else
        xstar[0] = xout;

      fprintf(fpout, "%le", xout);
      fprintf(fpout2, "%le", xout);
      for (k = 0; k < ny; k++) {
        setPredict(k);
        yout = predictor(xstar);
        ymse = MSE(xstar);
        fprintf(fpout, "\t%le", yout);
        fprintf(fpout2, "\t%le", ymse);
      }
      fprintf(fpout, "\n");
      fprintf(fpout2, "\n");
    }
  } else if (n_dim == 2) {
    for (j = 0; j < out_points; j++)
      for (i = 0; i < out_points; i++) {
        xout = xbound[0][0] + dx[0] * i;
        xout1 = xbound[1][0] + dx[1] * j;
        /*-------------------------------------------------------------
        | normalization of sampling data
        -------------------------------------------------------------*/
        if (norm == 1) {
          xstar[0] = (xout - xbound[0][0]) / (xbound[0][1] - xbound[0][0]);
          xstar[1] = (xout1 - xbound[1][0]) / (xbound[1][1] - xbound[1][0]);
        } else {
          xstar[0] = xout;
          xstar[1] = xout1;
        }

        fprintf(fpout, "%le\t%le", xout, xout1);
        fprintf(fpout2, "%le\t%le", xout, xout1);
        for (k = 0; k < ny; k++) {
          setPredict(k);
          yout = predictor(xstar);
          ymse = MSE(xstar);
          fprintf(fpout, "\t%le", yout);
          fprintf(fpout2, "\t%le", ymse);
        }
        fprintf(fpout, "\n");
        fprintf(fpout2, "\n");
      }
  } else {
    ; /** to be added ... **/
  }
  fclose(fpout);
  fclose(fpout2);

  free(xstar);
  free(dx);
}

void Kriging::setPredict(int k) {
  yvi = nyvi[k];  // Loop for n dims of Y
  sigma_sq = nsigma_sq[k];
  if (const_theta && !isGK)
    ;
  else {
    if (isGK) F = nF[k];
    v = nv[k];
    pk = npk[k];
    theta = ntheta[k];
    indx = nindx[k];
    diag = ndiag[k];
  }
}

/*******************************************************************************
* Func.  : Call destructor to free memory of kriging model
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/

void Kriging::destructor() {
  int i, j;
  /*----------------------------------------------------------------------------
  | free memory for sampled data
  ----------------------------------------------------------------------------*/
  for (i = 0; i < points; i++)
    free(xx[i]);
  free(xx);

  for (i = 0; i < n_dim; i++)
    free(xbound[i]);
  free(xbound);

  free(flag);
  for (i = 0; i < ny; i++)
    free(nyy[i]);  // Free nyy...
  free(nyy);
  /*----------------------------------------------------------------------------
  | free memory for y^T*v^(-1),v and r(x)^T*v^(-1)
  ----------------------------------------------------------------------------*/
  for (i = 0; i < ny; i++)
    free(nyvi[i]);
  free(nyvi);
  free(nsigma_sq);

  for (i = 0; i < ny; i++)
    for (j = 0; j < allpoints; j++)
      free(nv[i][j]);
  for (i = 0; i < ny; i++)
    free(nv[i]);
  free(nv);

  free(rvi);

  for (i = 0; i < ny; i++)
    free(ndiag[i]);
  free(ndiag);

  for (i = 0; i < ny; i++)
    free(nindx[i]);
  free(nindx);

  /*----------------------------------------------------------------------------
  | free memory for model fitting
  ----------------------------------------------------------------------------*/
  for (i = 0; i < ny; i++)
    free(ntheta[i]);
  free(ntheta);
  for (i = 0; i < ny; i++)
    free(npk[i]);
  free(npk);

  for (i = 0; i < allpoints; i++)
    free(F[i]);
  free(F);
  free(phi);
  free(phibar);
  if (isreadXinput) {
    for (i = 0; i < ninterp; i++)
      free(Xinput[i]);
    free(Xinput);
  }
}

/*******************************************************************************
* Func. :  kriging model initialization and fitting
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
void Kriging::fitting() {
  int i, j, k;
  int IS;

  for (i = 0; i < allpoints; i++)
    yvi[i] = 0.0;

  for (i = 0; i < allpoints; i++)
    rvi[i] = 0.0;

  /*----------------------------------------------------------------------------
  | Initialize regression  matrix F
  ----------------------------------------------------------------------------*/
  if (init_F == 0) {
    /*----------------------------------------------------------------
    | constant regression
    ---------------------------------------------------------------*/
    for (i = 0; i < points; i++) {
      if (flag[i] == 0)
        F[i][0] = 1.0;
      else
        F[i][0] = 0.0;
      if (cokrig > 0) {
        if (flag[i] == -1)
          F[i][1] = 1.0;
        else
          F[i][1] = 0.0;
      }
    } /* for(i=0;i< points; i++) */

    /*----------------------------------------------------------------
    | linear regression
    ----------------------------------------------------------------*/
    if (porder == 1) {
      for (j = 1; j <= n_dim; j++) {
        for (i = 0; i < points; i++) {
          if (flag[i] == 0)
            F[i][j] = xx[i][j - 1];
          else if (j == flag[i])
            F[i][j] = 1.0;
          else
            F[i][j] = 0.0;

        } /* for(i=0;i< points; i++) */

      } /* for(j=1;j<= n_dim;j++)*/

      if (cokrig > 0) {
        for (i = 0; i < points; i++) {
          if (flag[i] == -1)
            F[i][n_dim + 1] = 1.0;
          else
            F[i][n_dim + 1] = 0.0;
        }
        for (j = n_dim + 2; j < np; j++) {
          for (i = 0; i < points; i++) {
            if (flag[i] == -1)
              F[i][j] = xx[i][j - n_dim - 2];
            else if (j - np - n_dim == abs(flag[i] - 2))
              F[i][j] = 1.0;
            else
              F[i][j] = 0.0;
          }
        }
      } /* if( cokrig>0) */

    } /* if( porder==1) */

    for (i = points; i < allpoints; i++)
      for (j = 0; j < np; j++)
        F[i][j] = 0.0;
    init_F = 1;

  } /* if( init_F) */

  /*----------------------------------------------------------------------------
  | Initialize regression vector phi
  | Here only for porder = 0;for porder = 1,it will update for each
  interpolation
  ----------------------------------------------------------------------------*/
  if (init_phi == 0 && porder == 0) {
    phi[0] = 1.0;

    for (i = 1; i < np; i++)
      phi[i] = 0.0;
    init_phi = 1;
  }

  /*----------------------------------------------------------------------------
  | Initialize regression vector phibar for grad predictor
  | Here only for porder = 0;for porder = 1,it will update for each
  interpolation
  ----------------------------------------------------------------------------*/
  if (init_phibar == 0 && porder == 0) {
    phibar[0] = 0.0;

    for (i = 1; i < np; i++)
      phibar[i] = 0.0;
    init_phibar = 1;
  }

  /*----------------------------------------------------------------------------
  | Hyper parameter tuning
  ----------------------------------------------------------------------------*/
  if (const_theta) {
    /*----------------------------------------------------------------
    | use constant hyper parameters
    ----------------------------------------------------------------*/
    theta_init();
    // printf("\n The constant theta is:\n");
    // for (i = 0; i < n_dim; i++){
    //	printf(" theta[%d] = %lf\n", i, theta[i]);
    //	if (corr == 0)
    //		printf(" pk[%d] = %lf\n", i, pk[i]);
    //}
    // if (regular == 2 || regular == 1)
    //	printf(" mu = %le\n", mu);
    // printf("\n");

    /*----------------------------------------------------------------
    | call kriging_predictor_vector(paras) to cal. y^T*v^(-1)
    | store y^T*v^(-1) for predcitor y^T*v^(-1)*r(x)
    ----------------------------------------------------------------*/
    if (isGK) {
      for (i = 0; i < ny; i++) {
        F = nF[i];
        v = nv[i];
        yy = nyy[i];
        yvi = nyvi[i];
        theta = ntheta[i];
        pk = npk[i];
        diag = ndiag[i];
        indx = nindx[i];
        if (predictor_vector() != 0) {
          printf(" ---> Something wrong in kriging_fitting()! \n");
          printf(" ---> EXIT \n");
          exit(9928);
        }
      }
    } else {
      if (predictor_vector() != 0) {
        printf(" ---> Something wrong in kriging_fitting()! \n");
        printf(" ---> EXIT \n");
        exit(9928);
      }
    }
  } /* if( const_theta) */
  else {
    for (i = 0; i < ny; i++) {
      if (isGK) {
        F = nF[i];
      }
      v = nv[i];
      yy = nyy[i];
      yvi = nyvi[i];
      theta = ntheta[i];
      pk = npk[i];
      diag = ndiag[i];
      indx = nindx[i];
      /*--------------------------------------------------------------------------
      | BEGIN: hyper parameters opt.: fitting of kriging model
      --------------------------------------------------------------------------*/

      /*----------------------------------------------------------------
      | Initialize hyperparameters and their lower and uppper bounds
      | lbd : lower bound
      | ubd : upper bound
      ----------------------------------------------------------------*/
      /*----------------------------------------------------------------
      | *x; design varialbe array for optimization
      ----------------------------------------------------------------*/
      double *lbd, *ubd;
      double *lbd0, *ubd0;
      double* x;
      lbd0 = (double*)malloc((h_dim) * sizeof(double));
      ubd0 = (double*)malloc((h_dim) * sizeof(double));
      lbd = (double*)malloc((h_dim) * sizeof(double));
      ubd = (double*)malloc((h_dim) * sizeof(double));
      x = (double*)malloc((h_dim) * sizeof(double));

      theta_init();

      for (k = 0; k < n_dim; k++) {
        lbd[k] = 1.0e-8 * theta[k];
        ubd[k] = 1.0e3 * theta[k];
        lbd0[k] = 1.0e-8 * theta[k];
        ubd0[k] = 1.0e3 * theta[k];
        if (corr == 0) {
          lbd[k + n_dim] = 1.5;
          ubd[k + n_dim] = 2.0;
        }
      }

      if (corr == 10 || corr == 11) {
        lbd[h_dim - 1] = 0.001;
        ubd[h_dim - 1] = 100;
      }
      if (regular == 2) {
        mu = 1.0e-6;
        lbd[h_dim - 1] = 2.22e-16;
        ubd[h_dim - 1] = 1.0;
      }

      /*----------------------------------------------------------------
      | Begin to performe theta optimization
      | IS : number of iteration steps in Hooke and Jeeves method
      ----------------------------------------------------------------*/
      if (n_dim >= 4)
        IS = 4;
      else {
        if (n_dim >= 2)
          IS = n_dim;
        else
          IS = 2;
      }

      for (k = 0; k < n_dim; k++) {
        x[k] = theta[k];
        if (corr == 0) x[k + n_dim] = pk[k];
      }
      if (regular == 2) x[h_dim - 1] = mu;

      /*----------------------------------------------------------------
      | Hooke and Jeeves method for hyper parameters opt.
      ----------------------------------------------------------------*/
      int    niters = 10 * n_dim;
      double TRradius = 100.0;
      double MLEU = 999.0;
      // printf("\n");
      for (int ii = 0; ii < niters; ii++) {
        HookeAndJeeves(x, h_dim, lbd, ubd, IS);

        // printf(" %d-th trust region iteration, Maximum Likelihood is %f\n",
        //	ii, -mle);
        for (k = 0; k < n_dim; k++) {
          lbd[k] = max(lbd0[k], theta[k] - TRradius);
          ubd[k] = min(ubd0[k], theta[k] + TRradius);
        }
        double tol = fabs(MLEU - mle) / (fabs(mle) + 1e-16);
        if ((tol <= 1e-2) && (tol > 1.0e-6)) {
          TRradius *= 0.618 / 5.0;
        } else if (tol <= 1e-6)
          break;
        MLEU = mle;
      }
      /*--------------------------------------------------------------------------
      | End to performe hyper parameters optimization
      --------------------------------------------------------------------------*/

      /*--------------------------------------------------------------------------
      | free local memory
      --------------------------------------------------------------------------*/
      free(x);
      free(lbd);
      free(ubd);
      free(lbd0);
      free(ubd0);

      /*--------------------------------------------------------------------------
      | Output optimized theta
      --------------------------------------------------------------------------*/
      // printf("\n The optimized theta is:\n");
      // for (k = 0; k < n_dim; k++)
      //{
      //	printf(" theta[%d] = %lf\n", k, theta[i]);
      //	if (corr == 0)
      //		printf(" pk[%d] = %lf\n", k, pk[i]);
      //}
      // if (regular == 2 || regular == 1)
      //	printf(" mu = %le\n", mu);
      // printf("\n");
      nsigma_sq[i] = sigma_sq;
    }
  } /* if( const_theta) */
}

/*----------------------------------------------------------------------------
| Output fitted data
----------------------------------------------------------------------------*/
void Kriging::writeRestart(string fileRoute) {
  int         i, j;
  const char* file = fileRoute.c_str();
  FILE*       fpout;
  if ((fpout = fopen(file, "wb")) == NULL) {
    printf(" Failed to open file restart.dat\n");
    exit(0);
  }
  /*--------------------------------------------------------
  | Optimized hyper parameters
  --------------------------------------------------------*/
  fwrite(nsigma_sq, sizeof(double), ny, fpout);
  if (const_theta)
    fwrite(ntheta[0], sizeof(double), n_dim, fpout);
  else
    for (i = 0; i < ny; i++)
      fwrite(ntheta[i], sizeof(double), n_dim, fpout);
  if (corr == 0) {
    if (const_theta)
      fwrite(npk[0], sizeof(double), n_dim, fpout);
    else
      for (i = 0; i < ny; i++)
        fwrite(npk[i], sizeof(double), n_dim, fpout);
  }

  /*-----------------------------------------------------------------
  | regression vector phi
  -----------------------------------------------------------------*/
  fwrite(phi, sizeof(double), np, fpout);

  /*-----------------------------------------------------------------
  | Pre calculated yvi = y^T*v^(-1)
  -----------------------------------------------------------------*/

  for (i = 0; i < ny; i++)
    fwrite(nyvi[i], sizeof(double), allpoints, fpout);

  /*-----------------------------------------------------------------
  |  Pre decompled correlation matrix v and ralated data index or diag
  -----------------------------------------------------------------*/

  if (const_theta)
    for (i = 0; i < allpoints; i++)
      fwrite(nv[0][i], sizeof(double), allpoints, fpout);
  else
    for (j = 0; j < ny; j++)
      for (i = 0; i < allpoints; i++)
        fwrite(nv[j][i], sizeof(double), allpoints, fpout);

  if (dcmp == 0) {
    if (const_theta)
      fwrite(nindx[0], sizeof(int), allpoints, fpout);
    else
      for (i = 0; i < ny; i++)
        fwrite(nindx[i], sizeof(int), allpoints, fpout);
  }

  else {
    if (const_theta)
      fwrite(ndiag[0], sizeof(double), allpoints, fpout);
    else
      for (i = 0; i < ny; i++)
        fwrite(ndiag[i], sizeof(double), allpoints, fpout);
  }
  fclose(fpout);
}

/*---------------------------------------------------------------------------
| Read the trained rsm model data
----------------------------------------------------------------------------*/

void Kriging::readRestart(string fileRoute) {
  int         i, j;
  const char* file = fileRoute.c_str();
  FILE*       fpin;
  if ((fpin = fopen(file, "rb")) == NULL) {
    printf(" Failed to open file restart.dat\n");
    exit(0);
  }
  /*--------------------------------------------------------
  | Optimized hyper parameters
  --------------------------------------------------------*/

  fread(nsigma_sq, sizeof(double), ny, fpin);  // Binary Read

  if (const_theta)
    fread(ntheta[0], sizeof(double), n_dim, fpin);
  else
    for (i = 0; i < ny; i++)
      fread(ntheta[i], sizeof(double), n_dim, fpin);

  if (corr == 0) {
    if (const_theta)
      fread(npk[0], sizeof(double), n_dim, fpin);
    else
      for (i = 0; i < ny; i++)
        fread(npk[i], sizeof(double), n_dim, fpin);
  }

  /*-----------------------------------------------------------------
  | regression vector phi
  -----------------------------------------------------------------*/
  fread(phi, sizeof(double), np, fpin);

  /*-----------------------------------------------------------------
  | Pre calculated yvi = y^T*v^(-1)
  -----------------------------------------------------------------*/

  for (i = 0; i < ny; i++)
    fread(nyvi[i], sizeof(double), allpoints, fpin);

  /*-----------------------------------------------------------------
  |  Pre decompled correlation matrix v and ralated data index or diag
  -----------------------------------------------------------------*/

  if (const_theta)
    for (i = 0; i < allpoints; i++)
      fread(nv[0][i], sizeof(double), allpoints, fpin);
  else
    for (j = 0; j < ny; j++)
      for (i = 0; i < allpoints; i++)
        fread(nv[j][i], sizeof(double), allpoints, fpin);

  if (dcmp == 0) {
    if (const_theta)
      fread(nindx[0], sizeof(int), allpoints, fpin);
    else
      for (i = 0; i < ny; i++)
        fread(nindx[i], sizeof(int), allpoints, fpin);
  }

  else {
    if (const_theta)
      fread(ndiag[0], sizeof(double), allpoints, fpin);
    else
      for (i = 0; i < ny; i++)
        fread(ndiag[i], sizeof(double), allpoints, fpin);
  }
  fclose(fpin);
}
double Kriging::predictor_unNorm(double* x) {
  double* xx = new double[n_dim];
  if (norm == 1) {
    for (int i = 0; i < n_dim; ++i) {
      xx[i] = (x[i] - xbound[i][0]) / (xbound[i][1] - xbound[i][0]);
    }
  }
  double y = predictor(xx);
  delete[] xx;
  return y;
}

double Kriging::predictor(double* vec_interp) {
  int     i = 0, k = 0;
  double  lastval = 0.0;
  double* vstar;

  /*----------------------------------------------------------------------------
  | Initialize regression vector phi
  ----------------------------------------------------------------------------*/
  if (init_phi == 0 && porder == 1) {
    phi[0] = 1.0;
    for (i = 1; i <= n_dim; i++) {
      phi[i] = vec_interp[i - 1];
    }
    if (cokrig > 0) {
      for (i = n_dim + 1; i < np; i++) {
        phi[i] = 0.0;
      }
    }
    init_phi = 0;
  }

  /*----------------------------------------------------------------------------
  | allocate local memory for correlation vector r(x)
  ----------------------------------------------------------------------------*/
  vstar = (double*)malloc((allpoints) * sizeof(double));

  /*----------------------------------------------------------------------------
  | compute correlation vector r(x)
  ----------------------------------------------------------------------------*/
  for (i = 0; i < points; i++) {
    if (flag[i] == 0 || flag[i] == -1) {
      vstar[i] = correlation_of_r(&vec_interp[0], &xx[i][0]);
    } else if (flag[i] > 0 || flag[i] < -1) {
      if (flag[i] > 0)
        k = flag[i] - 1;
      else
        k = abs(flag[i]) - 2;
      vstar[i] = cross_correlation_of_r(&xx[i][0], &vec_interp[0], k);
    }
  }

  /*----------------------------------------------------------------------------
  | for augmented correlation vector
  ----------------------------------------------------------------------------*/
  for (i = points; i < allpoints; i++)
    vstar[i] = phi[i - points];

  /*----------------------------------------------------------------------------
  | compute yvi^T*r(x) = y^T*v^(-1)r(x)
  ----------------------------------------------------------------------------*/
  lastval = 0.0;
  for (i = 0; i < allpoints; i++) {
    lastval += yvi[i] * vstar[i];
  }

  /*----------------------------------------------------------------------------
  | free local memory
  ----------------------------------------------------------------------------*/
  free(vstar);

  return (lastval);

} /** kriging_predictor() **/

/*******************************************************************************
* Func.  : kriging gradient predictor, use it alfter model fitting
*          predict the gradient of unkown function at a untried site
*          paras is kriging structure
*          vec_interp is untried site x to be predicted
*          k is the k-th partial derivative
* Author : Zhong-Hua.Han
* Date   : 30.07.2009
*******************************************************************************/

double Kriging::grad_predictor(double* vec_interp, int k) {
  int     i = 0, l = 0;
  double  lastval = 0.0;
  double* vstar;

  /*----------------------------------------------------------------------------
  | Initialize regression vector phibar
  ----------------------------------------------------------------------------*/
  if (init_phibar == 0 && porder == 1) {
    phibar[0] = 0.0;
    for (i = 1; i <= n_dim; i++) {
      if (k == i)
        phibar[i] = 1.0;
      else
        phibar[i] = vec_interp[i - 1];
    }

    if (cokrig > 0) {
      for (i = n_dim + 1; i < np; i++) {
        phibar[i] = 0.0;
      }
    }

    init_phibar = 0;
  }

  /*----------------------------------------------------------------------------
  | allocate local memory for correlation vector r(x)
  ----------------------------------------------------------------------------*/
  vstar = (double*)malloc((allpoints) * sizeof(double));

  /*----------------------------------------------------------------------------
  | compute correlation vector r(x)
  ----------------------------------------------------------------------------*/
  for (i = 0; i < points; i++) {
    if (flag[i] == 0 || flag[i] == -1) {
      vstar[i] = cross_correlation_of_r(&vec_interp[0], &xx[i][0], k - 1);
    } else if (flag[i] > 0 || flag[i] < -1) {
      if (flag[i] > 0)
        l = flag[i] - 1;
      else
        l = abs(flag[i]) - 2;

      if (k - 1 == l) {
        vstar[i] = cross_correlation_of_r2(&vec_interp[0], &xx[i][0], k - 1);
      } else {
        vstar[i] =
            cross_correlation_of_r22(&vec_interp[0], &xx[i][0], k - 1, l);
      }
    }
  }

  /*----------------------------------------------------------------------------
  | for augmented correlation vector
  ----------------------------------------------------------------------------*/
  for (i = points; i < allpoints; i++)
    vstar[i] = phi[i - points];

  /*----------------------------------------------------------------------------
  | compute yvi^T*r(x) = y^T*v^(-1)r(x)
  ----------------------------------------------------------------------------*/
  lastval = 0.0;
  for (i = 0; i < allpoints; i++)
    lastval += yvi[i] * vstar[i];

  /*----------------------------------------------------------------------------
  | free local memory
  ----------------------------------------------------------------------------*/
  free(vstar);

  return (lastval);

} /** kriging_grad_predictor() **/

/*******************************************************************************
* Func.  : Mean Squared Error (MSE) estimation
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
double Kriging::MSE(double* vec_interp) {
  int     i, k;
  double  lastval = 0.0;
  double* vstar;

  /*----------------------------------------------------------------------------
  | Initialize regression vector phi
  ----------------------------------------------------------------------------*/
  if (init_phi == 0 && porder == 1) {
    phi[0] = 1.0;
    for (i = 1; i <= n_dim; i++) {
      phi[i] = vec_interp[i - 1];
    }

    if (cokrig > 0) {
      for (i = n_dim + 1; i < np; i++) {
        phi[i] = 0.0;
      }
    }

    init_phi = 0;
  }

  /*----------------------------------------------------------------------------
  | Define arraies *vstar,
  ----------------------------------------------------------------------------*/
  vstar = (double*)malloc((allpoints) * sizeof(double));

  /*----------------------------------------------------------------------------
  | Initialization of *vstar and rvi= r(x)^T*v^(-1)
  ----------------------------------------------------------------------------*/
  for (i = 0; i < points; i++) {
    if (flag[i] == 0 || flag[i] == -1) {
      vstar[i] = correlation_of_r(&vec_interp[0], &xx[i][0]);
    } else if (flag[i] > 0 || flag[i] < -1) {
      if (flag[i] > 0)
        k = flag[i] - 1;
      else
        k = abs(flag[i]) - 2;
      vstar[i] = cross_correlation_of_r(&xx[i][0], &vec_interp[0], k);
    }
  }

  /*----------------------------------------------------------------------------
  | for augmented correlation vector
  ----------------------------------------------------------------------------*/
  for (i = points; i < allpoints; i++)
    vstar[i] = phi[i - points];

  for (i = 0; i < allpoints; i++) {
    rvi[i] = vstar[i];
  }

  /*----------------------------------------------------------------------------
  | Computer v^(-1)*y through decomposition of v
  | LU decomp. Cholesky decomp. and a new cholesky decomp.
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    lubksb(v, allpoints, indx, rvi);
  } else if (dcmp == 1) {
    cholbksb(v, allpoints, diag, rvi, rvi);
  } else if (dcmp == 2) {
    ncholbksb(v, allpoints, diag, rvi, rvi);
  }

  lastval = 0.0;
  for (i = 0; i < allpoints; i++) {
    lastval += rvi[i] * vstar[i];
  }

  /*----------------------------------------------------------------------------
  | calculate value of correlation when distance equals to "0"
  ----------------------------------------------------------------------------*/
  double variance_r0 = correlation_of_r(&vec_interp[0], &vec_interp[0]);
  double mse = sqrt(sigma_sq * fabs(variance_r0 - lastval));

  /*----------------------------------------------------------------------------
  | free local memory
  ----------------------------------------------------------------------------*/
  free(vstar);

  return (mse);

} /** kriging_MSE() **/

/*******************************************************************************
* Func.  : Calculate and store yvi = y^T*v^(-1)
as result of kriging model fitting
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
int Kriging::predictor_vector() {
  int    i = 0;
  int    j = 0;
  double d = 0.0;

  /*----------------------------------------------------------------------------
  | calculate correlation matrix
  ----------------------------------------------------------------------------*/
  correlation_matrix();

  /*----------------------------------------------------------------------------
  | decomposition and solution of a set of linear equations
  | LU decomp. Cholesky Decomp. and a new Cholesky Decomp.
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    ludcmp(v, allpoints, indx, &d);
  } else if (dcmp == 1) {
    choldcmp(v, allpoints, diag);
  } else if (dcmp == 2) {
    ncholdcmp(v, allpoints, diag);
    psdf = 1;
    for (i = 0; i < points; i++) {
      if (diag[i] <= 0) {
        printf("None positive definite! correction...\n");
        psdf = 0;
      }
      diag[i] = max(diag[i], mu);
    }
  }

  /*----------------------------------------------------------------------------
  | Estimate optimal value of gamma (sigma_1/sigma_2) for cokriging
  ----------------------------------------------------------------------------*/

  for (i = 0; i < ny; i++) {
    yy = nyy[i];
    yvi = nyvi[i];

    if (cokrig > 0) {
      gamma = cokriging_gamma();
      printf("\n optimal value of gamma = : %lf\n\n", gamma);
    }

    /*----------------------------------------------------------------------------
    | Estimate optimal value of rho for hybrid bridge function
    ----------------------------------------------------------------------------*/
    if (rho_flag == 1.0) {
      rho = bridgefunction_rho();
      printf("\n clocsed-form optimal value of rho= : %lf\n\n", rho);
    }

    /*----------------------------------------------------------------------------
    | Calculate process variance to be used for MSE esitmation
    ----------------------------------------------------------------------------*/
    sigma_sq = process_variance();

    nsigma_sq[i] = sigma_sq;
  }
  return (0);

} /** kriging_predictor_vector() **/

/*******************************************************************************
* Func.  : Compute correlation matrix of kriging model
*          Augumented by extra column(s) and line(s) for unbiasedness constrains
*          Both for kriging, GECK, cokriging and GECK
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
int Kriging::correlation_matrix() {
  int i, j;
  int k, l;

  for (i = 0; i < points; i++) {
    for (j = i; j < points; j++) {
      if (flag[i] >= 0)
        k = flag[i] - 1;
      else
        k = abs(flag[i]) - 2;
      if (flag[j] >= 0)
        l = flag[j] - 1;
      else
        l = abs(flag[j]) - 2;

      if ((flag[i] == 0 && flag[j] == 0) || (flag[i] == -1 && flag[j] == -1) ||
          (flag[i] == -1 && flag[j] == 0) || (flag[i] == 0 && flag[j] == -1)) {
        v[i][j] = v[j][i] = correlation_of_r(&xx[i][0], &xx[j][0]);
      } else if ((flag[i] > 0 && flag[j] == 0) ||
                 (flag[i] > 0 && flag[j] == -1) ||
                 (flag[i] < -1 && flag[j] == 0) ||
                 (flag[i] < -1 && flag[j] == -1)) {
        v[i][j] = v[j][i] = cross_correlation_of_r(&xx[i][0], &xx[j][0], k);
      } else if ((flag[i] == 0 && flag[j] > 0) ||
                 (flag[i] == 0 && flag[j] < -1) ||
                 (flag[i] == -1 && flag[j] > 0) ||
                 (flag[i] == -1 && flag[j] < -1)) {
        v[i][j] = v[j][i] = cross_correlation_of_r(&xx[j][0], &xx[i][0], l);
      } else if ((flag[i] > 0 && flag[j] > 0) ||
                 (flag[i] < -1 && flag[j] < -1) ||
                 (flag[i] > 0 && flag[j] < -1) ||
                 (flag[i] < -1 && flag[j] > 0)) {
        if (k == l) {
          v[i][j] = v[j][i] = cross_correlation_of_r2(&xx[i][0], &xx[j][0], k);
        } else
          v[i][j] = v[j][i] =
              cross_correlation_of_r22(&xx[i][0], &xx[j][0], k, l);
      } /* else if( flag[i]>0&& flag[j]>0) */

    } /* for(j=i; j< points; j++) */

  } /* for(i=0; i< points; i++) */

  /*----------------------------------------------------------------------------
  | Augmented column(s) and line(s) for correlation matrix
  ----------------------------------------------------------------------------*/
  for (i = 0; i < allpoints; i++) {
    for (j = points; j < allpoints; j++) {
      v[i][j] = v[j][i] = F[i][j - points];
    }
  }

  /*----------------------------------------------------------------------------
  | regularization (R+mu*diag(1...1))
  ----------------------------------------------------------------------------*/
  for (i = 0; i < points; i++)
    v[i][i] += mu;

  /*----------------------------------------------------------------------------
  | sreen output for correlation matrix
  ----------------------------------------------------------------------------*/

  return (0);

} /**  kriging_correlation_matrix() **/

/*******************************************************************************
* Func   : Estimate process variance
* Author : By Zhonghua Han
* Date   : 06.08.2008
*******************************************************************************/
double Kriging::process_variance() {
  int     i;
  double  lastval = 0.0;
  double* vstar;

  /*----------------------------------------------------------------------------
  | Define arraies *vstar
  ----------------------------------------------------------------------------*/
  vstar = (double*)malloc((allpoints) * sizeof(double));

  /*----------------------------------------------------------------------------
  | Initialization of *vstar
  ----------------------------------------------------------------------------*/
  for (i = 0; i < allpoints; i++) {
    if (flag[i] >= 0) {
      if (rho_flag == 0.0) {
        vstar[i] = yy[i];
      } else if (rho_flag == 1.0) {
        yy[i] = vstar[i] = yhf[i] - rho * ylf[i];
      }

      yvi[i] = vstar[i];

    } /* if( flag[i]>=0) */
    else {
      vstar[i] = gamma * yy[i];
      yvi[i] = vstar[i];
    }
  } /* if( flag[i]>=0) */

  for (i = points; i < allpoints; i++) {
    yvi[i] = vstar[i] = 0.0;
  }

  /*----------------------------------------------------------------------------
  | Computer v^(-1)*y through decomposition of v
  | LU decomp. Cholesky decomp. and new Cholesky decomp.
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    lubksb(v, allpoints, indx, yvi);
  } else if (dcmp == 1) {
    cholbksb(v, allpoints, diag, yvi, yvi);
  } else if (dcmp == 2) {
    ncholbksb(v, allpoints, diag, yvi, yvi);
  }

  lastval = 0.0;
  for (i = 0; i < allpoints; i++) {
    lastval += yvi[i] * vstar[i];
  }

  /*----------------------------------------------------------------------------
  | free local memory
  ----------------------------------------------------------------------------*/
  free(vstar);

  /*----------------------------------------------------------------------------
  | return
  ----------------------------------------------------------------------------*/
  return (fabs(lastval / points));

} /** kriging_process_variance() **/

/*******************************************************************************
* Func.  : Spatial correlation function of two points with distance of r
* Author : By Zhonghua Han
* Date   : 06.08.2008
*******************************************************************************/
double Kriging::correlation_of_r(double* x1, double* x2) {
  int    i;
  double kesi, _theta;
  double temp;
  double r, _pk;

  temp = 1.0;
  if (corr < 10) {
    for (i = 0; i < n_dim; i++) {
      _theta = theta[i];
      r = fabs(x1[i] - x2[i]);
      switch (corr) {
        case 0:
          _pk = pk[i];
          temp *= exp(-_theta * pow(r, _pk));
          break;
        case 1:
          temp *= exp(-_theta * r * r);
          break;
        case 2:
          kesi = _theta * r;
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 1.0 - 15.0 * kesi * kesi * (1.0 - 2.0 * kesi);
          else if (kesi > 0.2 && kesi < 1)
            temp *= 1.25 * (1.0 - kesi) * (1.0 - kesi) * (1.0 - kesi);
          else
            temp *= 0.0;
          break;
        case 3:
          kesi = _theta * r;
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *=
                1.0 + kesi * kesi * (-15.0 + kesi * (35.0 + kesi * (-24.375)));
          else if (kesi > 0.4 && kesi < 1)
            temp *=
                (5.0 / 3.0) +
                kesi * ((-20.0 / 3.0) +
                        kesi * (10.0 +
                                kesi * ((-20.0 / 3.0) + kesi * (5.0 / 3.0))));
          else
            temp *= 0.0;
          break;
        default:
          printf("\n Currently there is no such correlation funciton !!!\n");

      } /* switch(corr) */

    } /* for(i=0;i< n_dim;i++) */
  }   /* if(corr<=10) */

  else {
    r = 0.0;
    for (i = 0; i < n_dim; i++) {
      r += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }

    r = sqrt(r);

    if (corr == 10) {
      temp = 1.0 / sqrt(1.0 + r * r);
    } else if (corr == 11) {
      temp = sqrt(1.0 + r * r);
    } else if (corr == 12) {
      r *= 10.0;
      if (r <= 1E-12)
        temp = 0;
      else
        temp = r * r * log(r);
    } else if (corr == 13) {
      temp = pow(r, 1.5);
    }
  }

  return (temp);

} /** kriging_correlation_of_r() **/

/*******************************************************************************
* Func.  : cross correlation,
i.e. first-oder patial derivative of correlation function
* Author : By Zhonghua Han
* Date   : 29.06.2009
*******************************************************************************/

double Kriging::cross_correlation_of_r(double* x1, double* x2, int k) {
  int    i;
  double kesi, _theta;
  double temp;
  double r, xi, xj;

  temp = 1.0;
  for (i = 0; i < n_dim; i++) {
    _theta = theta[i];
    r = fabs(x1[i] - x2[i]);
    xi = x1[i];
    xj = x2[i];

    switch (corr) {
      case 1:
        if (i == k)
          temp *= -2.0 * _theta * (xi - xj) * exp(-_theta * r * r);
        else
          temp *= exp(-_theta * r * r);
        break;
      case 2:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 30.0 * _theta * (3.0 * kesi - 1.0) * kesi * sign(xi - xj);
          else if (kesi > 0.2 && kesi < 1)
            temp *= -3.75 * _theta * (kesi - 1) * (kesi - 1) * sign(xi - xj);
          else
            temp *= 0.0;
        } else {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 1.0 - kesi * kesi * (15.0 - 30.0 * kesi);
          else if (kesi > 0.2 && kesi < 1)
            temp *= 1.25 * (1.0 - kesi) * (1.0 - kesi) * (1.0 - kesi);
          else
            temp *= 0.0;
        }
        break;
      case 3:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *= _theta *
                    (kesi * (-30.0 + kesi * (105.0 + kesi * (-97.5)))) *
                    sign(xi - xj);
          else if (kesi > 0.4 && kesi < 1)
            temp *= _theta *
                    ((-20.0 / 3.0) +
                     kesi * (20.0 + kesi * (-20.0 + kesi * (20.0 / 3.0)))) *
                    sign(xi - xj);
          else
            temp *= 0.0;
        } else {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *=
                1.0 + kesi * kesi * (-15.0 + kesi * (35.0 + kesi * (-24.375)));
          else if (kesi > 0.4 && kesi < 1)
            temp *=
                (5.0 / 3.0) +
                kesi * ((-20.0 / 3.0) +
                        kesi * (10.0 +
                                kesi * ((-20.0 / 3.0) + kesi * (5.0 / 3.0))));
          else
            temp *= 0.0;
        }
        break;
      default:
        printf("\n Currently there is no such correlation funciton !!!\n");

    } /* switch(corr) */

  } /* for(i=0;i< n_dim;i++) */

  return (temp);

} /** kriging_cross_correlation_of_r() **/

/*******************************************************************************
* Func.  : cross correlation,
*          i.e. second-order patial derivative of correlation function when k=l
* Author : By Zhonghua Han
* Date   : 29.06.2009
*******************************************************************************/
double Kriging::cross_correlation_of_r2(double* x1, double* x2, int k) {
  int    i;
  double kesi, _theta;
  double temp;

  double r, xi, xj;

  temp = 1.0;
  for (i = 0; i < n_dim; i++) {
    _theta = theta[i];
    r = fabs(x1[i] - x2[i]);
    xi = x1[i];
    xj = x2[i];

    switch (corr) {
      case 1:
        if (i == k)
          temp *= -2.0 * _theta * (2.0 * _theta * (xi - xj) * (xi - xj) - 1.0) *
                  exp(-_theta * r * r);
        else
          temp *= exp(-_theta * r * r);
        break;
      case 2:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= -30.0 * (6.0 * kesi - 1.0) * _theta * _theta;
          else if (kesi > 0.2 && kesi < 1)
            temp *= 7.5 * (kesi - 1.0) * _theta * _theta;
          else
            temp *= 0.0;
        } else {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 1.0 - kesi * kesi * (15.0 - 30.0 * kesi);
          else if (kesi > 0.2 && kesi < 1)
            temp *= 1.25 * (1.0 - kesi) * (1.0 - kesi) * (1.0 - kesi);
          else
            temp *= 0.0;
        }
        break;
      case 3:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *= (30.0 + kesi * (-210.0 + kesi * (292.5))) * _theta * _theta;
          else if (kesi > 0.4 && kesi < 1)
            temp *= (-20.0 + kesi * (40.0 + kesi * (-20.0))) * _theta * _theta;
          else
            temp *= 0.0;
        } else {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *=
                1.0 + kesi * kesi * (-15.0 + kesi * (35.0 + kesi * (-24.375)));
          else if (kesi > 0.4 && kesi < 1)
            temp *=
                (5.0 / 3.0) +
                kesi * ((-20.0 / 3.0) +
                        kesi * (10.0 +
                                kesi * ((-20.0 / 3.0) + kesi * (5.0 / 3.0))));
          else
            temp *= 0.0;
        }
        break;
      default:
        printf("\n Currently there is no such correlation funciton !!!\n");
    } /* switch(corr) */

  } /* for(i=0;i< n_dim;i++)*/

  return (temp);

} /** kriging_cross_correlation_of_r2() **/

/*******************************************************************************
* Func.  : cross correlation,
*          i.e. second-order patial derivative of correlation function when k!=l
* Author : By Zhonghua Han
* Date   : 29.06.2009
*******************************************************************************/

double Kriging::cross_correlation_of_r22(double* x1, double* x2, int k, int l) {
  int    i;
  double kesi, _theta;
  double temp;

  double r, xi, xj;

  temp = 1.0;
  for (i = 0; i < n_dim; i++) {
    _theta = theta[i];
    r = fabs(x1[i] - x2[i]);
    xi = x1[i];
    xj = x2[i];
    switch (corr) {
      case 1:
        if (i == k)
          temp *= -2.0 * _theta * (xi - xj) * exp(-_theta * r * r);
        else if (i == l)
          temp *= 2.0 * _theta * (xi - xj) * exp(-_theta * r * r);
        else
          temp *= exp(-_theta * r * r);
        break;
      case 2:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 30.0 * _theta * (3.0 * kesi - 1.0) * kesi * sign(xi - xj);
          else if (kesi > 0.2 && kesi < 1)
            temp *= -3.75 * _theta * (kesi - 1) * (kesi - 1) * sign(xi - xj);
          else
            temp *= 0.0;
        } else if (i == l) {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 30.0 * _theta * (3.0 * kesi - 1.0) * kesi * sign(xi - xj);
          else if (kesi > 0.2 && kesi < 1)
            temp *= -3.75 * _theta * (kesi - 1) * (kesi - 1) * sign(xi - xj);
          else
            temp *= 0.0;
          temp *= -1;
        } else {
          if (kesi >= 0.0 && kesi <= 0.2)
            temp *= 1.0 - kesi * kesi * (15.0 - 30.0 * kesi);
          else if (kesi > 0.2 && kesi < 1)
            temp *= 1.25 * (1.0 - kesi) * (1.0 - kesi) * (1.0 - kesi);
          else
            temp *= 0.0;
        }
        break;
      case 3:
        kesi = _theta * r;
        if (i == k) {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *= _theta *
                    (kesi * (-30.0 + kesi * (105.0 + kesi * (-97.5)))) *
                    sign(xi - xj);
          else if (kesi > 0.4 && kesi < 1)
            temp *= _theta *
                    ((-20.0 / 3.0) +
                     kesi * (20.0 + kesi * (-20.0 + kesi * (20.0 / 3.0)))) *
                    sign(xi - xj);
          else
            temp *= 0.0;
        } else if (i == l) {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *= _theta *
                    (kesi * (-30.0 + kesi * (105.0 + kesi * (-97.5)))) *
                    sign(xi - xj);
          else if (kesi > 0.4 && kesi < 1)
            temp *= _theta *
                    ((-20.0 / 3.0) +
                     kesi * (20.0 + kesi * (-20.0 + kesi * (20.0 / 3.0)))) *
                    sign(xi - xj);
          else
            temp *= 0.0;
          temp *= -1;
        } else {
          if (kesi >= 0.0 && kesi <= 0.4)
            temp *=
                1.0 + kesi * kesi * (-15.0 + kesi * (35.0 + kesi * (-24.375)));
          else if (kesi > 0.4 && kesi < 1)
            temp *=
                (5.0 / 3.0) +
                kesi * ((-20.0 / 3.0) +
                        kesi * (10.0 +
                                kesi * ((-20.0 / 3.0) + kesi * (5.0 / 3.0))));
          else
            temp *= 0.0;
        }
        break;
      default:
        printf("\n Currently there is no such correlation funciton !!!\n");

    } /* switch(corr) */

  } /* for(i=0;i< n_dim;i++) */

  return (temp);

} /** kriging_cross_correlation_of_r22() **/

/*******************************************************************************
* Func.  : sign of x
* Author : By Zhonghua Han
* Date   : 29.06.2009
*******************************************************************************/

double Kriging::sign(double x) {
  if (x > 0.0)
    return (1.0);
  else if (x < 0.0)
    return (-1.0);
  else if (x == 0.0)
    return (0.0);
  else {
    printf("\n error in function sign!\n");
    return (100);
  }
} /** sign() **/

/*******************************************************************************
* Determine function value of Mximum Likelihood Estimation (MLE)
* see eq. (2.5) "Apects Of The MATLAB Toolbox DACE", Lophaven,
* Compare to Chung/Alonso 2002 eq. (7), Laurenceau, Sageaut 2008 eq. (17), (18)
*
* Function value of sigma^2 is written to  sigma_sq
* MLE is written to  mle
* @ R. Zimmermann, Dec. 2008
********************************************************************************/

double Kriging::MLE(double* val) {
  int    i;
  double d;

  /*----------------------------------------------------------------------------
  | detcormatrix : determinate of correlation matrix
  ----------------------------------------------------------------------------*/

  double detcormatrix = 0.0;

  /*----------------------------------------------------------------------------
  | copy current parameter to kriging data structure
  | and reconsturct correlation matrix
  ----------------------------------------------------------------------------*/
  for (i = 0; i < n_dim; i++) {
    theta[i] = val[i];
    if (corr == 0) pk[i] = val[i + n_dim];
  }
  if (regular == 2) {
    mu = val[h_dim - 1];
  }

  correlation_matrix();
  ncal_mle += 1;

  /*----------------------------------------------------------------------------
  | decomposition of v and computation of its determinate
  | LU decomp. could be relpaced by Cholesky decomp.
  |            After call v is LU decomp of itself, diagonal of L is not stores
  | Cholesky dcomp. L is stored in lower triangle of v except the diagonal of L
  |            is stored in diag[]
  | Compute log of determinant by adding the log values of the diagonal of R
  | dcmp = 2 is strongly recommened in current implimentation
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    ludcmp(v, allpoints, indx, &d);
    detcormatrix = 0.0;
    for (i = 0; i < allpoints; i++)
      detcormatrix += log(fabs((double)(v[i][i])));
  } else if (dcmp == 1) {
    psdf = choldcmp(v, allpoints, diag);
    detcormatrix = 0.0;
    for (i = 0; i < points; i++)
      detcormatrix += log(diag[i] * diag[i]);
  } else if (dcmp == 2) {
    ncholdcmp(v, allpoints, diag);
    psdf = 1;
    for (i = 0; i < points; i++) {
      if (diag[i] <= 0) {
        printf("None positive definite! correction...\n");
        psdf = 0;
      }
      diag[i] = max(diag[i], mu);
    }
    detcormatrix = 0.0;
    for (i = 0; i < points; i++)
      detcormatrix += log(diag[i]);
  }

  /*----------------------------------------------------------------------------
  | Estimate optimal value of gamma (sigma_1/sigma_2) for cokriging
  ----------------------------------------------------------------------------*/
  if (cokrig > 0) {
    gamma = cokriging_gamma();
    printf("\n optimal value of gamma = : %lf\n\n", gamma);
  }

  /*----------------------------------------------------------------------------
  | Estimate optimal value of rho bybrid bridge function
  ----------------------------------------------------------------------------*/

  if (rho_flag == 1.0) {
    rho = bridgefunction_rho();
    printf("\n clsed form optimal value of rho= : %lf\n\n", rho);
  }

  sigma_sq = process_variance();

  mle = (double)points * log(sigma_sq) + detcormatrix;

  /*----------------------------------------------------------------------------
  | free local memory
  ----------------------------------------------------------------------------*/

  return (mle);

} /** kriging_MLE() **/

/*******************************************************************************
* Func.  : Estimate rho for hybrid bridge function
* Author : Zhong-Hua.Han
* Date   : 29.07.2009, 12:29 a.m.,  Wed.,July 29th, 2009
*******************************************************************************/
double Kriging::bridgefunction_rho() {
  int     i;
  double  sum1 = 0, sum2 = 0;
  double *y10, *y20, *y20_v_inv;
  double  rho;

  y10 = (double*)malloc((allpoints) * sizeof(double));
  y20 = (double*)malloc((allpoints) * sizeof(double));
  y20_v_inv = (double*)malloc((allpoints) * sizeof(double));

  for (i = 0; i < allpoints; i++) {
    y10[i] = yhf[i];
    y20[i] = ylf[i];
    y20_v_inv[i] = y20[i];
  }

  /*----------------------------------------------------------------------------
  | Computer v^(-1)*y20 through decomposition of v
  | LU decomp. Cholesky decomp. and a new cholesky decomp.
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    lubksb(v, allpoints, indx, y20_v_inv);
  } else if (dcmp == 1) {
    cholbksb(v, allpoints, diag, y20_v_inv, y20_v_inv);
  } else if (dcmp == 2) {
    ncholbksb(v, allpoints, diag, y20_v_inv, y20_v_inv);
  }

  for (i = 0; i < allpoints; i++) {
    sum1 += y20_v_inv[i] * y10[i];
    sum2 += y20_v_inv[i] * y20[i];
  }
  if (sum2 == 0)
    printf(" !!! something is wrong with gamma estimation in kriging_rho");

  rho = sum1 / sum2;
  // printf("sum1=%lf,sum2=%lf",sum1,sum2);

  free(y10);
  free(y20);
  free(y20_v_inv);

  return (rho);

} /** bridgefunction_rho() **/

/*******************************************************************************
* Func.  : Estimate gamma ( = sigma1/sigma2 )
*          opti. value is gamma = y20^T*v^(-1)*y10/(y20^T*v^(-1)*y20)
* Author : Zhong-Hua.Han
* Date   : 04.08.2010
*******************************************************************************/
double Kriging::cokriging_gamma() {
  int     i;
  double  sum1 = 0, sum2 = 0;
  double *y10, *y20, *y20_v_inv;
  double  gamma;

  y10 = (double*)malloc((allpoints) * sizeof(double));
  y20 = (double*)malloc((allpoints) * sizeof(double));
  y20_v_inv = (double*)malloc((allpoints) * sizeof(double));

  for (i = 0; i < allpoints; i++) {
    if (flag[i] >= 0) {
      y10[i] = -yy[i];
      y20[i] = 0;
    } else {
      y10[i] = 0;
      y20[i] = yy[i];
    }
    y20_v_inv[i] = y20[i];
  }

  /*----------------------------------------------------------------------------
  | Computer v^(-1)*y20 through decomposition of v
  | LU decomp. Cholesky decomp. and a new cholesky decomp.
  ----------------------------------------------------------------------------*/
  if (dcmp == 0) {
    lubksb(v, allpoints, indx, y20_v_inv);
  } else if (dcmp == 1) {
    cholbksb(v, allpoints, diag, y20_v_inv, y20_v_inv);
  } else if (dcmp == 2) {
    ncholbksb(v, allpoints, diag, y20_v_inv, y20_v_inv);
  }

  for (i = 0; i < allpoints; i++) {
    sum1 += y20_v_inv[i] * y10[i];
    sum2 += y20_v_inv[i] * y20[i];
  }
  if (sum2 == 0)
    printf(" !!! something is wrong with gamma estimation in kriging_gamma");

  gamma = sum1 / sum2;

  /*if(gamma<0)
  {
  printf("\n !!! gamma<0, improper value of gamma. Force it to zero \n");
  gamma = 0.0;
  }*/

  free(y10);
  free(y20);
  free(y20_v_inv);

  return (gamma);
} /** cokriging_gamma() **/
/*******************************************************************************
* Func.  : Theta initialization
* Author : R.Zimmermann
* Date   : Dec. 2008
*******************************************************************************/
void Kriging::theta_init() {
  int     i, j, k;
  double* max_k_dist;

  max_k_dist = (double*)malloc(n_dim * sizeof(double));

  /*----------------------------------------------------------------------------
  | compute initial distance for comparison
  | |x_k^0 - x_k^1|
  ----------------------------------------------------------------------------*/
  for (k = 0; k < n_dim; k++)
    max_k_dist[k] = fabs(xx[0][k] - xx[1][k]);

  /*----------------------------------------------------------------------------
  | Computation of maximal distance between sample points components
  ----------------------------------------------------------------------------*/
  for (k = 0; k < n_dim; k++) {
    for (i = 0; i < points; i++) {
      for (j = i + 1; j < points; j++) {
        if (max_k_dist[k] < fabs(xx[i][k] - xx[j][k]))
          max_k_dist[k] = fabs(xx[i][k] - xx[j][k]);
      }
    }
  }

  /*----------------------------------------------------------------------------
  | error treatment for max_k_dist[k]
  ----------------------------------------------------------------------------*/
  for (k = 0; k < n_dim; k++) {
    if (max_k_dist[k] == 0) {
      printf(" !!! inproper value for max_k_dist in %dth dimention\n", k);
      max_k_dist[k] = 1;
    }
  }
  /*----------------------------------------------------------------------------
  | Initial guess might be modified by other factors instead of log(5)
  | Initial theta values are stored in paras
  ----------------------------------------------------------------------------*/
  for (k = 0; k < n_dim; k++) {
    if (corr == 0) {
      theta[k] = 0.2 * log(5.0) / max_k_dist[k];
      pk[k] = 1.99;
    } else if (corr == 1) {
      theta[k] = 1.2 * log(5.0) / max_k_dist[k];
    } else {
      theta[k] = 0.1 / max_k_dist[k];
    }
  }
  free(max_k_dist);
}
void Kriging::ludcmp(double** a, int n, int* indx, double* d) {
  int     i, imax, j, k;
  double  big, dum, sum, temp;
  double* vv;

  vv = (double*)malloc(n * sizeof(double));
  *d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((temp = fabs(a[i - 1][j - 1])) > big) big = temp;
    if (big == 0.0) printf(" ---> Singular matrix in kriging_ludcmp");
    vv[i - 1] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = a[i - 1][j - 1];
      for (k = 1; k < i; k++)
        sum -= a[i - 1][k - 1] * a[k - 1][j - 1];
      a[i - 1][j - 1] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = a[i - 1][j - 1];
      for (k = 1; k < j; k++)
        sum -= a[i - 1][k - 1] * a[k - 1][j - 1];
      a[i - 1][j - 1] = sum;
      if ((dum = vv[i - 1] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    } /* for(i=j;i<=n;i++) */
    if (j != imax) {
      for (k = 1; k <= n; k++) {
        dum = a[imax - 1][k - 1];
        a[imax - 1][k - 1] = a[j - 1][k - 1];
        a[j - 1][k - 1] = dum;
      }
      *d = -(*d);
      vv[imax - 1] = vv[j - 1];
    }
    indx[j - 1] = imax;
    if (a[j - 1][j - 1] == 0.0) a[j - 1][j - 1] = TINYa;
    if (j != n) {
      dum = 1.0 / (a[j - 1][j - 1]);
      for (i = j + 1; i <= n; i++)
        a[i - 1][j - 1] *= dum;
    }
  } /*  for(j=1;j<=n;j++) */

  free(vv);

} /** kriging_ludcmp() **/

/*******************************************************************************
* Func.  : Soving the linear equations after performing LU decomposition
back substitution
* Author : By Zhonghua Han
* Date   : 29.06.2009
*******************************************************************************/
void Kriging::lubksb(double** a, int n, int* indx, double b[]) {
  int    i, ii = 0, ip, j;
  double sum;

  for (i = 1; i <= n; i++) {
    ip = indx[i - 1];
    sum = b[ip - 1];
    b[ip - 1] = b[i - 1];
    if (ii)
      for (j = ii; j <= i - 1; j++)
        sum -= a[i - 1][j - 1] * b[j - 1];
    else if (sum)
      ii = i;
    b[i - 1] = sum;
  }
  for (i = n; i >= 1; i--) {
    sum = b[i - 1];
    for (j = i + 1; j <= n; j++)
      sum -= a[i - 1][j - 1] * b[j - 1];
    b[i - 1] = sum / a[i - 1][i - 1];
  }
} /** kriging_lubksb() **/

/*******************************************************************************
* Func.  : Perform Cholesky decomposition of a square matrix
*          A = L*L^T
* Note   : only for symmetric, positive definite square matrix
* Input  : n is number of lines or columns of the matrix
* Output : diag[0...n-1] store the diagonal
*          Function itself will return a value to indicate whether the matrix
*          is positive definite (0/1)
* In/out : a[0...n-1][0..n-1], only the upper triangle of a[][] need be given
*          The Cholesky factor "L" is retruned in the lower triangle of a[][]
*          except the diagonal elements which are stored in diag
* Author : Zhonghua Han
* Date   : 13.07.2009
*******************************************************************************/
int Kriging::choldcmp(double** a, int n, double diag[]) {
  int    i, j, k;
  double sum;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      sum = a[i][j];
      for (k = i - 1; k >= 0; k--)
        sum -= a[i][k] * a[j][k];
      if (i == j) {
        if (sum <= 0.0) {
          printf(" choldcmp failed as the matrix is not positive definite!\n");
          return (0);
        } else {
          diag[i] = sqrt(sum);
        }
      } /* if(i == j) */
      else {
        a[j][i] = sum / diag[i];
      } /*if(i == j)*/

    } /* for(j=i;j<n,j++) */

  } /* for(i=0;i<n,i++) */

  return (1);
} /** kriging_choldcmp() **/

/*******************************************************************************
* Func.  : Sovle the set of n linear euqations A*x = b, where A is a positive
*          definite symmetric matrix.
* input  : a[0...n-1][0...n-1] and diag[0...n-1] are the output of function
*          kriging_choldc().Only lower triangle of a is  accessed.
*          b[0...n-1] is the right hand side vector
* Output : x[0...n-1] is the solution vector
* Author : Zhong-Hua Han
* Date   : 13.07.2009
*******************************************************************************/

void Kriging::cholbksb(double** a,
                       int      n,
                       double   diag[],
                       double   b[],
                       double   x[]) {
  int    i, k;
  double dum;
  /*----------------------------------------------------------------------------
  | Solve L*y=b,storing y in x
  ----------------------------------------------------------------------------*/
  for (i = 0; i < n; i++) {
    dum = b[i];
    for (k = i - 1; k >= 0; k--)
      dum -= a[i][k] * x[k];
    x[i] = dum / diag[i];
  }

  /*----------------------------------------------------------------------------
  | Solve L^T*x=y
  ----------------------------------------------------------------------------*/
  for (i = n - 1; i >= 0; i--) {
    dum = x[i];
    for (k = i + 1; k < n; k++)
      dum -= a[k][i] * x[k];
    x[i] = dum / diag[i];
  }

} /** kriging_cholbksb() **/

/*******************************************************************************
* Func.  : Perform new Cholesky decomposition of a square matrix
*          A = L*D*L^T
* Note   : for any symmetric square matrix
* Input  : n is number of lines or columns of the matrix
* Output : diag[0...n-1] store the diagonal
*          Function itself will return a value to indicate whether the matrix
*          is positive definite (0/1)
* In/out : a[0...n-1][0..n-1], only the upper triangle of a[][] need be given
*          The Cholesky factor "L" is retruned in the lower triangle of a[][]
*          except the diagonal elements which are stored in diag
* Author : Zhonghua Han
* Date   : 13.07.2009
*******************************************************************************/
void Kriging::ncholdcmp(double** a, int n, double diag[]) {
  int    i, j, k;
  double sum;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      sum = a[i][j];
      for (k = i - 1; k >= 0; k--)
        sum -= a[i][k] * a[j][k] * diag[k];
      if (i == j) {
        if (sum == 0.0) {
          printf(" the matrix is sigular\n");
          return;
        } else {
          diag[i] = sum;
        }
      } /* if(i == j) */
      else {
        a[j][i] = sum / diag[i];
      } /*if(i == j)*/

    } /* for(j=i;j<n,j++) */

  } /* for(i=0;i<n,i++) */

} /** kriging_ncholdcmp() **/

/*******************************************************************************
* Func.  : Sovle the set of n linear euqations A*x = b, where A is a
*          symmetric square matrix.
* input  : a[0...n-1][0...n-1] and diag[0...n-1] are the output of function
*          kriging_ncholdc().Only lower triangle of a is accessed.
*          b[0...n-1] is the right hand side vector
* Output : x[0...n-1] is the solution vector
* Author : Zhong-Hua Han
* Date   : 13.07.2009
*******************************************************************************/

void Kriging::ncholbksb(double** a,
                        int      n,
                        double   diag[],
                        double   b[],
                        double   x[]) {
  int    i, k;
  double dum;

  /*----------------------------------------------------------------------------
  | Solve L*z=b,storing z in x
  ----------------------------------------------------------------------------*/
  for (i = 0; i < n; i++) {
    dum = b[i];
    for (k = i - 1; k >= 0; k--)
      dum -= a[i][k] * x[k];
    x[i] = dum;
  }

  /*----------------------------------------------------------------------------
  | Solve D*y=b,storing y in x
  ----------------------------------------------------------------------------*/
  for (i = 0; i < n; i++) {
    x[i] = b[i] / diag[i];
  }

  /*----------------------------------------------------------------------------
  | Solve L^T*x=y
  ----------------------------------------------------------------------------*/
  for (i = n - 1; i >= 0; i--) {
    dum = x[i];
    for (k = i + 1; k < n; k++)
      dum -= a[k][i] * x[k];
    x[i] = dum;
  }

} /** kriging_cholbksb() **/
