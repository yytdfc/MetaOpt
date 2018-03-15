#include "gkrig.h"
#include "common.h"
#include <cstdio>
#include <fstream>
using namespace std;

template class GKrig<double>;

template <typename Real>
Real GKrig<Real>::GKpredictorEI(Real* x) {
  Real s, dy;
  GKsetPredict(0);
  dy = Kriging<Real>::EI - GKpredictor(x);
  s = Kriging<Real>::MSE(x);
  x[Kriging<Real>::n_dim] = -dy * Kriging<Real>::NormSDist(dy / s) - s * Kriging<Real>::Normal(dy / s);
  if (Kriging<Real>::EIcons) {
    for (int i = 1; i < Kriging<Real>::ny; i++) {
      GKsetPredict(i);
      x[Kriging<Real>::n_dim] =
          x[Kriging<Real>::n_dim] * Kriging<Real>::NormSDist(GKpredictor(x) / Kriging<Real>::MSE(x));
      x[Kriging<Real>::n_dim + i] = 1;
    }
  } else {
    for (int i = 1; i < Kriging<Real>::ny; i++) {
      GKsetPredict(i);
      x[Kriging<Real>::n_dim + i] = GKpredictor(x);
    }
  }
  return x[Kriging<Real>::n_dim];
}
template <typename Real>
Real GKrig<Real>::GKpredictorMP(Real* x) {
  for (int i = 0; i < Kriging<Real>::ny; i++) {
    GKsetPredict(i);
    x[Kriging<Real>::n_dim + i] = GKpredictor(x);
  }
  return x[Kriging<Real>::n_dim];
}

template <typename Real>
void GKrig<Real>::setFx(Real f(Real*)) {
  fx = f;
  Kriging<Real>::isGK = 1;
}
template <typename Real>
Real GKrig<Real>::fxIndex(Real* x, int index) {
  Real* xx = new Real[Kriging<Real>::n_dim + Kriging<Real>::ny];
  for (int i = 0; i < Kriging<Real>::n_dim; ++i) {
    xx[i] = x[i];
  }
  fx(xx);
  return (xx[Kriging<Real>::n_dim + index]);
}

template <typename Real>
void GKrig<Real>::GKinitialize(int    ncorr,
                               int    nconst_theta,
                               int    nporder,
                               int    nnorm,
                               int    ndcmp,
                               int    nParaOpt,
                               int    nregular,
                               int    ndim,
                               int    points,
                               int    nout_points,
                               int    nny,
                               Real** xx,
                               Real*  up,
                               Real*  low) {
  Kriging<Real>::initialize(ncorr, nconst_theta, nporder, nnorm, ndcmp, nParaOpt, nregular,
             ndim, points, nout_points, nny, xx, up, low);
}

template <typename Real>
void GKrig<Real>::GKsetPredict(int k) {
  Kriging<Real>::setPredict(k);
  if (Kriging<Real>::isGK) yIndex = k;
}

/*******************************************************************************
* Func. :  two-lelvels hierachy kriging model initialization and fitting
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
template <typename Real>
void GKrig<Real>::GKtraining() {
  if (Kriging<Real>::isGK) {
    int i, j;
    /*----------------------------------------------------------------------------
    | Initialize regression  matrix F
    ----------------------------------------------------------------------------*/
    Kriging<Real>::init_F = 1;
    Kriging<Real>::init_phi = 1;
    for (int k = 0; k < Kriging<Real>::ny; ++k) {
      Kriging<Real>::F = Kriging<Real>::nF[k];
      /*--------------------------------------------------------------------------
      | constant regression
      --------------------------------------------------------------------------*/
      for (i = 0; i < Kriging<Real>::points; i++) {
        if (Kriging<Real>::flag[i] == 0) {
          Kriging<Real>::F[i][0] = fxIndex(Kriging<Real>::xx[i], k);
        } else {
          cerr << "can't deal with gradient." << endl;
          Kriging<Real>::F[i][0] = 0;
        }
      } /*for(i=0;i< krighf->points; i++)*/

      /*----------------------------------------------------------------------------
      | linear regression
      ----------------------------------------------------------------------------*/
      if (Kriging<Real>::porder >= 1) {
        for (j = 1; j <= Kriging<Real>::n_dim; j++) {
          for (i = 0; i < Kriging<Real>::points; i++) {
            if (Kriging<Real>::flag[i] == 0) {
              Kriging<Real>::F[i][j] = Kriging<Real>::xx[i][j - 1] +
                                       fxIndex(Kriging<Real>::xx[i], k);
            } else if (j == Kriging<Real>::flag[i])
              Kriging<Real>::F[i][j] = 1.0;
            else
              Kriging<Real>::F[i][j] = 0.0;
          }
        }
      }

      for (i = Kriging<Real>::points; i < Kriging<Real>::allpoints; i++)
        for (j = 0; j < Kriging<Real>::np; j++)
          Kriging<Real>::F[i][j] = 0.0;
    }
  }
  Kriging<Real>::training();
}

/*******************************************************************************
* Func.  : two-level hierarchy kriging predictor, use it alfter model fitting
*          predict the respone at a untried site
*          paras is  kriging structure
*          vec_interp is untried site x to be predicted
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
template <typename Real>
Real GKrig<Real>::GKpredictor(Real* vec_interp) {
  if (Kriging<Real>::isGK) {
    int i;
    /*----------------------------------------------------------------------------
    | Initialize regression vector Kriging<Real>::phi
    ----------------------------------------------------------------------------*/
    Kriging<Real>::init_phi = 1;
    Kriging<Real>::phi[0] = fxIndex(vec_interp, yIndex);

    if (Kriging<Real>::porder == 1) {
      for (i = 1; i <= Kriging<Real>::n_dim; i++) {
        Kriging<Real>::phi[i] = fxIndex(vec_interp, yIndex) * vec_interp[i - 1];
      }
    }
  }
  Kriging<Real>::init_phi = 0;
  return (Kriging<Real>::predictor(vec_interp));
} /** krig2h_predictor() **/

/*******************************************************************************
* Func.  : Mean Squared Error (MSE) estimation for two-levels hierarchy krigings
* Author : Zhong-Hua.Han
* Date   : 17.07.2009
*******************************************************************************/
template <typename Real>
Real GKrig<Real>::GKMSE(Real* vec_interp) {
  if (Kriging<Real>::isGK) {
    int  i;
    Real rmse;
    Kriging<Real>::init_phi = 1;
    Kriging<Real>::phi[0] = fxIndex(vec_interp, yIndex);

    if (Kriging<Real>::porder == 1) {
      for (i = 1; i <= Kriging<Real>::n_dim; i++) {
        Kriging<Real>::phi[i] = vec_interp[i - 1];
      }
      for (i = Kriging<Real>::n_dim + 1; i < Kriging<Real>::np; i++) {
        Kriging<Real>::phi[i] = 0.0;
      }
    }
    Kriging<Real>::init_phi = 0;
  }
  return (Kriging<Real>::MSE(vec_interp));

} /** krig2h_MSE() **/

template <typename Real>
void GKrig<Real>::GKoutputRSM() {
  int i, j, k;
  //----------------write "output_rsm.dat"-------------------
  FILE* fpout;
  FILE* fpout2;
  if ((fpout = fopen("output_GK.dat", "w")) == NULL) {
    printf(" Failed to open file output_rsm.dat\n");
    exit(0);
  }
  if ((fpout2 = fopen("output_rsm_mse.dat", "w")) == NULL) {
    printf(" Failed to open file output_rsm_mse.dat\n");
    exit(0);
  }
  if (Kriging<Real>::n_dim == 1) {
    fprintf(fpout, "VARIABLES=x1");
    for (i = 0; i < Kriging<Real>::ny; i++)
      fprintf(fpout, ",y%d", i + 1);
    fprintf(fpout, "\n");
    fprintf(fpout, "ZONE T = \"Interpolated\", I=%d\n",
            Kriging<Real>::out_points);
    fprintf(fpout2, "VARIABLES=x1");
    for (i = 0; i < Kriging<Real>::ny; i++)
      fprintf(fpout2, ",y%d", i + 1);
    fprintf(fpout2, "\n");
    fprintf(fpout2, "ZONE T = \"Interpolated_MSE\", I=%d\n",
            Kriging<Real>::out_points);
  } else if (Kriging<Real>::n_dim == 2) {
    fprintf(fpout, "VARIABLES=x1,x2");
    for (i = 0; i < Kriging<Real>::ny; i++)
      fprintf(fpout, ",y%d", i + 1);
    fprintf(fpout, "\n");
    fprintf(fpout, "ZONE T = \"Interpolated\",I=%d,J=%d\n",
            Kriging<Real>::out_points, Kriging<Real>::out_points);
    fprintf(fpout2, "VARIABLES=x1,x2");
    for (i = 0; i < Kriging<Real>::ny; i++)
      fprintf(fpout2, ",y%d", i + 1);
    fprintf(fpout2, "\n");
    fprintf(fpout2, "ZONE T = \"Interpolated_MSE\",I=%d,J=%d\n",
            Kriging<Real>::out_points, Kriging<Real>::out_points);
  } else {
    ; /** to be added... **/
  }
  Real  xout, xout1, yout, ymse;
  Real* dx;
  dx = (Real*)malloc((Kriging<Real>::n_dim) * sizeof(Real));
  for (i = 0; i < Kriging<Real>::n_dim; i++) {
    if (Kriging<Real>::out_points == 1)
      dx[i] = 0.0;
    else
      dx[i] = (Kriging<Real>::xbound[i][1] - Kriging<Real>::xbound[i][0]) /
              (Kriging<Real>::out_points - 1.0);
  }

  Real* xstar;
  xstar = (Real*)malloc((Kriging<Real>::n_dim) * sizeof(Real));
  if (Kriging<Real>::n_dim == 1) {
    for (i = 0; i < Kriging<Real>::out_points; i++) {
      xout = Kriging<Real>::xbound[0][0] + dx[0] * i;
      /*-------------------------------------------------------------
      | normalization of sampling data
      -------------------------------------------------------------*/
      if (Kriging<Real>::norm == 1)
        xstar[0] = (xout - Kriging<Real>::xbound[0][0]) /
                   (Kriging<Real>::xbound[0][1] - Kriging<Real>::xbound[0][0]);
      else
        xstar[0] = xout;

      fprintf(fpout, "%le", xout);
      fprintf(fpout2, "%le", xout);
      for (k = 0; k < Kriging<Real>::ny; k++) {
        GKsetPredict(k);
        yout = GKpredictor(xstar);
        ymse = GKMSE(xstar);
        fprintf(fpout, "\t%le", yout);
        fprintf(fpout2, "\t%le", ymse);
      }
      fprintf(fpout, "\n");
      fprintf(fpout2, "\n");
    }
  } else if (Kriging<Real>::n_dim == 2) {
    for (j = 0; j < Kriging<Real>::out_points; j++)
      for (i = 0; i < Kriging<Real>::out_points; i++) {
        xout = Kriging<Real>::xbound[0][0] + dx[0] * i;
        xout1 = Kriging<Real>::xbound[1][0] + dx[1] * j;
        /*-------------------------------------------------------------
        | normalization of sampling data
        -------------------------------------------------------------*/
        if (Kriging<Real>::norm == 1) {
          xstar[0] =
              (xout - Kriging<Real>::xbound[0][0]) /
              (Kriging<Real>::xbound[0][1] - Kriging<Real>::xbound[0][0]);
          xstar[1] =
              (xout1 - Kriging<Real>::xbound[1][0]) /
              (Kriging<Real>::xbound[1][1] - Kriging<Real>::xbound[1][0]);
        } else {
          xstar[0] = xout;
          xstar[1] = xout1;
        }

        fprintf(fpout, "%le\t%le", xout, xout1);
        fprintf(fpout2, "%le\t%le", xout, xout1);
        for (k = 0; k < Kriging<Real>::ny; k++) {
          GKsetPredict(k);
          yout = GKpredictor(xstar);
          ymse = GKMSE(xstar);
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
template <typename Real>
void GKrig<Real>::GKprediction() {
  int i, j, k;

  /*----------------------------------------------------------------------------
  |  read parameters and sampled data from input file
  ----------------------------------------------------------------------------*/
  Kriging<Real>::readXinput();
  FILE* fpout4;
  if ((fpout4 = fopen("Youtput.dat", "w")) == NULL) {
    printf("Failed to open file Youtput.dat \n");
    exit(0);
  }
  Real **Youtput, **Yrmse, *xinterp;
  Youtput = (Real**)malloc((Kriging<Real>::ninterp) * sizeof(Real*));
  for (i = 0; i < Kriging<Real>::ninterp; i++)
    Youtput[i] = (Real*)malloc((Kriging<Real>::ny) * sizeof(Real));
  Yrmse = (Real**)malloc((Kriging<Real>::ninterp) * sizeof(Real*));
  for (i = 0; i < Kriging<Real>::ninterp; i++)
    Yrmse[i] = (Real*)malloc((Kriging<Real>::ny) * sizeof(Real));
  xinterp = (Real*)malloc((Kriging<Real>::n_dim) * sizeof(Real));
  for (k = 0; k < Kriging<Real>::ny; k++) {
    GKsetPredict(k);
    for (i = 0; i < Kriging<Real>::ninterp; i++) {
      for (j = 0; j < Kriging<Real>::n_dim; j++) {
        if (Kriging<Real>::norm == 1)
          xinterp[j] =
              (Kriging<Real>::Xinput[i][j] - Kriging<Real>::xbound[j][0]) /
              (Kriging<Real>::xbound[j][1] - Kriging<Real>::xbound[j][0]);
        else
          xinterp[j] = Kriging<Real>::Xinput[i][j];
      }
      Youtput[i][k] = GKpredictor(xinterp);
      Yrmse[i][k] = GKMSE(xinterp);
    }
  }

  fprintf(fpout4, "VARIABLES=x1");
  for (i = 1; i < Kriging<Real>::n_dim; i++)
    fprintf(fpout4, ",x%d", i + 1);
  for (i = 0; i < Kriging<Real>::ny; i++)
    fprintf(fpout4, ",y%d", i + 1);
  fprintf(fpout4, "\nZONE T=\"Y_output\",I=%d\n", Kriging<Real>::ninterp);

  for (i = 0; i < Kriging<Real>::ninterp; i++) {
    for (j = 0; j < Kriging<Real>::n_dim; j++)
      fprintf(fpout4, "%le\t", Kriging<Real>::Xinput[i][j]);
    for (j = 0; j < Kriging<Real>::ny; j++)
      fprintf(fpout4, "%le\t", Youtput[i][j]);
    fprintf(fpout4, "\n");
  }
  fprintf(fpout4, "\n\n");

  fprintf(fpout4, "VARIABLES=x1");
  for (i = 1; i < Kriging<Real>::n_dim; i++)
    fprintf(fpout4, ",x%d", i + 1);
  for (i = 0; i < Kriging<Real>::ny; i++)
    fprintf(fpout4, ",y%d", i + 1);
  fprintf(fpout4, "\nZONE T=\"Y_MSE\",I=%d\n", Kriging<Real>::ninterp);

  for (i = 0; i < Kriging<Real>::ninterp; i++) {
    for (j = 0; j < Kriging<Real>::n_dim; j++)
      fprintf(fpout4, "%le\t", Kriging<Real>::Xinput[i][j]);
    for (j = 0; j < Kriging<Real>::ny; j++)
      fprintf(fpout4, "%le\t", Yrmse[i][j]);
    fprintf(fpout4, "\n");
  }

  fclose(fpout4);
  for (i = 0; i < Kriging<Real>::ninterp; i++) {
    free(Youtput[i]);
    free(Yrmse[i]);
  }
  free(Youtput);
  free(Yrmse);
  free(xinterp);
}
