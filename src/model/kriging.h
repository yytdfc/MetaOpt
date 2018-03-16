#ifndef METAOPT_KRIGING_H
#define METAOPT_KRIGING_H

#include <string>
#include <cmath>

using namespace std;

template <typename Real>
class Kriging
{
 public:
  /*******************************************************************************
   * Prototypes for kriging/GEK model
   *******************************************************************************/
  void initialize(int    ncorr,
                  int    nconst_theta,
                  int    nporder,
                  int    nnorm,
                  int    ndcmp,
                  int    nParaOpt,
                  int    nregular,
                  int    ndim,
                  int    npoints,
                  int    nout_points,
                  int    nny,
                  Real** xxx,
                  Real*  up,
                  Real*  low);
  void readInput(string);
  void readRestart(string);
  void writeRestart(string);

  void training();
  void setPredict(int);
  void prediction();
  void predictionDB();
  void outputRSM();
  void destructor();
  Real predictor(Real*);
  Real predictor_unNorm(Real*);
  Real grad_predictor(Real*, int);

  Real predictorMP(Real*);
  Real predictorEI(Real*);
  Real predictorPI(Real*);
  Real predictorME(Real*);
  Real predictorLCB(Real*);

  Real EI;
  int  EIcons;

 protected:
  void fitting();
  void allocate();
  Real MSE(Real*);
  void readXinput();
  /*----------------------------------------------------------------------------
    | about sampling data
    ----------------------------------------------------------------------------*/
  int n_dim; /* number of dimensions                            */
  int ny;    //		number of y_dimensions

  int    out_points;
  int    points; /* number of sample points                         */
  Real** xx;     /* sampling sites [point][dimension]               */
  Real** xbound;
  Real*  yy;   /* observed reponse [point]                        */
  Real** nyy;  //
  int*   flag; /* internal variable for GEK and cokriging
                  indicate which kind of varialbe one sample include
                  ="0" for origional varialbe
                  >"0" for gradient information
                  indicates with repect which dimension
                  the patial derivative is taken             */
  int norm;    /* indicate if the data is normalized              */

  /*----------------------------------------------------------------------------
    | about correlation matrix
    ----------------------------------------------------------------------------*/
  Real*   yvi;    /* internal variable for y^T*v^(-1)                */
  Real**  nyvi;   /*n dim for yvi	*/
  Real*   rvi;    /* internal variable for r(x)^T*v^(-1)             */
  Real**  v;      /* internal variable for correlation matrix        */
  Real*** nv;     //
  Real*   diag;   /* diagonal of cholesky deompostion of v           */
  Real**  ndiag;  //
  int     ninterp;
  Real**  Xinput;
  int     isreadXinput;
  int     porder; /* order of regression part                        */
  /* 0 for constant regression (ordinary kriging )   */
  /* 1 for linear regression (universal kriging )    */
  int np; /* number of extra lines in correlation matrix     */
  /* for regression part                             */
  int allpoints; /* Total number of lines for correlation matrix    */
  int psdf;      /* indicate whether v is positive definite         */
  int dcmp;      /* indicator for the type of decompostion of v     */
  /* 0 for LU decompostion                           */
  /* 1 for Cholesky decompostion                     */
  /* 2 for new Cholesky decompostion                 */
  int*  indx;   /* used for lu decomposition                       */
  int** nindx;  //

  int corr;    /* choice of correlation functions
                  corr = 0 : GEXP
                  corr = 1 : GAUSS
                  corr = 2 : CSPL
                  corr = 10: IMQ
                  corr = 11: MQ
                  corr = 12: TPS
                  corr = 13: POW                                  */
  int regular; /* indicator for type of regularization            */
  /* 0 no regularization                             */
  /* 1 constant regularization                       */
  /* 2 optimized regularization                      */

  /*----------------------------------------------------------------------------
    | For global trend function and unbiasedness condition
    ----------------------------------------------------------------------------*/
  Real*** nF;
  Real**  F;      /* Regression model for design matrix              */
  Real*   phi;    /* Regression model for correlation vector         */
  Real*   phibar; /* Regression model for correlation vector         */
  /* only for grad predictor                         */
  int init_F; /* Indicate if F has been initialized              */
  /* 0 : not initialized                             */
  /* 1 : initialized                                 */
  int init_phi; /* Indicate if phi has been initialized            */
  /* 0 : not initialized                             */
  /* 1 : initialized                                 */
  int init_phibar; /* Indicate if phibar has been initialized         */
  /* 0 : not initialized                             */
  /* 1 : initialized                                 */
  int ParaOpt; /* method of hyper parameter tuning                */
  /* 1 Hooke Jeeve local search                      */
  /* 2 gradient-based Quasi-Newton algorithm         */
  int const_theta; /* indicate whether we use                         */
  int isGK = 0;
  /*----------------------------------------------------------------------------
    | For vfm
    ----------------------------------------------------------------------------*/
  int cokrig; /* Indicator for kriging and cokriging             */
  /* "0" for bridge function-based vfm               */
  /* "1" for cokriging                              */

  const Real PI = 3.14159265358979323846264338;

  // normal distribution
  Real Normal(Real z) { return exp((-1) * z * z / 2) / sqrt(2 * PI); }
  Real NormSDist(const Real z) {
    // this guards against overflow
    if (z > 6) return 1;
    if (z < -6) return 0;
    const Real gamma = 0.231641900, a1 = 0.319381530, a2 = -0.356563782,
               a3 = 1.781477973, a4 = -1.821255978, a5 = 1.330274429;
    Real k = 1.0 / (1 + fabs(z) * gamma);
    Real n = k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5))));
    n = 1 - Normal(z) * n;
    if (z < 0) return 1.0 - n;
    return n;
  }

 private:
  int  predictor_vector();
  int  correlation_matrix();
  Real process_variance();

  Real correlation_determine_alpha();

  /*------------------------------------------------------------------------------
    |  Cal. correlation function
    ------------------------------------------------------------------------------*/
  Real correlation_of_r(Real*, Real*);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | First derivative of correlation function
    ------------------------------------------------------------------------------*/
  Real cross_correlation_of_r(Real*, Real*, int);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | second derivative of correlation function
    ------------------------------------------------------------------------------*/
  Real cross_correlation_of_r2(Real*, Real*, int);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | second derivative of correlation function
    ------------------------------------------------------------------------------*/
  Real cross_correlation_of_r22(Real*, Real*, int, int);

  /*------------------------------------------------------------------------------
    | for hyper parameter opt.
    ------------------------------------------------------------------------------*/
  void theta_init();
  Real MLE(Real*);
  Real bridgefunction_rho();
  Real cokriging_gamma();

  /*------------------------------------------------------------------------------
    | extral functions
    ------------------------------------------------------------------------------*/
  Real sign(Real);
  /*------------------------------------------------------------------------------
    | LU decomposition and solution of linear equations
    ------------------------------------------------------------------------------*/
  void ludcmp(Real**, int, int*, Real*);
  void lubksb(Real**, int, int*, Real[]);

  /*------------------------------------------------------------------------------
    | Cholesky decomposition and solution of linear equations A*x = b
    ------------------------------------------------------------------------------*/
  int  choldcmp(Real**, int, Real[]);
  void cholbksb(Real**, int, Real[], Real[], Real[]);
  void ncholdcmp(Real**, int, Real[]);
  void ncholbksb(Real**, int, Real[], Real[], Real[]);

  /*******************************************************************************
    Prototypes of functions
   *******************************************************************************/

  int HookeAndJeeves(Real*, int, Real*, Real*, int);

  int output_theta_design_space(Real* x,
                                int   dim,
                                Real (*fpointer)(Real*),
                                Real* lbd,
                                Real* ubd);

  Real mle_val_and_diff(Real* val, Real* dmle);
  int mle_beta_sigma(Real* val);
  int mle_optimization(Real* x);
  int mle_initialization(Real* x);
  int permutation_of_vector(Real* array, int n);

  /*----------------------------------------------------------------------------
    | For hyper parameters optimization, @R.Zimmermann, Dec. 2008
    ----------------------------------------------------------------------------*/
  Real*  theta;   /* correlation parameter [n_dim]                   */
  Real** ntheta;  //
  Real*  pk;      /* correlation parameter [n_dim]                   */
  Real** npk;

  Real mu; /* factor added to the diagonal of v               */
  /* for regularization                              */
  Real  sigma_sq;  /* process variance, depending on theta and regbeta*/
  Real* nsigma_sq; /*n dim for sigmasq*/
  Real  mle;       /* conentrated likelihood function                 */
  int   ncal_mle;  /* number of calculate mle                         */
  int   h_dim;     /* number of total hyper parameters                */
  Real  rho;       /* The scaling change between yhf and ylf          */
  /* rho, ylf and yhf are only for bridge function   */
  Real rho_flag; /* indicator for scaling between yhf and ylf       */
  /* 0 : no scaling                                  */
  /* 1 : constant, closed form                       */
  Real* ylf; /* Low fidelity value at samping sites             */
  Real* yhf; /* high fidelity value at samping sites            */
  /*----------------------------------------------------------------------------
    | For cokriging
    ----------------------------------------------------------------------------*/
  Real** GAMMA;  /* regression matrix for linear scaling            */
  Real   rhobar; /* Relaxision constant for cross-correlation       */
  Real   gamma;  /* Ratio of variance of hi- and low fidelity model
                      for correlogram of cokriging                    */
};

#endif  // METAOPT_KRIGING_H