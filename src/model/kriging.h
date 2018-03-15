#include <string>
#include <cmath>

using namespace std;

class Kriging
{
 public:
  /*******************************************************************************
   * Prototypes for kriging/GEK model
   *******************************************************************************/
  void initialize(int      ncorr,
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
                  double*  low);
  void readInput(string);
  void readRestart(string);
  void writeRestart(string);

  void   training();
  void   setPredict(int);
  void   prediction();
  void   predictionDB();
  void   outputRSM();
  void   destructor();
  double predictor(double*);
  double predictor_unNorm(double*);
  double grad_predictor(double*, int);

  double predictorMP(double*);
  double predictorEI(double*);
  double predictorPI(double*);
  double predictorME(double*);
  double predictorLCB(double*);

  double EI;
  int    EIcons;

 protected:
  void   fitting();
  void   allocate();
  double MSE(double*);
  void   readXinput();
  /*----------------------------------------------------------------------------
    | about sampling data
    ----------------------------------------------------------------------------*/
  int n_dim; /* number of dimensions                            */
  int ny;    //		number of y_dimensions

  int      out_points;
  int      points; /* number of sample points                         */
  double** xx;     /* sampling sites [point][dimension]               */
  double** xbound;
  double*  yy;   /* observed reponse [point]                        */
  double** nyy;  //
  int*     flag; /* internal variable for GEK and cokriging
                    indicate which kind of varialbe one sample include
                    ="0" for origional varialbe
                    >"0" for gradient information
                    indicates with repect which dimension
                    the patial derivative is taken             */
  int norm;      /* indicate if the data is normalized              */

  /*----------------------------------------------------------------------------
    | about correlation matrix
    ----------------------------------------------------------------------------*/
  double*   yvi;    /* internal variable for y^T*v^(-1)                */
  double**  nyvi;   /*n dim for yvi	*/
  double*   rvi;    /* internal variable for r(x)^T*v^(-1)             */
  double**  v;      /* internal variable for correlation matrix        */
  double*** nv;     //
  double*   diag;   /* diagonal of cholesky deompostion of v           */
  double**  ndiag;  //
  int       ninterp;
  double**  Xinput;
  int       isreadXinput;
  int       porder; /* order of regression part                        */
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
  double*** nF;
  double**  F;      /* Regression model for design matrix              */
  double*   phi;    /* Regression model for correlation vector         */
  double*   phibar; /* Regression model for correlation vector         */
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

  const double PI = 3.14159265358979323846264338;

  // normal distribution
  double Normal(double z) { return exp((-1) * z * z / 2) / sqrt(2 * PI); }
  double NormSDist(const double z) {
    // this guards against overflow
    if (z > 6) return 1;
    if (z < -6) return 0;
    const double gamma = 0.231641900, a1 = 0.319381530, a2 = -0.356563782,
                 a3 = 1.781477973, a4 = -1.821255978, a5 = 1.330274429;
    double k = 1.0 / (1 + fabs(z) * gamma);
    double n = k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5))));
    n = 1 - Normal(z) * n;
    if (z < 0) return 1.0 - n;
    return n;
  }

 private:
  int    predictor_vector();
  int    correlation_matrix();
  double process_variance();

  double correlation_determine_alpha();

  /*------------------------------------------------------------------------------
    |  Cal. correlation function
    ------------------------------------------------------------------------------*/
  double correlation_of_r(double*, double*);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | First derivative of correlation function
    ------------------------------------------------------------------------------*/
  double cross_correlation_of_r(double*, double*, int);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | second derivative of correlation function
    ------------------------------------------------------------------------------*/
  double cross_correlation_of_r2(double*, double*, int);

  /*------------------------------------------------------------------------------
    | for cross correlation of GEK
    | second derivative of correlation function
    ------------------------------------------------------------------------------*/
  double cross_correlation_of_r22(double*, double*, int, int);

  /*------------------------------------------------------------------------------
    | for hyper parameter opt.
    ------------------------------------------------------------------------------*/
  void   theta_init();
  double MLE(double*);
  double bridgefunction_rho();
  double cokriging_gamma();

  /*------------------------------------------------------------------------------
    | extral functions
    ------------------------------------------------------------------------------*/
  double sign(double);
  /*------------------------------------------------------------------------------
    | LU decomposition and solution of linear equations
    ------------------------------------------------------------------------------*/
  void ludcmp(double**, int, int*, double*);
  void lubksb(double**, int, int*, double[]);

  /*------------------------------------------------------------------------------
    | Cholesky decomposition and solution of linear equations A*x = b
    ------------------------------------------------------------------------------*/
  int  choldcmp(double**, int, double[]);
  void cholbksb(double**, int, double[], double[], double[]);
  void ncholdcmp(double**, int, double[]);
  void ncholbksb(double**, int, double[], double[], double[]);

  /*******************************************************************************
    Prototypes of functions
   *******************************************************************************/

  int HookeAndJeeves(double*, int, double*, double*, int);

  int output_theta_design_space(double* x,
                                int     dim,
                                double (*fpointer)(double*),
                                double* lbd,
                                double* ubd);

  double mle_val_and_diff(double* val, double* dmle);
  int mle_beta_sigma(double* val);
  int mle_optimization(double* x);
  int mle_initialization(double* x);
  int permutation_of_vector(double* array, int n);

  /*----------------------------------------------------------------------------
    | For hyper parameters optimization, @R.Zimmermann, Dec. 2008
    ----------------------------------------------------------------------------*/
  double*  theta;   /* correlation parameter [n_dim]                   */
  double** ntheta;  //
  double*  pk;      /* correlation parameter [n_dim]                   */
  double** npk;

  double mu; /* factor added to the diagonal of v               */
  /* for regularization                              */
  double  sigma_sq;  /* process variance, depending on theta and regbeta*/
  double* nsigma_sq; /*n dim for sigmasq*/
  double  mle;       /* conentrated likelihood function                 */
  int     ncal_mle;  /* number of calculate mle                         */
  int     h_dim;     /* number of total hyper parameters                */
  double  rho;       /* The scaling change between yhf and ylf          */
  /* rho, ylf and yhf are only for bridge function   */
  double rho_flag; /* indicator for scaling between yhf and ylf       */
  /* 0 : no scaling                                  */
  /* 1 : constant, closed form                       */
  double* ylf; /* Low fidelity value at samping sites             */
  double* yhf; /* high fidelity value at samping sites            */
  /*----------------------------------------------------------------------------
    | For cokriging
    ----------------------------------------------------------------------------*/
  double** GAMMA;  /* regression matrix for linear scaling            */
  double   rhobar; /* Relaxision constant for cross-correlation       */
  double   gamma;  /* Ratio of variance of hi- and low fidelity model
                      for correlogram of cokriging                    */
};
