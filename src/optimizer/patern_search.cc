#include "kriging.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <Windows.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
/*******************************************************************************
* Func. : HOOKE AND JEEVES METHOD
*         algorithm for multidimensional MINIMIZATION without using gradients
* Notes : SEE: Kowalik, J,Osborne M.R."Methods for unconstrained
*         optimization problems" in "modern analytic and computational
*         methods in science and mathematics"
*         R. Bellman, (ed.)
*         Elsevier New York 1968
*
*         implementation here according to
*         LOPHAVEN, S, NIELSEN, HB, SOENDERGAARD, J.
*         DACE - A MATLAB Kriging Toolbox,
*         Technical Report IMM-TR-2002-12
*         Technical University Of Denmark
*         http://www2.imm.dtu.dk/~hbn/dace/
*         refered as LNS 2002
*
*         return value is adress of final parameter array x
*
* Author: R. Zimmermann, Nov. 2008
*******************************************************************************/

int Kriging::HookeAndJeeves(double *x, int dim,
	double *lbd, double *ubd,
	int IS)

	/*------------------------------------------------------------------------------
	| Input
	|   dim      : dimension of hyper parameter *x
	|   *fpointer: the pointer for ln-likeilihood funtion
	|   *lbd     : the lower bounds of design variables
	|   *ubd     : the upper bounds of design variables
	|   *paras   : kriging structure
	|   IS       : number of interration
	| Output
	|   x        : design variables
	------------------------------------------------------------------------------*/
{
	int i, j, k, l, J;
	int d = dim;

	/*----------------------------------------------------------------------------
	| new_x : bool variable indicating if function value has to be calculated anew
	| f_x : function value in x
	| f_phi: function value in phi
	| f_x_temp2: function value in x_temp
	| *x_temp: emporary parameter for function value comparison
	| *x_temp2: temporary parameter for function value comparison
	| *phi : temporary parameter for function value comparison
	| *delta : array storing the step size in each dimension
	| *v : used for adjusting bad initial values
	| *bad_init: corresponds to N in LNS 2008,
	|            stores indices of bad starting values
	----------------------------------------------------------------------------*/
	int new_x = 1;
	double alpha;
	double swapvar;
	double f_x;
	double f_phi;
	double f_x_temp2;

	double *x_temp;
	double *x_temp2;
	double *phi;
	double *delta;
	double *v;
	int *bad_init;

	x_temp = (double *)malloc(d * sizeof(double));
	x_temp2 = (double *)malloc(d * sizeof(double));
	phi = (double *)malloc(d *sizeof(double));
	v = (double *)malloc(d * sizeof(double));

	delta = (double *)malloc(d * sizeof(double));
	bad_init = (int *)malloc(d * sizeof(int));

	/*----------------------------------------------------------------------------
	|  Corresponds to Alg. "start", LNS 2008, 6
	---------------------------------------------------------------------------*/

	/*------------------------------------------------------------------
	|  initialize arrays
	| "-1" indicates an empty cell
	-----------------------------------------------------------------*/
	for (k = 0; k < d; k++)
	{
		bad_init[k] = -1;
		v[k] = 1.0;
	}
	for (k = 0; k < d; k++)
	{
		if (lbd[k] == ubd[k])
		{
			delta[k] = 1;
			/*----------------------------------------------------
			| this is to keep x[k] fixed at the boundary
			---------------------------------------------------*/
			x[k] = ubd[k];
		}
		else
		{
			/*----------------------------------------------------
			| initial step size
			| initial value outside boundary|
			---------------------------------------------------*/
			delta[k] = pow(2.0, (double)(k + 1) / (d + 2));
			if ((x[k]<lbd[k]) || (x[k]>ubd[k]))
			{
				/*----------------------------------------
				| 0.125 = 1/8
				---------------------------------------*/
				x[k] = pow(lbd[k] * ubd[k] * ubd[k] * ubd[k] * ubd[k] * ubd[k] * ubd[k] * ubd[k],
					0.125);
				bad_init[k] = 1;
			}
		}/* if(lbd[k]==ubd[k]) */
	}/* for(k=0;k<d;k++) */

	/*------------------------------------------------------------------
	| treatment of bad initial values
	-----------------------------------------------------------------*/
	for (k = 0; k < d; k++)
	{
		if (bad_init[k] == 1)
		{
			printf("bad initial value k= %d \n", k);
			for (i = 0; i < d; i++)
			{
				x_temp[i] = x[i];
				x_temp2[i] = x_temp[i];
			}
			/*----------------------------------------
			|  0.0625 = 1/16
			---------------------------------------*/
			v[k] = 0.0625;
			alpha = log(lbd[0] / x[0]) / log(v[0]);
			/*----------------------------------------
			| find minimal alpha
			---------------------------------------*/
			for (j = 0; j < d; j++)
			{
				if (bad_init[j] == 1)
				{
					if (alpha < (log(lbd[j] / x[j]) / log(v[j])))
					{
						alpha = log(lbd[j] / x[j]) / log(v[j]);
					}
				}
			}

			for (j = 0; j < d; j++)
			{
				if (bad_init[j] == 1)
				{
					v[j] = pow(v[j], alpha*0.2);
				}
			}

			for (l = 1; l <= 4; l++)
			{
				for (j = 0; j < d; j++)
				{
					phi[j] = pow(v[j], (double)l)*x_temp[j];
				}

				f_phi = MLE(phi);
				f_x_temp2 = MLE(x_temp2);
				f_x = MLE(x);
				if (f_phi <= f_x_temp2)
				{
					for (j = 0; j < d; j++)
					{
						x_temp2[j] = phi[j];
					}

					if (f_phi < f_x)
					{
						for (j = 0; j < d; j++)
						{
							new_x = 1;
							x[j] = phi[j];
						}
						new_x = 1;
						J = j;
					}
					else
						/*------------------------------------
						| Stops the l-loop
						-----------------------------------*/
						l = 5;

				}/* if(f_phi<= f_x_temp2) */

			}/* for(l=1;l<=4;l++) */
			/*------------------------------------------
			| swaps step size values in array "delta"
			------------------------------------------*/
			swapvar = delta[0];
			delta[0] = delta[J];
			delta[J] = swapvar;

		} /* if(bad_init[k]==1) */

	}/* for(k=0; k<d ;k++ */
	/*----------------------------------------------------------------------------
	| End Alg. "start"
	---------------------------------------------------------------------------*/

	/*----------------------------------------------------------------------------
	| Start interation loop
	---------------------------------------------------------------------------*/
	for (i = 1; i <= IS; i++)
	{
		/*----------------------------------------------------------------
		| copy current x to x_temp
		---------------------------------------------------------------*/
		for (k = 0; k < d; k++)
		{
			x_temp[k] = x[k];
		}

		/*----------------------------------------------------------------
		| Alg. "explore", LNS 2008, 6
		| note if parameter is at boundary
		---------------------------------------------------------------*/

		int atbd = 0;

		for (j = 0; j < d; j++)
		{
			if (lbd[j] < ubd[j])
			{
				/*--------------------------------------------------
				| copy current x to phi
				-------------------------------------------------*/
				for (k = 0; k < d; k++)
					phi[k] = x[k];
				if (x[j] == lbd[j])
				{
					phi[j] = lbd[j] * sqrt(delta[j]);
					atbd = 1;
				}
				else if (x[j] == ubd[j])
				{
					phi[j] = ubd[j] / (sqrt(delta[j]));
					atbd = 1;
				}
				else
				{
					/*--------------------------------------------------
					| phi[j]:=min( x[j]*delta[j], ubd[j])
					-------------------------------------------------*/
					if (x[j] * delta[j] < ubd[j])
					{
						phi[j] = x[j] * delta[j];
					}
					else
					{
						phi[j] = ubd[j];
					}
					atbd = 0;
				}
			}/* if(lbd[j]<ubd[j]) */

			f_phi = MLE(phi);
			if (new_x)
				f_x = MLE(x);

			/*----------------------------------------------------
			| current phi is better solution than current x
			----------------------------------------------------*/
			if (f_phi < f_x)
			{
				for (k = 0; k < d; k++)
				{
					x[k] = phi[k];
				}
				new_x = 1;
			}
			else
			{
				if (atbd == 0)
				{
					/*------------------------------------------------
					| phi[j]:=max( x[j]/delta[j], lbd[j])
					------------------------------------------------*/
					if (x[j] / delta[j]>lbd[j])
						phi[j] = x[j] / delta[j];
					else
						phi[j] = lbd[j];

					f_phi = MLE(phi);

					if (new_x)
						f_x = MLE(x);
					/*------------------------------------------------
					| current phi is better solution than current x
					------------------------------------------------*/
					if (f_phi < f_x)
					{
						new_x = 1;
						for (k = 0; k < d; k++)
						{
							x[k] = phi[k];
						}
					}
					else
						new_x = 0;

				}/* if(atbd==0) */

			}/* if(f_phi< f_x) */

		}/* for(j=0;j<d;j++) */
		/*----------------------------------------------------------------
		| end Alg. "explore"
		---------------------------------------------------------------*/

		/*----------------------------------------------------------------
		| Alg. "move" LNS 2008, 6
		---------------------------------------------------------------*/

		int x_NOTEQUAL_x_temp = 0;
		int notstop = 0;
		/*----------------------------------------------------------------
		| check if array x is equal to array x_temp
		---------------------------------------------------------------*/
		for (k = 0; k < d; k++)
		{
			if (x[k] != x_temp[k])
			{
				x_NOTEQUAL_x_temp += 1;
				/*--------------------
				| exit loop
				-------------------*/
				break;
			}
		}
		/*----------------------------------------------------------------
		| i.e. x = x_temp
		---------------------------------------------------------------*/
		if (x_NOTEQUAL_x_temp == 0)
		{
			for (k = 0; k < d; k++)
				delta[k] = pow(delta[k], 0.2);
		}
		else
		{
			for (k = 0; k < d; k++)
			{
				v[k] = x[k] / x_temp[k];
			}
			notstop = 1;
			while (notstop)
			{
				for (k = 0; k < d; k++)
					phi[k] = x[k] * v[k];

				for (j = 0; j < d; j++)
				{
					if (phi[j] <= lbd[j])
					{
						phi[j] = lbd[j];
						/*--------------------
						| stopping criterion
						-------------------*/
						notstop = 0;
					}
					else if (phi[j] >= ubd[j])
					{
						phi[j] = ubd[j];
						notstop = 0;
					}
				}/* for(j=0;j<d;j++) * /

				 /*------------------------------------------------------------
				 | theta might have changed the loop step before
				 | so new computation is necessary
				 | phi might have changed, so new computation is necessary
				 | current phi is better solution than current x
				 -----------------------------------------------------------*/
				f_x = MLE(x);
				f_phi = MLE(phi);
				if (f_phi < f_x)
				{
					for (k = 0; k < d; k++)
					{
						x[k] = phi[k];
						v[k] = v[k] * v[k];
					}
				}
				else
					notstop = 0;
			}
			for (k = 0; k < d; k++)
				delta[k] = pow(delta[k], 0.25);
		}/* if(x_NOTEQUAL_x_temp ==0) */

		/*----------------------------------------------------------------
		| end Alg. "move"
		---------------------------------------------------------------*/

		/*----------------------------------------------------------------
		| Rotate delta
		---------------------------------------------------------------*/

		swapvar = delta[0];
		for (k = 0; k < (d - 1); k++)
			delta[k] = delta[k + 1];
		delta[d - 1] = swapvar;

	} /* for(i=1;i<=IS;i++ * /

	  /*----------------------------------------------------------------
	  | Compute function value after optimization process
	  ---------------------------------------------------------------*/
	f_x = MLE(x);

	/*----------------------------------------------------------------------------
	| free local memory
	----------------------------------------------------------------------------*/

	free(bad_init);
	free(x_temp);
	free(x_temp2);
	free(phi);
	free(v);
	free(delta);

	return(0); /*adress of array x is returned*/

}/** HookeAndJeeves() **/

/*******************************************************************************
* Func. : Output the design space of theta for 1 and 2 dimentional problem
* Author : Zhong-hua.Han
* Date : 30.06.2009
*******************************************************************************/

int output_theta_design_space(double *x, int dim,
	double(*fpointer)(double *),
	double *lbd, double *ubd)
{
	int i, j, k;

	int d = dim;

	int n = 10001;

	double *theta, mle;
	theta = (double *)malloc(d *sizeof(double));

	double *dtheta;
	dtheta = (double *)malloc(d*sizeof(double));
	for (i = 0; i < d; i++)
		dtheta[i] = (ubd[i] - lbd[i]) / (n - 1.0);

	FILE *fpout;
	if ((fpout = fopen("output_theta_design_space.dat", "w")) == NULL)
	{
		printf("Failed to open file output_theta_design_space.dat.dat\n");
		exit(0);
	}
	if (d == 1)
	{
		fprintf(fpout, "Variables=theta,MLE\n");

		for (k = 1; k < n; k++)
		{
			for (i = 0; i < d; i++)
			{
				theta[i] = lbd[i] + dtheta[i] * k;
				fprintf(fpout, " %lf ", theta[i]);
			}
			mle = (*fpointer)(theta);
			fprintf(fpout, " %lf\n", mle);
		}
	}
	else if (d == 2)
	{
		fprintf(fpout, "Variables=theat1,theta2,MLE\n");
		fprintf(fpout, "zone t = \"\", i=%d,j=%d\n", n - 1, n - 1);
		for (k = 1; k < n; k++)
			for (j = 1; j < n; j++)
			{
			theta[0] = lbd[0] + dtheta[0] * j;
			theta[1] = lbd[1] + dtheta[1] * k;
			for (i = 0; i < d; i++)
			{
				fprintf(fpout, " %lf ", theta[i]);
			}
			mle = (*fpointer)(theta);
			fprintf(fpout, " %lf\n", mle);
			}
	}

	fclose(fpout);

	free(theta);
	free(dtheta);

	return (0);

}/* output_theta_design_space() */

/*******************************************************************************
* Function: Determine starting value x=(theta,beta,sigma^2) for MLE
* Note    : all input data (x-sites and y-values@x-sites) has to be transformed
*           to [0,1], gradient-information accordingly, BEFORE function call!
* Date    : 2010-01-25  $
* Author  : Benjamin Rosenbaum $
*******************************************************************************/
int Kriging::mle_initialization(double *x)
{
	/*----------------------------------------------------------------------------
	| check for wrong choice of correlation function in GEK
	----------------------------------------------------------------------------*/

	int i;

	if (corr == 2)
	{
		for (i = 0; i < points; i++)
		{
			if (flag[i]>0)
			{
				printf("ParaOpt=2 can't be used for corr=2 in GEK.\n");
				printf("Use ParaOpt=1 or corr=3 instead!\n");
				exit(1);
			}
		}
	}


	int j, k, l;

	double *dummygrad = (double *)malloc((n_dim + 2) * sizeof(double));
	double mle;

	/*----------------------------------------------------------------------------
	| lbd/ubd = lower/upper boundary for initial theta-space
	----------------------------------------------------------------------------*/

	double *lbd = (double *)malloc(n_dim * sizeof(double));
	double *ubd = (double *)malloc(n_dim * sizeof(double));

	switch (corr)
	{
	case 1:
	{
		for (i = 0; i < n_dim; i++)
		{
			lbd[i] = 1.0;
			ubd[i] = 100.0;
		}
		break;
	}
	case 2:
	{
		for (i = 0; i < n_dim; i++)
		{
			lbd[i] = .1;
			ubd[i] = 10.0;
		}
		break;
	}
	case 3:
	{
		for (i = 0; i < n_dim; i++)
		{
			lbd[i] = .1;
			ubd[i] = 10.0;
		}
		break;
	}
	default:
	{
		printf("\n Currently there is no such correlation function !!!\n");
		return(1);
	}
	}/*end switch( corr)*/

	/*----------------------------------------------------------------------------
	| initial value for hyperparameter-optimization
	| evaluate mle over theta-grid and save best mle-value
	----------------------------------------------------------------------------*/

	int ctheta = 5;
	double *xtheta = (double *)malloc((n_dim + 2) * sizeof(double));

	for (i = 0; i < n_dim; i++)
		x[i] = lbd[i];

	switch (n_dim)
	{
	case 1:
	{
		double mle_min = 1.0e+12;
		for (i = 0; i < ctheta; i++)
		{
			xtheta[0] = lbd[0] + i*(ubd[0] - lbd[0]) / (ctheta - 1);
			mle_beta_sigma(xtheta);
			mle = mle_val_and_diff(xtheta, dummygrad);
			if (mle < mle_min)
			{
				mle_min = mle;
				for (k = 0; k < n_dim + 2; k++)
					x[k] = xtheta[k];
			}
		}
		break;
	}
	case 2:
	{
		double mle_min = 1.0e+12;
		for (i = 0; i < ctheta; i++)
		{
			for (j = 0; j < ctheta; j++)
			{
				xtheta[0] = lbd[0] + i*(ubd[0] - lbd[0]) / (ctheta - 1);
				xtheta[1] = lbd[1] + j*(ubd[1] - lbd[1]) / (ctheta - 1);
				mle_beta_sigma(xtheta);
				mle = mle_val_and_diff(xtheta, dummygrad);
				if (mle < mle_min)
				{
					mle_min = mle;
					for (k = 0; k < n_dim + 2; k++)
						x[k] = xtheta[k];
				}
			}
		}
		break;
	}
	case 3:
	{
		double mle_min = 1.0e+12;
		for (i = 0; i < ctheta; i++)
		{
			for (j = 0; j < ctheta; j++)
			{
				for (l = 0; l < ctheta; l++)
				{
					xtheta[0] = lbd[0] + i*(ubd[0] - lbd[0]) / (ctheta - 1);
					xtheta[1] = lbd[1] + j*(ubd[1] - lbd[1]) / (ctheta - 1);
					xtheta[2] = lbd[2] + l*(ubd[2] - lbd[2]) / (ctheta - 1);
					mle_beta_sigma( xtheta);
					mle = mle_val_and_diff(xtheta, dummygrad);
					if (mle < mle_min)
					{
						mle_min = mle;
						for (k = 0; k < n_dim + 2; k++)
							x[k] = xtheta[k];
					}
				}
			}
		}
		break;
	}
	default:
	{
		printf("( n_dim)>3 still to be implemented in ");
		printf("kriging_mle_initialisation\n");
		printf("Or use LHC sampling over theta-space!\n");
		return(1);
	}
	}/* end switch( n_dim)*/

	/*----------------------------------------------------------------------------
	| initial value for hyperparameter-optimization
	| evaluate mle over Latin Hypercube sampling and save best mle-value
	----------------------------------------------------------------------------*/
	/*
	int nLHC=50;

	double **xthetaLHC = (double**)malloc( n_dim*sizeof(double*));
	for(k=0;k< n_dim;k++)
	xthetaLHC[k]=(double*)malloc(nLHC*sizeof(double));

	for(k=0;k< n_dim;k++)
	{
	for(i=0;i<nLHC;i++)
	xthetaLHC[k][i]=lbd[k]+i*(ubd[k]-lbd[k])/(nLHC-1);
	permutation_of_vector(xthetaLHC[k],nLHC);
	}

	double *xtheta     = (double *)malloc(( n_dim+2) * sizeof(double));
	double mle_min=1.0e+12;

	for(i=0;i< n_dim;i++)
	x[i]=lbd[i];

	for(i=0;i<nLHC;i++)
	{
	for(k=0;k< n_dim;k++)
	xtheta[k]=xthetaLHC[k][i];
	kriging_mle_beta_sigma(paras, xtheta);
	mle=kriging_mle_val_and_diff(paras,xtheta,dummygrad);
	if( mle < mle_min)
	{
	mle_min=mle;
	for(k=0;k< n_dim;k++)
	x[k]=xtheta[k];
	}
	}

	kriging_mle_beta_sigma(x);
	*/

	/*----------------------------------------------------------------------------
	| free local memory
	----------------------------------------------------------------------------*/
	/*
	for(k=0;k< n_dim;k++)
	free(xthetaLHC[k]);
	free(xthetaLHC);
	*/
	free(dummygrad);
	free(ubd);
	free(lbd);
	free(xtheta);

	return(0);

}/** kriging_mle_initalization **/


/*******************************************************************************
* Function: Optimization algorithm for MLE problem
*           uses gradient-based Quasi-Newton-Method and Trust-Region-Ansatz
* Note    : All input data (x-sites and y-values@x-sites) has to be transformed
*           to [0,1], gradient-information accordingly, BEFORE function call!
*           Beta and sigma^2 are NOT being implicitely computed, but are
*           optimization variables and are determined by the optimization alg.
* Date    : 2010-06-09  $
* Author  : Benjamin Rosenbaum $
*******************************************************************************/
int Kriging::mle_optimization(double *x)
{

	/*----------------------------------------------------------------------------
	| check for wrong choice of correlation function in GEK
	----------------------------------------------------------------------------*/

	int i;

	if (corr == 2)
	{
		for (i = 0; i < points; i++)
		{
			if (flag[i]>0)
			{
				printf("ParaOpt=2 can't be used for corr=2 in GEK.\n");
				printf("Use ParaOpt=1 or corr=3 instead!\n");
				exit(1);
			}
		}
	}

	/*----------------------------------------------------------------------------
	| initialization for optimization algorithm
	----------------------------------------------------------------------------*/

	int j, k;
	int iteration;

	double mle;
	double mleneu;

	double *lbd = (double *)malloc(n_dim * sizeof(double));
	double *ubd = (double *)malloc(n_dim * sizeof(double));

	double* grad = (double *)malloc((n_dim + 2) * sizeof(double));
	double* gradneu = (double *)malloc((n_dim + 2) * sizeof(double));
	double* xneu = (double *)malloc((n_dim + 2) * sizeof(double));
	double* d = (double *)malloc((n_dim + 2) * sizeof(double));

	double t;  /* stepsize */
	double c = 1.0e-8;
	double gradnorm, dTgrad, sTy;

	double** H = (double **)malloc((n_dim + 2) * sizeof(double*));
	double** Hinv = (double **)malloc((n_dim + 2) * sizeof(double*));
	for (i = 0; i < (n_dim + 2); i++)
	{
		H[i] = (double *)calloc((n_dim + 2), sizeof(double));
		Hinv[i] = (double *)calloc((n_dim + 2), sizeof(double));
	}

	double* xtheta = (double *)malloc(n_dim * sizeof(double));

	int *indx_lu = (int *)malloc((n_dim + 2)* sizeof(int));
	double d_lu;
	double ared, pred;
	double dummy1, dummy2;
	double* pB = (double *)malloc((n_dim + 2) * sizeof(double));
	double* pU = (double *)malloc((n_dim + 2) * sizeof(double));
	double* p = (double *)malloc((n_dim + 2) * sizeof(double));
	double normpU, normpB, normp;
	double deltaTR;
	double tauTR;
	double rhoTR;
	double aa, ab, bb;
	int TRboundary;

	double deltaTRneu;
	double deltaTR_ubnd = 1.0;
	double deltaTR_lbnd = 1.0e-14;
	double etaTR = 0.125;

	deltaTR = 0.1;

	double* sBFGS = (double *)malloc((n_dim + 2) * sizeof(double));
	double* yBFGS = (double *)malloc((n_dim + 2) * sizeof(double));
	double* BsBFGS = (double *)malloc((n_dim + 2) * sizeof(double));

	int lbnd_violate;

	normp = 1.0;
	normpB = 1.0;
	normpU = 1.0;

	int fall = 1;

	FILE     *OptimizationFile;
	OptimizationFile = fopen("./OptimizationFile.dat", "w");
	fprintf(OptimizationFile, "   i  ");
	for (k = 0; k < n_dim; k++)
		fprintf(OptimizationFile, "theta%d        ", k + 1);
	fprintf(OptimizationFile, "beta          sigma^2       ");
	fprintf(OptimizationFile, "funval        gradnorm\n");

	/*----------------------------------------------------------------------------
	| initialization for optimization algorithm
	----------------------------------------------------------------------------*/

	iteration = 0;

	/* initialize Hessian approx. as identity */
	for (j = 0; j < (n_dim + 2); j++)
		H[j][j] = 1.0;

	mle = mle_val_and_diff( x, grad);
	gradnorm = 0.0;
	for (j = 0; j < (n_dim + 2); j++)
		gradnorm += grad[j] * grad[j];
	gradnorm = sqrt(gradnorm);

	/*----------------------------------------------------------------------------
	| start optimization loop
	----------------------------------------------------------------------------*/

	while (gradnorm > 1.0e-3)
	{
		iteration++;
		mle = mle_val_and_diff(x, grad);
		gradnorm = 0.0;
		for (j = 0; j < (n_dim + 2); j++)
			gradnorm += grad[j] * grad[j];
		gradnorm = sqrt(gradnorm);

		fprintf(OptimizationFile, "% 4d ", iteration);
		for (k = 0; k < n_dim + 2; k++)
			fprintf(OptimizationFile, "% 1.6le ", x[k]);
		fprintf(OptimizationFile, "% 1.6le % 1.6le \n", mle, gradnorm);

		deltaTRneu = deltaTR;

		/*-------------------------------------------------------------------------
		| solve quadratic subproblem:
		|   min {f + g^T p + 0.5 p^T H p}
		|   s.t. p^T p <= delta^2
		| according to dogleg-method
		-------------------------------------------------------------------------*/

		/* compute LU-decomp of Hessian H in matrix Hinv */
		for (j = 0; j < n_dim + 2; j++)
			for (k = 0; k < n_dim + 2; k++)
				Hinv[j][k] = H[j][k];
		ludcmp(Hinv, n_dim + 2, indx_lu, &d_lu);

		/* (a) perform full Quasi-Newton-step pB=-H^(-1)*g */
		for (j = 0; j < n_dim + 2; j++)
			pB[j] = -grad[j];
		lubksb(Hinv, n_dim + 2, indx_lu, pB);

		normpB = 0.0;
		for (j = 0; j < n_dim + 2; j++)
			normpB += pB[j] * pB[j];
		normpB = sqrt(normpB);

		/* (b) perform full Gradient-step pU=- (g^T g)/(g^T H g) * g  */
		dummy1 = 0.0;
		for (j = 0; j < n_dim + 2; j++)
			dummy1 += grad[j] * grad[j];
		dummy2 = 0.0;
		for (j = 0; j < n_dim + 2; j++)
			for (k = 0; k < n_dim + 2; k++)
				dummy2 += grad[j] * H[j][k] * grad[k];
		dummy1 /= dummy2;
		for (j = 0; j < n_dim + 2; j++)
			pU[j] = -(dummy1)*grad[j];

		normpU = 0.0;
		for (j = 0; j < n_dim + 2; j++)
			normpU += pU[j] * pU[j];
		normpU = sqrt(normpU);

		/* dogleg-method: descent p as a combination of pB and pU */
		if (normpB < deltaTR) /*case 1: pB in TR -> full Quasi-Newton step pB*/
		{
			for (j = 0; j < n_dim + 2; j++)
				p[j] = pB[j];
			normp = normpB;
			TRboundary = 0;
			fall = 1;
		}
		else
		{
			if (normpU >= deltaTR) /*case 2: pU out of TR -> Cauchy-point*/
			{
				for (j = 0; j < n_dim + 2; j++)
					p[j] = -(deltaTR / gradnorm)*grad[j];
				normp = deltaTR;
				TRboundary = 1;
				fall = 2;
			}
			else /*case3: dogleg-path -> go from pU to pB until TR is reached */
			{
				tauTR = 1.0;
				normp = 0.0;
				for (j = 0; j < n_dim + 2; j++)
				{
					p[j] = pU[j] + tauTR*(pB[j] - pU[j]);
					normp += p[j] * p[j];
				}
				normp = sqrt(normp);
				while (normp > deltaTR)
				{
					tauTR = tauTR*0.95;
					normp = 0.0;
					for (j = 0; j < n_dim + 2; j++)
					{
						p[j] = pU[j] + tauTR*(pB[j] - pU[j]);
						normp += p[j] * p[j];
					}
				}
				TRboundary = 1;
				fall = 3;
			}
		}/*endif*/

		/*-------------------------------------------------------------------------
		| compute next iterate xneu, projection to admissable set if needed
		-------------------------------------------------------------------------*/

		for (j = 0; j < n_dim + 2; j++)
			xneu[j] = x[j] + p[j];

		lbnd_violate = 0;
		for (k = 0; k < n_dim; k++)
			if (xneu[k] < 1.0e-12)            /* theta_k */
				lbnd_violate = 1;
		if (xneu[n_dim + 1] < 1.0e-12) /* sigma^2 */
			lbnd_violate = 1;

		while (lbnd_violate == 1)
		{
			for (j = 0; j < n_dim + 2; j++)
			{
				p[j] = 0.95*p[j];
				xneu[j] = x[j] + p[j];
			}
			lbnd_violate = 0;
			for (k = 0; k < n_dim; k++)
				if (xneu[k] < 1.0e-12)
					lbnd_violate = 1;
			if (xneu[n_dim + 1] < 1.0e-12)
				lbnd_violate = 1;
		}

		/*-------------------------------------------------------------------------
		| perfom TR-step,
		| if actual_reduction/predicted_reduction=rhoTR is sufficient large
		-------------------------------------------------------------------------*/

		pred = 0.0;
		for (j = 0; j < n_dim + 2; j++)
			for (k = 0; k < n_dim + 2; k++)
				pred -= 0.5*p[j] * H[j][k] * p[k];
		for (j = 0; j < n_dim + 2; j++)
			pred -= grad[j] * p[j];

		ared = mle - mle_val_and_diff(xneu, gradneu);

		rhoTR = (double)ared / (double)pred;

		if ((rhoTR<0.25 || ared<0.0) && normp>deltaTR_lbnd)
		{
			deltaTRneu = 0.25*normp;
		}
		else
		{
			if (rhoTR>0.75 && TRboundary > 0)
			{
				if (2.0*deltaTR<deltaTR_ubnd)
					deltaTRneu = 2.0*deltaTR;
				else
					deltaTRneu = deltaTR_ubnd;
			}
		}
		if (rhoTR>etaTR && ared >= 0.0)
		{
			/*-----------------------------------------------------------------------
			| update Hessian BFGS
			-----------------------------------------------------------------------*/
			if (normp > deltaTR_lbnd)
			{
				for (j = 0; j < n_dim + 2; j++)
				{
					sBFGS[j] = xneu[j] - x[j];
					yBFGS[j] = gradneu[j] - grad[j];
				}
				dummy1 = 0.0;
				for (j = 0; j < n_dim + 2; j++)
					for (k = 0; k < n_dim + 2; k++)
						dummy1 += sBFGS[j] * H[j][k] * sBFGS[k];
				dummy2 = 0.0;
				for (j = 0; j < n_dim + 2; j++)
					dummy2 += yBFGS[j] * sBFGS[j];
				for (j = 0; j < n_dim + 2; j++)
				{
					BsBFGS[j] = 0.0;
					for (k = 0; k < n_dim + 2; k++)
					{
						BsBFGS[j] += H[j][k] * sBFGS[k];
					}
				}
				for (j = 0; j < n_dim + 2; j++)
					for (k = 0; k < n_dim + 2; k++)
						H[j][k] = H[j][k] - (BsBFGS[j] * BsBFGS[k]) / dummy1
						+ (yBFGS[j] * yBFGS[k]) / dummy2;
			}

			/*-----------------------------------------------------------------------
			| save new iterate xneu in old iterate x
			-----------------------------------------------------------------------*/

			for (j = 0; j < n_dim + 2; j++)
				x[j] = xneu[j];
		}

		/* restart Hessian as Identity every 100 iterations */
		if (iteration % 100 == 0)
			for (j = 0; j < n_dim + 2; j++)
			{
			for (k = 0; k < n_dim + 2; k++)
				H[j][k] = 0.0;
			H[j][j] = 1.0;
			}

		deltaTR = deltaTRneu;

		if (deltaTR < deltaTR_lbnd)
		{
			/*deltaTR=deltaTR_lbnd;*/
			printf("\n optimization algorithm ended: no further decrease possible\n");
			break;
		}

		if (iteration > 1.0e3)
		{
			printf("\n optimization algorithm ended: max. iterations\n");
			break;
		}
	}

	if (gradnorm <= 1.0e-3)
		printf("\n optimization algorithm ended: norm(gradient)<1.0e-3\n");

	printf("\n optimal solution: theta =");
	for (k = 0; k < n_dim; k++)
		printf("% 1.6le ", x[k]);
	printf("\n                    beta =% 1.6le", x[n_dim]);
	printf("\n                 sigma^2 =% 1.6le\n\n", x[n_dim + 1]);

	fclose(OptimizationFile);

	for (k = 0; k < n_dim; k++)
		xtheta[k] = x[k];

	/*----------------------------------------------------------------------------
	|  Initialize additional components of paras (F,phi,phibar,yvi,rvi)
	----------------------------------------------------------------------------*/

	for (i = 0; i < allpoints; i++)
		yvi[i] = 0.0;

	for (i = 0; i < allpoints; i++)
		rvi[i] = 0.0;

	/*----------------------------------------------------------------------------
	| Initialize regression  matrix F
	----------------------------------------------------------------------------*/
	if (init_F == 0)
	{
		/*----------------------------------------------------------------
		| constant regression
		---------------------------------------------------------------*/
		for (i = 0; i < points; i++)
		{
			if (flag[i] == 0)
				F[i][0] = 1.0;
			else
				F[i][0] = 0.0;

		}/* for(i=0;i< points; i++) */

		/*----------------------------------------------------------------
		| linear regression
		----------------------------------------------------------------*/
		if (porder == 1)
		{
			for (j = 1; j <= n_dim; j++)
			{
				for (i = 0; i < points; i++)
				{
					if (flag[i] == 0)
					{
						F[i][j] = xx[i][j - 1];
					}
					else if (j == flag[i])
						F[i][j] = 1.0;
					else
						F[i][j] = 0.0;

				} /* for(i=0;i< points; i++) */

			} /* for(j=1;j<= n_dim;j++)*/

		} /* if( porder==1) */

		for (i = points; i < allpoints; i++)
			for (j = 0; j < np; j++)
				F[i][j] = 0.0;
		init_F = 1;

	}/* if( init_F) */

	/*----------------------------------------------------------------------------
	| Initialize regression vector phi
	----------------------------------------------------------------------------*/
	if (init_phi == 0 && porder == 0)
	{
		phi[0] = 1.0;

		for (i = 1; i < np; i++)
		{
			phi[i] = 0.0;
		}
		init_phi = 1;
	}

	/*----------------------------------------------------------------------------
	| Initialize regression vector phibar
	----------------------------------------------------------------------------*/
	if (init_phibar == 0 && porder == 0)
	{
		phibar[0] = 0.0;

		for (i = 1; i < np; i++)
		{
			phibar[i] = 0.0;
		}
		init_phibar = 1;
	}


	mle = MLE(xtheta); /* writes add. necessary infos to paras */

	/*----------------------------------------------------------------------------
	| free local memory
	----------------------------------------------------------------------------*/
	free(lbd);
	free(ubd);
	free(grad);
	free(gradneu);
	free(xneu);
	free(d);
	for (i = 0; i < n_dim + 2; i++)
	{
		free(H[i]);
		free(Hinv[i]);
	}
	free(xtheta);
	free(indx_lu);
	free(pB);
	free(pU);
	free(p);
	free(sBFGS);
	free(yBFGS);
	free(BsBFGS);
	return(0);

}/** kriging_mle_optimization**/


/*******************************************************************************
* Function: Evaluation and Derivation of MLE-function
*           w.r.t. val = [theta,beta,sigma^2]
*           mle-value is returned value
*           Gradient is stored in dmle
* Note:     Extension from former function
kriging_MLE(krig_paras *paras, double *val) from kriging_MAIN.c.
Uses arbitrary beta and sigma^2.
Computation of correlations according to
Koehler/Owen "Computer Experiments" in "Handbook of Statistics",p279
* Date:     2010-01-11
* Author:   Benjamin Rosenbaum
********************************************************************************/
double Kriging::mle_val_and_diff(double *val, double *dmle)
{

	int i, j, k, l, n, ns;
	int *indx;
	double d;

	n = points;

	double* dummyvector = (double*)malloc(n*sizeof(double));
	double* z = (double*)malloc(n*sizeof(double));
	double* vz = (double*)malloc(n*sizeof(double));

	double* residual = (double*)malloc(n*sizeof(double));

	/*----------------------------------------------------------------------------
	| logdetv : log of determinate of correlation matrix
	----------------------------------------------------------------------------*/

	double logdetv = 0.0;

	/*----------------------------------------------------------------------------
	| **v : A copy of(points x points) block of   v
	----------------------------------------------------------------------------*/
	double **v = (double **)malloc(n*sizeof(double *));
	double **v_save = (double **)malloc(n*sizeof(double *));
	double ***dv = (double***)malloc(n_dim*sizeof(double**));

	indx = (int *)malloc(n* sizeof(int));

	for (i = 0; i < n; i++)
		v[i] = (double *)malloc(n*sizeof(double));

	for (i = 0; i < n; i++)
		v_save[i] = (double *)malloc(n*sizeof(double));

	for (k = 0; k < n_dim; k++)
	{
		dv[k] = (double**)malloc(n*sizeof(double*));
		for (i = 0; i < n; i++)
			dv[k][i] = (double*)malloc(n*sizeof(double));
	}

	/*----------------------------------------------------------------------------
	| copy current parameter to kriging data structure
	| and reconsturct correlation matrix
	----------------------------------------------------------------------------*/
	for (i = 0; i < n_dim; i++)
		theta[i] = val[i];
	correlation_matrix();

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			v_save[i][j] = v[i][j];

	if (corr == 1)/* matrix regularization in case of gauss-corr */
		for (i = 0; i < n; i++)
			v[i][i] += ((10.0 + n)*1.0e-8);

	/*----------------------------------------------------------------------------
	| copy beta and sigma_sq
	----------------------------------------------------------------------------*/

	double beta = val[n_dim];
	double sigma_sq = val[n_dim + 1];

	/*----------------------------------------------------------------------------
	| copy  y to y
	----------------------------------------------------------------------------*/

	double* y = (double*)malloc(n*sizeof(double));

	for (i = 0; i < n; i++)
		y[i] = yy[i];

	/*----------------------------------------------------------------------------
	| copy corrmatrix to **v
	| dvk = derivation of v resp. to theta_k
	----------------------------------------------------------------------------*/

	double* dist = (double*)malloc(n_dim*sizeof(double));
	double* R = (double*)malloc(n_dim*sizeof(double));
	int* ai = (int*)malloc(n_dim*sizeof(int));
	int* aj = (int*)malloc(n_dim*sizeof(int));

	double* xi = (double*)malloc(n_dim*sizeof(double));

	switch (corr)
	{
	case 1: /* gauss-correlation */
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				v[i][j] = v[i][j];

				for (k = 0; k < n_dim; k++)
				{
					dist[k] = xx[j][k] - xx[i][k];
					R[k] = exp(-theta[k] * dist[k] * dist[k]);
				}

				for (k = 0; k < n_dim; k++)
				{
					ai[k] = 0;
					aj[k] = 0;
					if (flag[i] == (k + 1))
						ai[k] = 1;
					if (flag[j] == (k + 1))
						aj[k] = 1;
				}

				/* initialize dv[k][][]: -1 if yi is derivation, else 1 */
				if (flag[i] == 0)
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = 1.0;
				else
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = -1.0;

				/* compute derivations dv[k][][] */
				for (k = 0; k < n_dim; k++)/* product of (dim) 1d-correlations */
				{
					for (l = 0; l < n_dim; l++)/* derivation resp. to (dim) thetas */
					{
						switch (ai[k] + aj[k])
						{
						case 0: /* neither yi nor yj is derivation */
						{
							if (k == l)
								dv[l][i][j] *= R[k] * (-1.0)*dist[k] * dist[k];
							else
								dv[l][i][j] *= R[k];
							break;
						}
						case 1: /* either yi or yj is derivation */
						{
							if (k == l)
								dv[l][i][j] *= R[k] * (2.0* theta[k]
								* dist[k] * dist[k] * dist[k]
								- 2.0*dist[k]);
							else
								dv[l][i][j] *= R[k] * (-2.0* theta[k] * dist[k]);
							break;
						}
						case 2: /* both yi and yj are derivations */
						{
							if (k == l)
								dv[l][i][j] *= R[k] * (-4.0* theta[k] * theta[k]
								* dist[k] * dist[k] * dist[k] * dist[k]
								+ 10.0* theta[k] * dist[k] * dist[k]
								- 2.0);
							else
								dv[l][i][j] *= R[k] * (-2.0* theta[k])
								*(1.0
								- 2.0* theta[k] * dist[k] * dist[k]);
							break;
						}
						}/* end switch */
					}
				}

			}
		}
		break;

	case 2: /* spline-correlation */
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{

				v[i][j] = v[i][j];

				for (k = 0; k < n_dim; k++)
				{
					dist[k] = fabs(xx[j][k] - xx[i][k]);
					xi[k] = dist[k] * theta[k];
				}

				for (k = 0; k < n_dim; k++)
				{
					ai[k] = 0;
					aj[k] = 0;
					if (flag[i] == (k + 1))
						ai[k] = 1;
					if (flag[j] == (k + 1))
						aj[k] = 1;
				}
				/* initialize dv[k][][]: -1 if yi is derivation, else 1 */
				if (flag[i] == 0)
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = 1.0;
				else
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = -1.0;

				/* compute derivations dv[k][][] */
				for (k = 0; k < n_dim; k++)/* product of (dim) 1d-correlations */
				{
					for (l = 0; l < n_dim; l++)/* derivation resp. to (#dim) thetas */
					{
						switch (ai[k] + aj[k])
						{
						case 0: /* neither yi nor yj is derivation */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.2)
							{
								if (k == l)
									dv[l][i][j] *= (-30.0 + 90.0*xi[k])*xi[k] * dist[k];
								else
									dv[l][i][j] *= 1.0 + (-15.0 + 30.0*xi[k])*xi[k] * xi[k];
							}
							else if (xi[k] > 0.2&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= -3.75*(1.0 - xi[k])*(1.0 - xi[k])*dist[k];
								else
									dv[l][i][j] *= 1.25*(1.0 - xi[k])*(1.0 - xi[k])*(1.0 - xi[k]);
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						case 1: /* either yi or yj is derivation */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.2)
							{
								if (k == l)
									dv[l][i][j] *= (-60.0 + 270.0*xi[k])*xi[k]
									* sign(xx[j][k] - xx[i][k]);
								else
									dv[l][i][j] *= (-30.0 + 90.0*xi[k])*xi[k] * theta[k]
									* sign(xx[j][k] - xx[i][k]);
							}
							else if (xi[k] > 0.2&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= (1.0 - xi[k])*(-3.75*(1.0 - xi[k]) + 7.5*xi[k])
									*sign(xx[j][k] - xx[i][k]);
								else
									dv[l][i][j] *= -3.75*(1.0 - xi[k])*(1.0 - xi[k])* theta[k]
									* sign(xx[j][k] - xx[i][k]);
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						case 2: /* both yi and yj are derivations */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.2)
							{
								if (k == l)
									dv[l][i][j] *= (-60.0 + 540.0*xi[k])* theta[k];
								else
									dv[l][i][j] *= (-30.0 + 180.0*xi[k])* theta[k]
									* theta[k];
							}
							else if (xi[k] > 0.2&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= (15.0 - 22.5*xi[k])* theta[k];
								else
									dv[l][i][j] *= 7.5*(1.0 - xi[k])* theta[k]
									* theta[k];
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						}
					}
				}
			}
		}
		break;

	case 3: /* spline-correlation */
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{

				v[i][j] = v[i][j];

				for (k = 0; k < n_dim; k++)
				{
					dist[k] = fabs(xx[j][k] - xx[i][k]);
					xi[k] = dist[k] * theta[k];
				}

				for (k = 0; k < n_dim; k++)
				{
					ai[k] = 0;
					aj[k] = 0;
					if (flag[i] == (k + 1))
						ai[k] = 1;
					if (flag[j] == (k + 1))
						aj[k] = 1;
				}
				/* initialize dv[k][][]: -1 if yi is derivation, else 1 */
				if (flag[i] == 0)
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = 1.0;
				else
					for (k = 0; k < n_dim; k++)
						dv[k][i][j] = -1.0;

				/* compute derivations dv[k][][] */
				for (k = 0; k < n_dim; k++)/* product of (dim) 1d-correlations */
				{
					for (l = 0; l < n_dim; l++)/* derivation resp. to (#dim) thetas */
					{
						switch (ai[k] + aj[k])
						{
						case 0: /* neither yi nor yj is derivation */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.4)
							{
								if (k == l)
									dv[l][i][j] *= (-30.0 + (105.0 + (-97.5)*xi[k])*xi[k])*xi[k]
									* dist[k];
								else
									dv[l][i][j] *= 1.0 + (-15.0 + (35.0 + (-24.375)*xi[k])*xi[k])
									*xi[k] * xi[k];
							}
							else if (xi[k] > 0.4&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= ((-20.0 / 3.0) + (20.0 + (-20.0 + (20.0 / 3.0)
									*xi[k])*xi[k])*xi[k])*dist[k];
								else
									dv[l][i][j] *= (5.0 / 3.0)
									+ ((-20.0 / 3.0) + (10.0 + ((-20.0 / 3.0) + (5.0 / 3.0)
									*xi[k])*xi[k])*xi[k])*xi[k];
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						case 1: /* either yi or yj is derivation */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.4)
							{
								if (k == l)
									dv[l][i][j] *= (-60.0 + (315.0 + (-390.0)*xi[k])*xi[k])*xi[k]
									* sign(xx[j][k] - xx[i][k]);
								else
									dv[l][i][j] *= (-30.0 + (105.0 + (-97.5)*xi[k])*xi[k])*xi[k]
									* theta[k]
									* sign(xx[j][k] - xx[i][k]);
							}
							else if (xi[k] > 0.4&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= ((-20.0 / 3.0)
									+ (40.0 + (-60.0 + (80.0 / 3.0)*xi[k])*xi[k])*xi[k])
									*sign(xx[j][k] - xx[i][k]);
								else
									dv[l][i][j] *= ((-20.0 / 3.0)
									+ (20.0 + (-20.0 + (20.0 / 3.0)*xi[k])*xi[k])*xi[k])
									* theta[k]
									* sign(xx[j][k] - xx[i][k]);
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						case 2: /* both yi and yj are derivations */
						{
							if (xi[k] >= 0.0&&xi[k] <= 0.4)
							{
								if (k == l)
									dv[l][i][j] *= (-60.0 + (630.0 + (-1170.0)*xi[k])*xi[k])
									* theta[k];
								else
									dv[l][i][j] *= (-30.0 + (210.0 + (-292.5)*xi[k])*xi[k])
									* theta[k] * theta[k];
							}
							else if (xi[k] > 0.4&&xi[k] < 1.0)
							{
								if (k == l)
									dv[l][i][j] *= (40.0 + (-120.0 + (80.0)*xi[k])*xi[k])
									* theta[k];
								else
									dv[l][i][j] *= (20.0 + (-40.0 + (20.0)*xi[k])*xi[k])
									* theta[k] * theta[k];
							}
							else
							{
								if (k == l)
									dv[l][i][j] *= 0.0;
								else
									dv[l][i][j] *= 0.0;
							}
							break;
						}
						}
					}
				}
			}
		}
		break;
	default:
		printf("\n Currently there is no such correlation funciton !!!\n");
	}/* end switch */


	/*----------------------------------------------------------------------------
	| decomposition of v and computation of its determinate
	| LU decomp. could be relpaced by Cholesky decomp.
	| After call v is LU decomp of itself, diagonal of L is not stores
	| Compute log of determinant by adding the log values of the diagonal of R
	----------------------------------------------------------------------------*/
	ludcmp(v, n, indx, &d);

	/*----------------------------------------------------------------------------
	| detcormatrix = log(det(v))
	| log(det(v))' = 1/det(v) * det(v)'
	|              = 1/det(v) * tr(v'*v^(-1)) * det(v)
	|              = trace(v'*v^(-1))
	----------------------------------------------------------------------------*/

	logdetv = 0.0;
	for (i = 0; i < n; i++)
		logdetv += log(fabs((double)(v[i][i])));

	double* dlogdetv = (double*)malloc(n_dim*sizeof(double));

	for (k = 0; k < n_dim; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
				dummyvector[j] = dv[k][j][i];
			lubksb(v, n, indx, dummyvector);
			dlogdetv[k] += dummyvector[i];
		}
	}

	/*----------------------------------------------------------------------------
	| compute z=y-f*beta and vz=v^(-1)*z=v^(-1)*(y-f*beta)
	----------------------------------------------------------------------------*/

	for (i = 0; i < n; i++)
	{
		if (flag[i] == 0)
			z[i] = y[i] - 1.0*beta;
		else
			z[i] = y[i];
		vz[i] = z[i];
	}
	lubksb(v, n, indx, vz);

	/*----------------------------------------------------------------------------
	| compute MLE-value
	----------------------------------------------------------------------------*/

	double mle = 0.0;
	for (i = 0; i < n; i++)
		mle += z[i] * vz[i];
	mle /= (double)sigma_sq;
	mle += n*log(sigma_sq);
	mle += logdetv;

	/*----------------------------------------------------------------------------
	| derivation resp. to theta ==> dmle[k]
	----------------------------------------------------------------------------*/

	for (k = 0; k < n_dim; k++)
	{
		dmle[k] = 0.0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				dmle[k] -= vz[i] * dv[k][i][j] * vz[j];
		dmle[k] /= (double)sigma_sq;
		dmle[k] += dlogdetv[k];
	}

	/*----------------------------------------------------------------------------
	| derivation resp. to beta ==> dmle[n_dim]
	----------------------------------------------------------------------------*/

	dmle[n_dim] = 0.0;
	for (i = 0; i < n; i++)
	{
		if (flag[i] == 0)
			dmle[n_dim] -= vz[i];
	}
	dmle[n_dim] *= (2.0 / (double)sigma_sq);

	/*----------------------------------------------------------------------------
	| derivation resp. to sigma ==> dmle[n_dim+1]
	----------------------------------------------------------------------------*/

	dmle[n_dim + 1] = 0.0;
	for (i = 0; i < n; i++)
		dmle[n_dim + 1] -= z[i] * vz[i];
	dmle[n_dim + 1] /= (double)(sigma_sq*sigma_sq);
	dmle[n_dim + 1] += (n / (double)sigma_sq);

	/*----------------------------------------------------------------------------
	| free local memory
	----------------------------------------------------------------------------*/

	for (i = 0; i < n; i++)
	{
		free(v[i]);
		free(v_save[i]);
	}
	free(v);
	free(v_save);

	for (k = 0; k < n_dim; k++)
	{
		for (i = 0; i < n; i++)
		{
			free(dv[k][i]);
		}
		free(dv[k]);
	}
	free(dv);

	free(R);
	free(dist);

	free(ai);
	free(aj);

	free(indx);
	free(dummyvector);
	free(z);
	free(vz);
	free(y);

	free(dlogdetv);

	free(residual);

	free(xi);

	return(mle);

} /** kriging_mle_val_and_diff() **/

/*******************************************************************************
* Function: Compute beta(theta) and sigma^2(beta,theta)
* Note    : Is NOT used in the optimation-routine,
*           only for initial value and checking.
* Date    : 2010-01-11  $
* Author  : Benjamin Rosenbaum $
*******************************************************************************/
int Kriging::mle_beta_sigma(double *val)
{

	int i, j;
	int *indx;
	double d;

	int n = points;

	double* rz = (double *)malloc(n* sizeof(double));
	double dummy1, dummy2;

	/*----------------------------------------------------------------------------
	| **v : A copy of(points x points) block of   v
	----------------------------------------------------------------------------*/
	double **v;

	indx = (int *)malloc(n* sizeof(int));

	v = (double **)malloc(n*sizeof(double *));
	for (i = 0; i < n; i++)
		v[i] = (double *)malloc(n*sizeof(double));

	/*----------------------------------------------------------------------------
	| copy current parameter to kriging data structure
	| and reconsturct correlation matrix
	----------------------------------------------------------------------------*/
	for (i = 0; i < n_dim; i++)
		theta[i] = val[i];

	correlation_matrix();

	if (corr == 1)
		for (i = 0; i < n; i++)
			v[i][i] += ((10.0 + n)*1.0e-8);

	/*----------------------------------------------------------------------------
	| copy  y to y
	----------------------------------------------------------------------------*/
	double* y = (double*)malloc(n*sizeof(double));

	for (i = 0; i < n; i++)
		y[i] = yy[i];

	/*----------------------------------------------------------------------------
	| copy corrmatrix to **v
	----------------------------------------------------------------------------*/
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			v[i][j] = v[i][j];
	}
	/*----------------------------------------------------------------------------
	| decomposition of v and computation of its determinate
	| After call v is LU decomp of itself, diagonal of L is not stores
	----------------------------------------------------------------------------*/
	ludcmp(v, n, indx, &d);

	/*----------------------------------------------------------------------------
	| compute rf=R^(-1)*f, ry=R^(-1)*y
	----------------------------------------------------------------------------*/
	double* rf = (double*)malloc(n*sizeof(double));
	double* ry = (double*)malloc(n*sizeof(double));

	for (i = 0; i < n; i++)
	{
		if (flag[i] == 0)
			rf[i] = 1.0;
		else
			rf[i] = 0.0;
		ry[i] = y[i];
	}
	lubksb(v, n, indx, rf);
	lubksb(v, n, indx, ry);

	/*----------------------------------------------------------------------------
	| compute beta=val[n_dim]
	----------------------------------------------------------------------------*/
	dummy1 = 0.0;
	for (i = 0; i < n; i++)
		if (flag[i] == 0)
			dummy1 += 1.0*rf[i];
	dummy2 = 0.0;
	for (i = 0; i < n; i++)
		if (flag[i] == 0)
			dummy2 += 1.0*ry[i];
	val[n_dim] = dummy2 / dummy1;

	/*----------------------------------------------------------------------------
	| compute sigma^2=val[n_dim+1]
	----------------------------------------------------------------------------*/
	for (i = 0; i < n; i++)
		rz[i] = ry[i] - val[2] * rf[i];

	val[n_dim + 1] = 0.0;
	for (i = 0; i < n; i++)
	{
		if (flag[i] == 0)
			val[n_dim + 1] += (y[i] - val[n_dim])*rz[i];
		else
			val[n_dim + 1] += (y[i])*rz[i];
	}

	val[n_dim + 1] = val[n_dim + 1] / ((double)(n));

	/*----------------------------------------------------------------------------
	| free local memory
	----------------------------------------------------------------------------*/

	for (i = 0; i < n; i++)
		free(v[i]);
	free(v);
	free(rf);
	free(ry);
	free(rz);
	free(y);
	free(indx);

	return(0);

}/** kriging_mle_beta_sigma **/

/*******************************************************************************
* Function: random permutation of vector entries
* Date    : 2009-10-21  $
* Author  : Benjamin Rosenbaum $
*******************************************************************************/
int permutation_of_vector(double* array, int n)
{
	int i, j;
	double dummy;
	struct timeval zeit;
	int     gettimeofday(struct timeval *tp, void *tzp);

	for (i = 0; i < n; i++)
	{
		gettimeofday(&zeit, NULL);
		srand(zeit.tv_usec);
		j = rand() % n;
		dummy = array[i];
		array[i] = array[j];
		array[j] = dummy;
	}
	return 0;
}/** permutation_of_vector **/
int gettimeofday(struct timeval *tp, void *tzp)
{
	time_t clock;
	struct tm tm;
	SYSTEMTIME wtm;

	GetLocalTime(&wtm);
	tm.tm_year = wtm.wYear - 1900;
	tm.tm_mon = wtm.wMonth - 1;
	tm.tm_mday = wtm.wDay;
	tm.tm_hour = wtm.wHour;
	tm.tm_min = wtm.wMinute;
	tm.tm_sec = wtm.wSecond;
	tm.tm_isdst = -1;
	clock = mktime(&tm);
	tp->tv_sec = clock;
	tp->tv_usec = wtm.wMilliseconds * 1000;

	return (0);
}

