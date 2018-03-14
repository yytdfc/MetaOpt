#include "GKrig.h"
#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

double GKrig::GKpredictorEI(double *x){
	double s, dy;
	GKsetPredict(0);
	dy = EI - GKpredictor(x);
	s = MSE(x);
	x[n_dim] = -dy*NormSDist(dy / s) - s*Normal(dy / s);
	if (EIcons){
		for (int i = 1; i < ny; i++){
			GKsetPredict(i);
			x[n_dim] = x[n_dim] * NormSDist(GKpredictor(x) / MSE(x));
			x[n_dim + i] = 1;
		}
	}
	else{
		for (int i = 1; i < ny; i++){
			GKsetPredict(i);
			x[n_dim + i] = GKpredictor(x);
		}
	}
	return x[n_dim];
}
double GKrig::GKpredictorMP(double *x){
	for (int i = 0; i < ny; i++){
		GKsetPredict(i);
		x[n_dim + i] = GKpredictor(x);
	}
	return x[n_dim];
}

void GKrig::setFx(double f(double*)){
	fx = f;
	isGK = 1;
}
double GKrig::fxIndex(double *x, int index){
	double *xx = new double[n_dim + ny];
	for (int i = 0; i < n_dim; ++i){
		xx[i] = x[i];
	}
	fx(xx);
	return(xx[n_dim + index]);
}

void GKrig::GKinitialize(int ncorr, int nconst_theta, int nporder,
	int nnorm, int ndcmp, int nParaOpt, int nregular, int ndim,
	int points, int nout_points, int nny, double **xx, double *up, double *low){
	initialize(ncorr, nconst_theta, nporder,
		nnorm, ndcmp, nParaOpt, nregular, ndim,
		points, nout_points, nny, xx, up, low);

}

void GKrig::GKsetPredict(int k){
	setPredict(k);
	if (isGK)
		yIndex = k;
}

/*******************************************************************************
* Func. :  two-lelvels hierachy kriging model initialization and fitting
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
void GKrig::GKtraining(){
	if (isGK){
		int i, j;
		/*----------------------------------------------------------------------------
		| Initialize regression  matrix F
		----------------------------------------------------------------------------*/
		init_F = 1;
		init_phi = 1;
		for (int k = 0; k < ny; ++k){
			F = nF[k];
			/*--------------------------------------------------------------------------
			| constant regression
			--------------------------------------------------------------------------*/
			for (i = 0; i < points; i++){
				if (flag[i] == 0){
					F[i][0] = fxIndex(xx[i], k);
				}
				else{
					cerr << "can't deal with gradient." << endl;
					F[i][0] = 0;
				}
			}/*for(i=0;i< krighf->points; i++)*/

			/*----------------------------------------------------------------------------
			| linear regression
			----------------------------------------------------------------------------*/
			if (porder >= 1){
				for (j = 1; j <= n_dim; j++){
					for (i = 0; i < points; i++){
						if (flag[i] == 0){
							F[i][j] = xx[i][j - 1]
								+ fxIndex(xx[i], k);
						}
						else if (j == flag[i])
							F[i][j] = 1.0;
						else
							F[i][j] = 0.0;
					}
				}
			}

			for (i = points; i < allpoints; i++)
				for (j = 0; j < np; j++)
					F[i][j] = 0.0;
		}
	}
	training();

}

/*******************************************************************************
* Func.  : two-level hierarchy kriging predictor, use it alfter model fitting
*          predict the respone at a untried site
*          paras is  kriging structure
*          vec_interp is untried site x to be predicted
* Author : Zhong-Hua.Han
* Date   : 29.06.2009
*******************************************************************************/
double GKrig::GKpredictor(double *vec_interp)
{
	if (isGK){
		int i;
		/*----------------------------------------------------------------------------
		| Initialize regression vector phi
		----------------------------------------------------------------------------*/
		init_phi = 1;
		phi[0] = fxIndex(vec_interp, yIndex);

		if (porder == 1)
		{
			for (i = 1; i <= n_dim; i++)
			{
				phi[i] = fxIndex(vec_interp, yIndex)
					*vec_interp[i - 1];
			}
		}
	}
	init_phi = 0;
	return(predictor(vec_interp));
} /** krig2h_predictor() **/

/*******************************************************************************
* Func.  : Mean Squared Error (MSE) estimation for two-levels hierarchy krigings
* Author : Zhong-Hua.Han
* Date   : 17.07.2009
*******************************************************************************/
double GKrig::GKMSE(double *vec_interp)
{
	if (isGK){
		int i;
		double rmse;
		init_phi = 1;
		phi[0] = fxIndex(vec_interp, yIndex);

		if (porder == 1)
		{
			for (i = 1; i <= n_dim; i++)
			{
				phi[i] = vec_interp[i - 1];
			}
			for (i = n_dim + 1; i < np; i++)
			{
				phi[i] = 0.0;
			}
		}
		init_phi = 0;
	}
	return(MSE(vec_interp));

} /** krig2h_MSE() **/

void GKrig::GKoutputRSM(){
	int i, j, k;
	//----------------write "output_rsm.dat"-------------------
	FILE *fpout;
	FILE *fpout2;
	if ((fpout = fopen("output_GK.dat", "w")) == NULL)
	{
		printf(" Failed to open file output_rsm.dat\n");
		exit(0);
	}
	if ((fpout2 = fopen("output_rsm_mse.dat", "w")) == NULL)
	{
		printf(" Failed to open file output_rsm_mse.dat\n");
		exit(0);
	}
	if (n_dim == 1)
	{
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
	}
	else if (n_dim == 2)
	{
		fprintf(fpout, "VARIABLES=x1,x2");
		for (i = 0; i < ny; i++)
			fprintf(fpout, ",y%d", i + 1);
		fprintf(fpout, "\n");
		fprintf(fpout, "ZONE T = \"Interpolated\",I=%d,J=%d\n", out_points, out_points);
		fprintf(fpout2, "VARIABLES=x1,x2");
		for (i = 0; i < ny; i++)
			fprintf(fpout2, ",y%d", i + 1);
		fprintf(fpout2, "\n");
		fprintf(fpout2, "ZONE T = \"Interpolated_MSE\",I=%d,J=%d\n", out_points, out_points);
	}
	else
	{
		; /** to be added... **/
	}
	double xout, xout1, yout, ymse;
	double *dx;
	dx = (double *)malloc((n_dim)* sizeof(double));
	for (i = 0; i < n_dim; i++)
	{
		if (out_points == 1)
			dx[i] = 0.0;
		else
			dx[i] = (xbound[i][1] - xbound[i][0]) / (out_points - 1.0);
	}

	double *xstar;
	xstar = (double *)malloc((n_dim)* sizeof(double));
	if (n_dim == 1)
	{
		for (i = 0; i < out_points; i++)
		{
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
			for (k = 0; k < ny; k++)
			{
				GKsetPredict(k);
				yout = GKpredictor(xstar);
				ymse = GKMSE(xstar);
				fprintf(fpout, "\t%le", yout);
				fprintf(fpout2, "\t%le", ymse);
			}
			fprintf(fpout, "\n");
			fprintf(fpout2, "\n");
		}
	}
	else if (n_dim == 2)
	{
		for (j = 0; j < out_points; j++)
			for (i = 0; i < out_points; i++)
			{
			xout = xbound[0][0] + dx[0] * i;
			xout1 = xbound[1][0] + dx[1] * j;
			/*-------------------------------------------------------------
			| normalization of sampling data
			-------------------------------------------------------------*/
			if (norm == 1)
			{
				xstar[0] = (xout - xbound[0][0]) / (xbound[0][1] - xbound[0][0]);
				xstar[1] = (xout1 - xbound[1][0]) / (xbound[1][1] - xbound[1][0]);
			}
			else
			{
				xstar[0] = xout;
				xstar[1] = xout1;
			}

			fprintf(fpout, "%le\t%le", xout, xout1);
			fprintf(fpout2, "%le\t%le", xout, xout1);
			for (k = 0; k < ny; k++)
			{
				GKsetPredict(k);
				yout = GKpredictor(xstar);
				ymse = GKMSE(xstar);
				fprintf(fpout, "\t%le", yout);
				fprintf(fpout2, "\t%le", ymse);
			}
			fprintf(fpout, "\n");
			fprintf(fpout2, "\n");
			}
	}
	else
	{
		;/** to be added ... **/
	}
	fclose(fpout);
	fclose(fpout2);

	free(xstar);
	free(dx);
}
void GKrig::GKprediction(){
	int i, j, k;

	/*----------------------------------------------------------------------------
	|  read parameters and sampled data from input file
	----------------------------------------------------------------------------*/
	readXinput();
	FILE *fpout4;
	if ((fpout4 = fopen("Youtput.dat", "w")) == NULL)
	{
		printf("Failed to open file Youtput.dat \n");
		exit(0);
	}
	double **Youtput, **Yrmse, *xinterp;
	Youtput = (double **)malloc((ninterp)* sizeof(double *));
	for (i = 0; i < ninterp; i++)
		Youtput[i] = (double *)malloc((ny)* sizeof(double));
	Yrmse = (double **)malloc((ninterp)* sizeof(double *));
	for (i = 0; i < ninterp; i++)
		Yrmse[i] = (double *)malloc((ny)* sizeof(double));
	xinterp = (double *)malloc((n_dim)* sizeof(double));
	for (k = 0; k < ny; k++)
	{
		GKsetPredict(k);
		for (i = 0; i < ninterp; i++)
		{
			for (j = 0; j < n_dim; j++)
			{
				if (norm == 1)
					xinterp[j] = (Xinput[i][j] - xbound[j][0])
					/ (xbound[j][1] - xbound[j][0]);
				else
					xinterp[j] = Xinput[i][j];
			}
			Youtput[i][k] = GKpredictor(xinterp);
			Yrmse[i][k] = GKMSE(xinterp);
		}
	}

	fprintf(fpout4, "VARIABLES=x1");
	for (i = 1; i < n_dim; i++)
		fprintf(fpout4, ",x%d", i + 1);
	for (i = 0; i < ny; i++)
		fprintf(fpout4, ",y%d", i + 1);
	fprintf(fpout4, "\nZONE T=\"Y_output\",I=%d\n", ninterp);

	for (i = 0; i < ninterp; i++)
	{
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

	for (i = 0; i < ninterp; i++)
	{
		for (j = 0; j < n_dim; j++)
			fprintf(fpout4, "%le\t", Xinput[i][j]);
		for (j = 0; j < ny; j++)
			fprintf(fpout4, "%le\t", Yrmse[i][j]);
		fprintf(fpout4, "\n");
	}

	fclose(fpout4);
	for (i = 0; i < ninterp; i++){
		free(Youtput[i]);
		free(Yrmse[i]);
	}
	free(Youtput);
	free(Yrmse);
	free(xinterp);

}