#include "optimizer/surro.h"
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
double prob_1_Rosenbrock(double *x){
	x[2] = (100 * (x[1] - x[0] * x[0])*(x[1] - x[0] * x[0]) + (1 - x[0])*(1 - x[0]));
	return x[2];
}
double prob_4(double *x){
	x[3] = x[0] + sin(5 * x[1]) - x[1] * x[1] + cos(7 * x[0]) - x[0] * x[1] + x[1] * x[2] - x[0] * x[2];
	x[4] = 2 * x[0] - 7 * x[1];
	x[5] = 3 * x[1] * x[1] - x[2];
	return  x[3];
}
double prob_2_G_9(double *x){
	x[8] = -(2 * pow(x[0], 2) + 3 * pow(x[1], 4) + x[2] + 4 * pow(x[3], 2) + 5 * x[4] - 127);
	x[9] = -(7 * x[0] + 3 * x[1] + 10 * pow(x[2], 2) + x[3] - x[4] - 282);
	x[10] = -(23 * x[0] + pow(x[1], 2) + 6 * pow(x[5], 2) - 8 * x[6] - 196);
	x[11] = -(4 * pow(x[0], 2) + pow(x[1], 2) - 3 * x[0] * x[1] + 2 * pow(x[2], 2) + 5 * x[5] - 11 * x[6]);
	x[7] = pow((x[0] - 10), 2) + 5 * pow((x[1] - 12), 2) + pow(x[2], 4) + 3 * pow((x[3] - 11), 2)
		+ 10 * pow(x[4], 6) + 7 * pow(x[5], 2) + pow(x[6], 4) - 4 * x[5] * x[6]
		- 10 * x[5] - 8 * x[6];
	return x[7];
}

double prob_3_dim2(double *x){
	const double pi = 3.1415926;
	x[3] = x[0] * x[1] - 0.2;
	x[2] = pow(15 * x[1] - 5.1 / 4 / pi / pi*pow(15 * x[0] - 5, 2) + 5 / pi*(15 * x[0] - 5) - 6, 2)
		+ 10 * ((1 - 1 / pi / 8)*cos(15 * x[0] - 5) + 1) + 5 * x[0];
	return  x[2];
}

double prob_100_dim100(double *x){
	x[100] = x[0] * x[0];
	for (int i = 1; i < 100; ++i){
		x[100] += x[i] * x[i] / (i + 1);
	}
	return  x[100];
}
double prob_101_dim100(double *x){
	x[100] = x[0] * x[0] / 2;
	for (int i = 1; i < 100; ++i){
		x[100] += x[i] * x[i] / 2;
	}
	return  x[100];
}

double xx(double *x){
	x[2] = x[1] + x[0];
	return x[2];
}

int main(int argc, char *argv[]){
	int prob;
	cout.setf(ios::scientific);
	Surro obj1;
	if (argc < 2){
		prob = obj1.readInput("KGopt2.in");
	}
	else{
		prob = obj1.readInput(argv[1]);
	}
	switch (prob){
	case 1:
		obj1.setFx(prob_1_Rosenbrock);
		break;
	case 2:
		obj1.setFx(prob_2_G_9);
		break;
	case 3:
		obj1.setFx(prob_3_dim2);
		break;
	case 4:
		obj1.setFx(prob_4);
		break;
	case 100:
		obj1.setFx(prob_100_dim100);
		break;
	case 101:
		obj1.setFx(prob_101_dim100);
		break;
	}
	obj1.opt();
}