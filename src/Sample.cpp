#include "Sample.h"
#include <iostream>
#include <fstream>
using namespace std;

Sample::Sample(int n, int num){
	ndim = n;
	nSample = num;
	sample = new double*[ndim];
	for (int i = 0; i < ndim; ++i){
		sample[i] = new double[nSample];
	}
	ifadd = 0;
	hasSample = 0;
	gen.seed(rd());
}

Sample::~Sample(){
	for (int i = 0; i < ndim; ++i){
		delete [] sample[i];
	}
	delete [] sample;
	if (ifadd){
		for (int i = 0; i < ndim; ++i){
			delete[] addSample[i];
		}
		delete[] addSample;
	}
}

void Sample::printSample(string str){
	if (!hasSample){
		cerr << "doesn't have sample!" << endl;
		return;
	}
	int i, j;
	ofstream fout(str);
	fout.setf(ios::scientific);
	for (i = 0; i < nSample; ++i){
		for (j = 0; j < ndim; ++j){
			fout << sample[j][i] << "\t";
		}
		fout << endl;
	}
	fout.close();
}

void Sample::printAdd(string str){
	if (!ifadd){
		cerr << "doesn't have sample!" << endl;
		return;
	}
	int i, j;
	ofstream fout(str);
	fout.setf(ios::scientific);
	for (i = 0; i < addNum; ++i){
		for (j = 0; j < ndim; ++j){
			fout << addSample[j][i] << "\t";
		}
		fout << endl;
	}
	fout.close();
}

void Sample::genLHS(){
	int i, j;
	int *list = new int [nSample];
	int rand;
	for (i = 0; i < ndim; ++i){
		for (j = 0; j < nSample; ++j){
			list[j] = j;
		}
		for (j = 0; j < nSample; ++j){
			rand = randomInt(j, nSample - 1);
			sample[i][j] = list[rand];
			sample[i][j] = (sample[i][j] + randomDouble()) / nSample;
			list[rand] = list[j];
		}
	}
	hasSample = 1;
	delete[]list;
}

void Sample::genMC(){
	int i, j;
	for (i = 0; i < ndim; ++i){
		for (j = 0; j < nSample; ++j){
			sample[i][j] = randomDouble();
		}
	}
	hasSample = 1;
}

void Sample::readSample(string str){
	int i, j;
	ifstream fin(str);
	for (i = 0; i<nSample; ++i){
		for (j = 0; j<ndim; ++j){
			fin >> sample[j][i];
		}
	}
	fin.close();
	hasSample = 1;
}

void Sample::addLHS(int num){
	if (!hasSample){
		cerr << "doesn't have sample!" << endl;
		return;
	}	
	int i, j, k;
	int position, rand;
	addNum = num;
	addSample = new double *[ndim];
	for (i = 0; i < ndim; ++i){
		addSample[i] = new double[addNum];
	}
	int *map = new int[nSample + addNum];
	int *list = new int [addNum];
	for (i = 0; i < ndim; ++i){
		// initial map
		for (j = 0; j < nSample + addNum; ++j){
			map[j] = 0;
		}
		// mark exist points
		for (j = 0; j < nSample; ++j){
			position = (int)(sample[i][j] * (nSample + addNum));
			for (k = 0; k < nSample + addNum; ++k){
				if (map[position + k] == 0 && (position + k) < (nSample + addNum)){
					map[position + k] = 1;
					break;
				}
				else if (map[position - k] == 0 && (position - k) >= 0){
					map[position - k] = 1;
					break;
				}

			}
		}
		// get the blank points sequence
		for (j = 0, k = 0; j < nSample + addNum; ++j){
			if (map[j] == 0){
				list[k] = j;
				++k;
			}
		}
		// get LHS for add samples
		for (j = 0; j < addNum; ++j){
			rand = randomInt(j, addNum - 1);
			addSample[i][j] = list[rand];
			addSample[i][j] = (addSample[i][j] + randomDouble()) / (nSample + addNum);
			list[rand] = list[j];
		}
	}
	ifadd = 1;
	delete[] map;
	delete[] list;
}
