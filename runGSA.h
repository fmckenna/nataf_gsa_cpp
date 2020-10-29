#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <vector>
#include <numeric> // std::iota
#include <armadillo>
extern std::ofstream theErrorFile;
using namespace arma;
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
using std::vector;

class runGSA
{
public:
	runGSA();
	runGSA(vector<vector<double>> xval,
			vector<vector<double>> gval,
			vector<vector<int>> combs_tmp,
			int Kos);
	~runGSA();
	//vector<double> Si;

	vector<vector<double>> xval;
	vector<vector<int>> combs_tmp;
	char Opt;
	//int Kos;
	vector<vector<double>> Simat;
	vector<vector<double>> Stmat;

private:
	double mvnPdf(mat x, mat mu, mat cov);
	double calMean(vector<double> x);
	double calVar(vector<double> x);
	const double PI = 3.1415926535897932384626433;
	vector<double> doGSA(vector<double> gval, int Kos, char Opt);
	int nrv;
	int ncombs;
	int nmc;
};