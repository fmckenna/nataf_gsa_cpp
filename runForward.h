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
#include "jsonInput.h"
extern std::ofstream theErrorFile;

using std::vector;

class runForward
{
public:
	runForward();
	runForward(vector<vector<double>> xval,
			vector<vector<double>> gval);
	~runForward();
	double writeOutputs(jsonInput inp);

	//vector<double> Si;

	vector<vector<double>> xval;
	vector<vector<double>> gval;

private:
	double calMean(vector<double> x);
	double calStd(vector<double> x, double m);
	double calSkewness(vector<double> x, double m, double s);
	double calKurtosis(vector<double> x, double m, double s);
	std::vector<double> mean;
	std::vector<double> stdDev;
	std::vector<double> skewness;
	std::vector<double> kurtosis;


	const double PI = 3.1415926535897932384626433;
	int nrv;
	int nmc;
};