

#include "chiSquaredDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlopt.hpp"
#include "boost/math/distributions/normal.hpp" // for normal_distribution

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlChiSq(unsigned n, const double* x, double* grad, void* my_func_data);

chiSquaredDist::chiSquaredDist(string opt, vector<double> val, vector<double> add): chiSqDist(1.0)
{
	name = "chisquared";
	int npa = 1;
	if (opt.compare("PAR") == 0)
	{
		if (val.size() != npa)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= 0) || (val[0]-floor(val[0])>1.e-10)) // not integer
		{
			theErrorFile << "Error running UQ engine: parameter of " << name << " distribution must a positive integer " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			k = val[0];
		}
	}
	else if (opt.compare("MOM") == 0)
	{

		if (val.size() != npa)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= 0) || (val[0] - floor(val[0]) > 1.e-10)) // not integer
		{
			theErrorFile << "Error running UQ engine: mean of " << name << " distribution must a positive integer " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			k = val[0];
		}

	}
	else if (opt.compare("DAT") == 0)
	{
		theErrorFile << "The Chisquare distribution is not supported in DATA input type" << std::endl;
		theErrorFile.close();
		exit(-1);
		/*
		const int np = 1;
		//double lb[np] = { -HUGE_VAL, 0 };
		//const double  *lb = { 0 };

		// initial lambda
		double mu = 0.0, var = 0.0, ns=0, std;
		for(double value : val)
		{
			mu += value;
			//var += value*value;
			ns += 1;
		}
		k = mu;

		vector<double> xs = val;

		// MLE optimization
		double lb[1] = { 1.e-10 };
		nlopt_opt optim;
		optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		nlopt_set_min_objective(optim, nnlChiSq, &xs);
		nlopt_set_xtol_rel(optim, 1e-6);
		double minf; 

		if (nlopt_optimize(optim, &k, &minf) < 0) {
			printf("nlopt failed!\n");
			theErrorFile << "Error running UQ engine: MLE optimization filed" << std::endl;
			theErrorFile.close();
			exit(1);
		}
		else {
			//lambda = x[0];
			printf("found minimum at f(%g) = %0.10g\n", k, minf);
		}
		*/
	}

	chi_squared chiSqDist1(k);
	chiSqDist = chiSqDist1;
	checkParams();

}

chiSquaredDist::~chiSquaredDist() {}
//==

void chiSquaredDist::checkParams()
{
	double std = getStd();
	double mean = getMean();
	vector<double> par = getParam();

	if (isnan(std) || isinf(std) || std <= 0)
	{
		theErrorFile << "Error running UQ engine: stdandard deviation of " << name << " must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	if ((mean <= 0) || (mean - floor(mean) > 1.e-10)) // not integer
	{
		theErrorFile << "Error running UQ engine: mean of " << name << " must a positive integer " << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	if ((par[0] <= 0) || (par[0] - floor(par[0]) > 1.e-10)) // not integer
	{
		theErrorFile << "Error running UQ engine: parameter of " << name << " must a positive integer " << std::endl;
		theErrorFile.close();
		exit(-1);
	}
}


double chiSquaredDist::getPdf(double x)
{
	//double result;
	//if (0.0 < x) {
	//	result = lambda * chiSq(-lambda * x);
	//}
	//else {
	//	result = 0.0;
	//}
	//return result;
	return pdf(chiSqDist,x);

}

double chiSquaredDist::getCdf(double x)
{
	/*
	double result;
	if (0.0 < x) {
		result = 1.0-chiSq(-(lambda) * x);
	}
	else {
		result = 0.0;
	}
	return result;
	*/
	return cdf(chiSqDist, x);

}


double chiSquaredDist::getMean(void)
{
	//return moms[0];
	double a =mean(chiSqDist);
	return mean(chiSqDist);
}

double chiSquaredDist::getStd(void)
{
	//return moms[1];
	return standard_deviation(chiSqDist);
}

double chiSquaredDist::getQuantile(double p)
{

	return quantile(chiSqDist, p);

}

vector<double> chiSquaredDist::getParam(void)
{
	return { k };
}


string chiSquaredDist::getName(void)
{
	return name;
}


double nnlChiSq(unsigned n, const double* x, double* grad, void* my_func_data)
{
	boost::math::chi_squared_distribution<> chiSqDist1(x[0]);

	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(chiSqDist1, xs[i]));
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;

}