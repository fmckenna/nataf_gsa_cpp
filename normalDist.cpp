

#include "normalDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlopt.hpp"


using std::vector;

extern std::ofstream theErrorFile; // Error log
//double nnlNormal(unsigned n, const double* x, double* grad, void* my_func_data);

normalDist::normalDist(string opt, vector<double> val, vector<double> add )
{
	name = "normal";
	if ((opt.compare("PAR") == 0) || (opt.compare("MOM") == 0))
	{

		if (val.size() != 2)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[1] <= 0))
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			mu = val[0];
			sigma = val[1];
		}
	}
	else if (opt.compare("DAT") == 0)
	{
		const int np = 1;

		mu = 0.0;
		double var = 0.0, ns = 0;
		for (double value : val)
		{
			mu += value;
			var += value * value;
			ns += 1;
		}
		mu = mu / ns;
		sigma = sqrt(var / ns - mu * mu);
	}

	normal normDist1(mu,sigma);
	normDist = normDist1;
}

normalDist::~normalDist() {}

//==

double normalDist::getPdf(double x)
{

	//return (1.0 / (sigma * sqrt(2.0 * PI))) * exp(-0.5 * (x - mu)*(x - mu) / (sigma * sigma));
	return pdf(normDist, x);
}

double normalDist::getCdf(double x)
{
	//return erfc(-x / sqrt(2.0)) / 2.0;
	return cdf(normDist, x);
}

double normalDist::getMean(void)
{
	return mean(normDist);
}

double normalDist::getStd(void)
{
	return standard_deviation(normDist);
}

double normalDist::getQuantile(double p)
{
	return quantile(normDist, p);
}

string normalDist::getName(void)
{
	return name;
}

vector<double> normalDist::getParam(void)
{
	return { mu, sigma };
}


/*
double nnlNormal(unsigned n, const double* x, double* grad, void* my_func_data)
{
	double lamb = x[0];
	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(lamb) - (-lamb * xs[i]);
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}
*/


