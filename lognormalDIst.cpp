

#include "lognormalDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "nlopt.hpp"

#include <iomanip>

using std::vector;

extern std::ofstream theErrorFile; // Error log
//double nnlLogn(unsigned n, const double* x, double* grad, void* my_func_data);

lognormalDist::lognormalDist(string opt, vector<double> val, vector<double> add) 
{
	name = "lognormal";
	int npa = 2;
	if (opt.compare("PAR") == 0)
	{

		if (val.size()!= npa)
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
			lambda = val[0]; // log mean
			zeta = val[1]; // log std
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
		else if ((val[0] <= 0) || (val[1] <= 0))
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			double mean = val[0];
			double std = val[1];

			lambda = log(mean) - log(sqrt(1 + (std / mean)* (std / mean))); // mean normal
			zeta = sqrt(log(1 + (std / mean)* (std / mean))); // sigma normal

		}

	}
	else if (opt.compare("DAT") == 0)
	{
		{
			const int np = 1;

			lambda = 0.0;
			double var = 0.0, ns = 0;
			for (double value : val)
			{
				double lnVal = log(value);
				lambda += lnVal;
				var += lnVal* lnVal;
				ns += 1;
			}
			lambda = lambda / ns;
			zeta = sqrt(var / ns - lambda * lambda);
		}

	}
	boost::math::lognormal_distribution<> lognDist1(lambda, zeta); // shape k, scale theta=1/lambda (or 1/lognormal)
	lognDist = lognDist1;

}

lognormalDist::~lognormalDist() {}

//==

double lognormalDist::getPdf(double x)
{
	return pdf(lognDist,x);
}

double lognormalDist::getCdf(double x)
{
	return cdf(lognDist, x);
}

double lognormalDist::getMean(void)
{

	//std::cout << (logn * a + alp * b) / (alp + logn) << std::endl;
	//std::cout << mean(lognDist)*(b-a)+a << std::endl;
	return mean(lognDist);
}

double lognormalDist::getStd(void)
{
	return standard_deviation(lognDist) ;
}

double lognormalDist::getQuantile(double p)
{
	return quantile(lognDist, p) ;
}

string lognormalDist::getName(void)
{
	return name;
}


vector<double> lognormalDist::getParam(void)
{
	return { lambda, zeta };
}

/*
double nnlLogn(unsigned n, const double* x, double* grad, void* my_func_data)
{

	boost::math::lognormal_distribution<> lognDist1(x[0], x[1]);
	
	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(lognDist1,(xs[i])));
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}
*/

