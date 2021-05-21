

#include "weibullDist.h"
#include "boost/math/distributions/normal.hpp" // for normal_distribution
#include <algorithm>
#include <iostream>
#include <fstream>
#include "nlopt.hpp"
#include <iomanip>
// bet is lamb, an, - follow ERA notation
// alp is k
// boost recieves in order (k,an)

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlWeib(unsigned n, const double* x, double* grad, void* my_func_data);
double paramWeibObj(unsigned n, const double* x, double* grad, void* my_func_data);

weibullDist::weibullDist(string opt, vector<double> val, vector<double> add): weibDist(1.,1.)
{
	name = "weibull";
	int npa = 2;
	if (opt.compare("PAR") == 0)
	{
		if (val.size() != npa)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if (!((val[0] > 0) && (val[1] > 0)))
		{
			theErrorFile << "Error running UQ engine: parameters of " << name << " distribution must be greater than 0" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			an = val[0];
			k = val[1];
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
		else if (!((val[0] > 0) && (val[1] > 0))) 
		{
			theErrorFile << "Error running UQ engine: parameters of " << name << " distribution must be greater than 0" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			// need to solve optimization

			// Param optimization
			k = { 0.5 };
			double lb[1] = { 1.e-10 };
			nlopt_opt optimP;
			optimP = nlopt_create(NLOPT_LN_COBYLA, 1); // derivative-free algorithm
			nlopt_set_lower_bounds(optimP, lb);
			nlopt_set_min_objective(optimP, paramWeibObj, &val);
			nlopt_set_xtol_rel(optimP, 1e-6);
			double minf;

			if (nlopt_optimize(optimP, &k, &minf) < 0) {
				printf("nlopt failed!\n");
				theErrorFile << "Error running UQ engine: parameter optimization filed" << std::endl;
				theErrorFile.close();
				exit(1);
			}
			else {
				//lambda = x[0];
				printf("found minimum at f(%g) = %0.10g\n", k, minf);
			}

			an = val[0] / std::tgamma(1 + 1 / k);
		}

	}
	else if (opt.compare("DAT") == 0)
	{

		double minSmp = *std::min_element(std::begin(val), std::end(val));
		if (minSmp < 0)
		{
			theErrorFile << "Error running UQ engine: samples of " << name << " distribution exceeds the range [0,inf]" << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		const int np = 2;
		//double lb[np] = { -HUGE_VAL, 0 };
		//const double  *lb = { 0 };

		// initial lambda
		double mu = 0.0, var = 0.0, ns=0, std;
		for(double value : val)
		{
			mu += value;
			ns ++;
		}
		an = mu/ns/2;
		k = 0.5;

		vector<double> xs = val;

		// MLE optimization
		double lb[2] = { 1.e-10,1.e-10 };
		nlopt_opt optim;
		optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		nlopt_set_min_objective(optim, nnlWeib, &xs);
		nlopt_set_xtol_rel(optim, 1e-6);
		double minf;
		double params[2] = { an,k };

		if (nlopt_optimize(optim, params, &minf) < 0) {
			printf("nlopt failed!\n");
			theErrorFile << "Error running UQ engine: MLE optimization filed" << std::endl;
			theErrorFile.close();
			exit(1);
		}
		else {
			//lambda = x[0];
			printf("found minimum at f(%g) = %0.10g\n", k, minf);
		}

	}

	weibull weibDist1(k,an);
	weibDist = weibDist1;
	checkParams();

}

weibullDist::~weibullDist() {}
//==


void weibullDist::checkParams()
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

	if (mean<=0)
	{
		theErrorFile << "Error running UQ engine: mean of " << name << " distribution must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	if (!(par[0] > 0 && (par[1] > 0)))
	{
		theErrorFile << "Error running UQ engine: parameters of " << name << " distribution must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}
}


double weibullDist::getPdf(double x)
{
	return pdf(weibDist,x);

}

double weibullDist::getCdf(double x)
{
	return cdf(weibDist, x);

}


double weibullDist::getMean(void)
{
	double a =mean(weibDist);
	return mean(weibDist);
}

double weibullDist::getStd(void)
{
	return standard_deviation(weibDist);
}

double weibullDist::getQuantile(double p)
{

	return quantile(weibDist, p);

}

vector<double> weibullDist::getParam(void)
{
	return { an,k };
}


string weibullDist::getName(void)
{
	return name;
}


double nnlWeib(unsigned n, const double* x, double* grad, void* my_func_data)
{
	weibull weibDist(x[1], x[0]);

	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(weibDist, xs[i]));
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}

double paramWeibObj(unsigned n, const double* x, double* grad, void* my_func_data)
{
	my_samples* values = (my_samples*)my_func_data;
	vector<double> val = values->xs;

	if (grad) {
		grad[0] = 0;
	}
	return abs(std::sqrt(tgamma(1 + 2 / x[0]) - (tgamma(1 + 1 / x[0]))* tgamma(1 + 1 / x[0])) / tgamma(1 + 1 / x[0]) - val[1] / val[0]);

}