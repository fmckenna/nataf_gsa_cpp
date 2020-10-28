

#include "exponentialDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlopt.hpp"

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlExponential(unsigned n, const double* x, double* grad, void* my_func_data);

exponentialDist::exponentialDist(string opt, vector<double> val, vector<double> add)
{
	name = "Exponential";
	if (opt.compare("PAR") == 0)
	{
		if (val.size() != 1)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if (val[0] <= 0)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			lambda = val[0];
		}
	}
	else if (opt.compare("MOM") == 0)
	{

		if (val.size() != 1)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if (val[0] <= 0)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			lambda = 1/val[0];
		}

	}
	else if (opt.compare("DAT") == 0)
	{
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
		mu = mu / ns;
		lambda = 1 / mu;
		vector<double> xs = val;

		// MLE optimization
		double lb[1] = { 1.e-10 };
		//double x[1] = { lambda };  //`*`some` `initial` `guess`*` 
		nlopt_opt optim;
		//optim = nlopt_create(NLOPT_LD_MMA, np); // gradient-based algorithm
		optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		//vector<double> xs = { 10.0, 20.0, 30.0 };
		nlopt_set_min_objective(optim, nnlExponential, &xs);
		nlopt_set_xtol_rel(optim, 1e-6);
		double minf; // `*`the` `minimum` `objective` `value,` `upon` `return`*` 

		if (nlopt_optimize(optim, &lambda, &minf) < 0) {
			printf("nlopt failed!\n");
			theErrorFile << "Error running UQ engine: MLE optimization filed" << std::endl;
			theErrorFile.close();
			exit(1);
		}
		else {
			//lambda = x[0];
			printf("found minimum at f(%g) = %0.10g\n", lambda, minf);
		}

	}

	exponential expDist1(lambda);
	expDist = expDist1;

}

exponentialDist::~exponentialDist() {}
//==


double exponentialDist::getPdf(double x)
{
	//double result;
	//if (0.0 < x) {
	//	result = lambda * exp(-lambda * x);
	//}
	//else {
	//	result = 0.0;
	//}
	//return result;
	return pdf(expDist,x);

}

double exponentialDist::getCdf(double x)
{
	/*
	double result;
	if (0.0 < x) {
		result = 1.0-exp(-(lambda) * x);
	}
	else {
		result = 0.0;
	}
	return result;
	*/
	return cdf(expDist, x);

}


double exponentialDist::getMean(void)
{
	//return moms[0];
	//double a =mean(expDist);
	return mean(expDist);
}

double exponentialDist::getStd(void)
{
	//return moms[1];
	return standard_deviation(expDist);
}

double exponentialDist::getQuantile(double p)
{

	return quantile(expDist, p);

}

vector<double> exponentialDist::getParam(void)
{
	return { lambda };
}


string exponentialDist::getName(void)
{
	return name;
}


double nnlExponential(unsigned n, const double* x, double* grad, void* my_func_data)
{
	double lamb = x[0];
	double nll=0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for(int i=0;i< xs.size();i++)
	{
		nll += - log(lamb) - (-lamb * xs[i]);
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}



