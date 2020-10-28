

#include "truncExponentialDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "nlopt.hpp"
#include <iomanip>

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlTruncExponential(unsigned n, const double* x, double* grad, void* my_func_data);
double paramTruncExpbObj(unsigned n, const double* x, double* grad, void* my_func_data);

truncExponentialDist::truncExponentialDist(string opt, vector<double> val, vector<double> add)
{
	name = "TruncatedExponential";
	int np=3;
	if (opt.compare("PAR") == 0)
	{
		if (val.size() != np)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= 0) || (val[1]<=0) || (val[1]>val[2]))
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			lambda = val[0];
			a = val[1];
			b = val[2];
		}
	}
	else if (opt.compare("MOM") == 0)
	{

		if (val.size() != np)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ( (val[0] <= 0) || (val[1]<=0) || (val[1]>val[2]) || (2*val[0]>=(val[1]+val[2])) )
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			a = val[1];
			b = val[2];

			lambda = { 1/val[0] };
			double lb[1] = { 1.e-10 };
			nlopt_opt optimP;
			optimP = nlopt_create(NLOPT_LN_COBYLA, 1); // derivative-free algorithm
			nlopt_set_lower_bounds(optimP, lb);
			nlopt_set_min_objective(optimP, paramTruncExpbObj, &val);
			nlopt_set_xtol_rel(optimP, 1e-6);
			double minf;

			if (nlopt_optimize(optimP, &lambda, &minf) < 0) {
				printf("nlopt failed!\n");
				theErrorFile << "Error running UQ engine: parameter optimization filed" << std::endl;
				theErrorFile.close();
				exit(1);
			}
			else {
				//lambda = x[0];
				printf("found minimum at f(%g) = %0.10g\n", lambda, minf);
			}


		}

	}
	else if (opt.compare("DAT") == 0)
	{

		if (add.size() == 0)
		{
			theErrorFile << "Error running UQ engine: Provide a valid range" << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		a = add[0];
		b = add[1];

		const int np = 1;

		// data
		my_samples info;
		info.xs = val;
		info.add = add;
		// initial value
		double mu = 0.0, ns = 0;
		for (double value : val)
		{
			mu += value;
			ns += 1;
		}
		mu = mu / ns;
		lambda = { 1 / mu };
		// MLE optimization
		double lb[1] = { 1.e-10 };
		nlopt_opt optim;
		optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		nlopt_set_min_objective(optim, nnlTruncExponential, &info);
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

	normConst=cdf(expDist,b)-cdf(expDist,a);
}

truncExponentialDist::~truncExponentialDist() {}


double truncExponentialDist::getPdf(double x)
{
	
	double pdfval;
	if ((x<a) || (x>b))
	{
		pdfval = 0;
	} else {
		pdfval=pdf(expDist,x)/normConst;
	}
	return pdfval;

}

double truncExponentialDist::getCdf(double x)
{
	double cdfval;
	if (x>b){
		cdfval=1; 
	}
	else if (x>a)
	{
		cdfval=(cdf(expDist,x)-cdf(expDist,a))/normConst;
	} else {
		cdfval = 0;
	}
	return cdfval;

}


double truncExponentialDist::getMean(void)
{
	double up = std::min(b,getQuantile(0.999));
	int ncount=1.e5;
  	double step = (up - a) / ncount;
  	double integral = 0.5*(a*getPdf(a)+ up* getPdf(up));

  	for (int i=1; i<ncount+1;i++) {
  		double xval = a+step*i;
  		integral += xval* getPdf(xval);
  	}

  	integral *= step;
	return integral;
}

double truncExponentialDist::getStd(void)
{
	double up = std::min(b, getQuantile(0.999));
	int ncount = 1.e5;
	double step = (up - a) / ncount;
	double mean = getMean();
	double integral = 0.5 * ((a* a) * getPdf(a) + (up*up) * getPdf(up));

	for (int i = 1; i < ncount + 1; i++) {
		double xval = a + step * i;
		integral += (xval*xval)* getPdf(xval);
	}

	integral *= step;
	return std::sqrt(integral-mean* mean);
}

double truncExponentialDist::getQuantile(double p)
{

	double minpf = cdf(expDist, a);
	return quantile(expDist, minpf + p*normConst);

}

vector<double> truncExponentialDist::getParam(void)
{
	return { lambda, a, b };
}


string truncExponentialDist::getName(void)
{
	return name;
}


double nnlTruncExponential(unsigned n, const double* x, double* grad, void* my_func_data)
{
	double nll=0;

	exponential expDist1(x[0]);

	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;
	vector<double> add = samp->add;

	double lb = add[0];
	double ub = add[1];

	double normConst = cdf(expDist1, ub) - cdf(expDist1, lb);
	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(expDist1, xs[i]))+log(normConst);
	}

	if (grad) {
		grad[0] = 0;
	}

	return nll;
}



double paramTruncExpbObj(unsigned n, const double* x, double* grad, void* my_func_data)
{
	my_samples* values = (my_samples*)my_func_data;
	vector<double> val = values->xs;
	double mean = val[0];
	double lb = val[1];
	//double ub = val[2];
	exponential expDist(x[0]);

	double ub = std::min(val[2], quantile(expDist,0.999));

	int ncount = 1.e5;
	double step = (ub - lb) / ncount;
	double integral = 0.5 * (lb * pdf(expDist,lb) + ub * pdf(expDist, lb));

	for (int i = 1; i < ncount + 1; i++) {
		double xval = lb + step * i;
		integral += xval * pdf(expDist, xval);
	}

	integral *= step;
	
	// let us use trapzoidal integration to find the mean


	if (grad) {
		grad[0] = 0;
	}
	return abs(integral / (cdf(expDist, ub) - cdf(expDist, lb))  - mean);

}

