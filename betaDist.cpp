

#include "betaDist.h"
#include "nlopt.hpp"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlBeta(unsigned n, const double* x, double* grad, void* my_func_data);

betaDist::betaDist(string opt, vector<double> val, vector<double> add) : betDist(1.0,1.0)
{
	name = "beta";
	int npa = 4;
	if (opt.compare("PAR") == 0)
	{

		if (val.size()!= npa)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= 0) || (val[1] <= 0) || (val[2] >= val[3]))
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			alp = val[0];
			bet = val[1];
			a = val[2];
			b = val[3];
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
		else if (val[1] <= 0)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= val[2]) || (val[0] >= val[3]))
		{
			theErrorFile << "Error running UQ engine: mean of " << name << " distribution must be in a valid range " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			double mu = val[0];
			double sig = val[1];
			a = val[2];
			b = val[3];

			alp = ((b - mu) * (mu - a) / (sig*sig) - 1) * (mu - a) / (b - a);
			bet = alp * (b - mu) / (mu - a);
		}
		if (!((alp > 0) && (bet > 0)))
		{
			theErrorFile << "Error running UQ engine: parameters of " << name << " distribution must be greater than 0 " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
	}
	else if (opt.compare("DAT") == 0)
	{

		if (add.size() == 0)
		{
			theErrorFile << "Error running UQ engine: Provide a valid rangef for " << name << " distribution." << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		a = add[0];
		b = add[1];

		double maxSmp = *std::max_element(std::begin(val), std::end(val));
		double minSmp = *std::min_element(std::begin(val), std::end(val));
		if ((maxSmp > b) || (minSmp < a))
		{
			theErrorFile << "Error running UQ engine: samples of " << name << " distribution exceeds the range [min,max]" << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		const int np = 2;

		double mu = 0.0;
		double var = 0.0, ns = 0;
		for (double value : val)
		{
			mu += value;
			var += value * value;
			ns += 1;
		}
		mu = mu / ns;
		double sig = sqrt(var / ns - mu * mu);
		//modify
		alp = ((b - mu) * (mu - a) / (sig * sig) - 1) * (mu - a) / (b - a);
		bet = mu * (b - mu) / (mu - a);

		vector<double> xs;
		for (int i = 0; i < xs.size(); i++)
		{
			xs.push_back((xs[i] - a) / (b - a));
		}

		
		double params[2] = { alp, bet };

		// MLE optimization
		double lb[2] = { 1.e-10 ,  1.e-10 };
		nlopt_opt optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		
		nlopt_set_min_objective(optim, nnlBeta, &xs);
		nlopt_set_xtol_rel(optim, 1e-6);
		double minf; // `*`the` `minimum` `objective` `value,` `upon` `return`*` 

		//double params;
		if (nlopt_optimize(optim, params, &minf) < 0) {
			printf("nlopt failed!\n");
			theErrorFile << "Error running UQ engine: MLE optimization filed" << std::endl;
			theErrorFile.close();
			exit(1);
		}
		else {
			//lambda = x[0];
			printf("found minimum at f(%g) = %0.10g\n", params[0], minf);
		}

		alp = params[0];
		bet = params[1];

	}
	boost::math::beta_distribution<> betDist1(alp, bet); // shape k, scale theta=1/lambda (or 1/beta)
	betDist = betDist1;
	checkParams();
}


betaDist::~betaDist() {}

void betaDist::checkParams()
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
	
	if ((mean <= par[2]) || (mean >=par[3]))
	{
		theErrorFile << "Error running UQ engine: mean of " << name << " distribution must be in a valid range " << std::endl;
		theErrorFile.close();
		exit(-1);
	}
	
	if (par[2] > par[3])
	{ 
		theErrorFile << "Error running UQ engine: range of " << name << "distribution is not valid " << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	if (!(par[0]>0 && (par[1]>0)))
	{
		theErrorFile << "Error running UQ engine: parameters of " << name << " distribution must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}
}


//==

double betaDist::getPdf(double x)
{
	return pdf(betDist, (x - a) / (b - a)) / (b - a);
}

double betaDist::getCdf(double x)
{
	return cdf(betDist, (x - a) / (b - a));
}

double betaDist::getMean(void)
{

	//std::cout << (bet * a + alp * b) / (alp + bet) << std::endl;
	//std::cout << mean(betDist)*(b-a)+a << std::endl;
	return mean(betDist) * (b - a) + a;
}

double betaDist::getStd(void)
{
	return standard_deviation(betDist) * (b - a);
}

double betaDist::getQuantile(double p)
{
	return quantile(betDist, p) * (b - a) + a;
}

vector<double> betaDist::getParam(void)
{
	return { alp,bet,a,b };
}

string betaDist::getName(void)
{
	return name;
}

double nnlBeta(unsigned n, const double* x, double* grad, void* my_func_data)
{
	//k:x[0]
	//lambda:x[1]
	//theta: 1/x[1]
	//Bet(k,theta)
	boost::math::beta_distribution<> BetDist1(x[0], x[1]);
	
	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(BetDist1,(xs[i])));
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}


