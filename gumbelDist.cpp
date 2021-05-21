

#include "gumbelDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "nlopt.hpp"

#include <iomanip>

using std::vector;

extern std::ofstream theErrorFile; // Error log
double nnlGumb(unsigned n, const double* x, double* grad, void* my_func_data);

gumbelDist::gumbelDist(string opt, vector<double> val, vector<double> add) 
{
	name = "gumbel";
	int npa = 2;
	if (opt.compare("PAR") == 0)
	{

		if (val.size()!= npa)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else if ((val[0] <= 0))
		{
			theErrorFile << "Error running UQ engine: scale parameter of " << name << " distribution must be greater than 0 " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			alp = val[0]; // 1/an
			bet = val[1];
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
			theErrorFile << "Error running UQ engine: stdandard deviation of " << name << " distribution must be greater than 0 " << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			double an = val[1] * sqrt(6) / PI;// ; // scale parameter
			bet = val[0] - EC * an; // location parameter, EC: euler constant
			alp = 1 / an;
		}

	}
	else if (opt.compare("DAT") == 0)
	{
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
		double sigma = sqrt(var / ns - mu * mu);

		double an = sigma * sqrt(6) / PI;// ; // scale parameter
		bet = mu - EC * an; // location parameter, EC: euler constant
		alp = 1 / an;


		vector<double> xs = val;
		double params[2] = { alp,bet };


		// MLE optimization
		double lb[2] = { 1.e-10 , -INFINITY };
		nlopt_opt optim = nlopt_create(NLOPT_LN_COBYLA, np); // derivative-free algorithm
		nlopt_set_lower_bounds(optim, lb);
		
		nlopt_set_min_objective(optim, nnlGumb, &xs);
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

		alp = params[0]; // 1/scale
		bet = params[1]; // loca


	}

	boost::math::extreme_value_distribution<> gumbDist1(bet, 1 / alp);
	gumbDist = gumbDist1;

	// (location, scale) parameters
	// =(a, b) in boost-doc
	// =(bn, an) in ERADist
	// =(bet, 1/alp) Dakota manual
	// User input of quoFEM & Dakota : alp, bet
	checkParams();
}

gumbelDist::~gumbelDist() {}

//==

void gumbelDist::checkParams()
{
	double std = getStd();
	double mean = getMean();
	vector<double> par = getParam();

	if (isnan(std) || isinf(std) || std <= 0)
	{
		theErrorFile << "Error running UQ engine: stdandard deviation of " << name << " distribution must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	if (par[0] <= 0)
	{
		theErrorFile << "Error running UQ engine: scale parameter of " << name << " distribution must be greater than 0 " << std::endl;
		theErrorFile.close();
		exit(-1);
	}
}


double gumbelDist::getPdf(double x)
{
	return pdf(gumbDist, x);
}

double gumbelDist::getCdf(double x)
{
	return cdf(gumbDist, x);
}

double gumbelDist::getMean(void)
{
	//std::cout << bet + EC / alp << std::endl;
	//std::cout << mean(gumbDist) << std::endl;
	return mean(gumbDist);
}

double gumbelDist::getStd(void)
{
	//std::cout << PI / alp / sqrt(6) << std::endl;
	//std::cout << standard_deviation(gumbDist) << std::endl;
	return standard_deviation(gumbDist);
}

double gumbelDist::getQuantile(double p)
{

	//printf("BOOST: %.5f\n", bet - 1 / alp*log(-log(p)));
	//printf("BOOST: %.5f\n", quantile(gumbDist, p));

	return quantile(gumbDist, p); // slow
	//return  bet - 1 / alp * log(-log(p));
}

vector<double> gumbelDist::getParam(void)
{
	return { alp, bet};
}

string gumbelDist::getName(void)
{
	return name;
}


double nnlGumb(unsigned n, const double* x, double* grad, void* my_func_data)
{
	boost::math::extreme_value_distribution<> gumbDist1(x[1], 1 / x[0]);
	
	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	double loc = x[1];
	double sca = 1/x[0];

	for (int i = 0; i < xs.size(); i++)
	{
		//nll += -log(pdf(gumbDist1,xs[i]));
		nll += log(sca)+(xs[i]- loc)/sca + exp(-(xs[i] - loc) / sca);
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;

}


