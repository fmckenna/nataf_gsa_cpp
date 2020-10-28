

#include "uniformDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
//#include "nlopt.hpp"


#include <iomanip>
//#include <cmath>

using std::vector;

extern std::ofstream theErrorFile; // Error log
//double nnlUnif(unsigned n, const double* x, double* grad, void* my_func_data);

uniformDist::uniformDist(string opt, vector<double> val, vector<double> add) 
{
	name = "uniform";
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
			a = val[0]; // log mean
			b = val[1]; // log std
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


			a = mean - sqrt(12) * std / 2;
			b = mean + sqrt(12) * std / 2;
		}

	}
	else if (opt.compare("DAT") == 0)
	{
		int ns = val.size();
		auto min_value = *std::min_element(val.begin(), val.end());
		auto max_value = *std::max_element(val.begin(), val.end());

		// Estimating the range
		a = max_value + ((double) ns + 1.0) * (min_value - max_value) / ns;
		b = min_value + ((double) ns + 1.0) * (max_value - min_value) / ns;

	}
	boost::math::uniform_distribution<> unifDist1(a, b); // shape k, scale theta=1/lambda (or 1/uniform)
	unifDist = unifDist1;

}

uniformDist::~uniformDist() {}

//==

double uniformDist::getPdf(double x)
{
	return pdf(unifDist,x);
}

double uniformDist::getCdf(double x)
{
	return cdf(unifDist, x);
}

double uniformDist::getMean(void)
{

	//std::cout << (unif * a + alp * b) / (alp + unif) << std::endl;
	//std::cout << mean(unifDist)*(b-a)+a << std::endl;
	return mean(unifDist);
}

double uniformDist::getStd(void)
{
	return standard_deviation(unifDist) ;
}

double uniformDist::getQuantile(double p)
{
	return quantile(unifDist, p) ;
}

string uniformDist::getName(void)
{
	return name;
}


vector<double> uniformDist::getParam(void)
{
	return { a, b };
}

/*
double nnlUnif(unsigned n, const double* x, double* grad, void* my_func_data)
{

	boost::math::uniform_distribution<> unifDist1(x[0], x[1]);
	
	double nll = 0;
	my_samples* samp = (my_samples*)my_func_data;
	vector<double> xs = samp->xs;

	for (int i = 0; i < xs.size(); i++)
	{
		nll += -log(pdf(unifDist1,(xs[i])));
	}
	if (grad) {
		grad[0] = 0;
	}
	return nll;
}

*/
