

#include "discreteDist.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlopt.hpp"



using std::vector;

extern std::ofstream theErrorFile; // Error log
//double nnldiscrete(unsigned n, const double* x, double* grad, void* my_func_data);

discreteDist::discreteDist(string opt, vector<double> val, vector<double> add )
{
	name = "discrete";
	int np = val.size();
	if ((opt.compare("PAR") == 0))
	{
		if (np%2==1)
		{
			theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		else
		{
			nv = np / 2;
			value.reserve(nv);
			weight.reserve(nv);
			double sumWei = 0;
			for (int i = 0 ; i < nv; i++)
			{
				double w = val[i * 2 + 1];
				if (w <0)
				{
					theErrorFile << "Error running UQ engine: The " << name << " distribution is not defined for your parameters" << std::endl;
					theErrorFile.close();
					exit(-1);
				}
				value.push_back(val[i * 2]);
				weight.push_back(w);
				sumWei += w;
			}
			if (sumWei/=1)
			{
				for (int i = 0; i < np / 2; i++)
				{
					weight[i]= weight[i] / sumWei;
				}
			}
		}
	}
	else if (opt.compare("MOM") == 0)
	{
		theErrorFile << "Error running UQ engine: Moment input is not supported for " << name << " distribution" << std::endl;
		theErrorFile.close();
		exit(-1);
	}
	else if (opt.compare("DAT") == 0)
	{
		value = val;
		value.erase(unique(value.begin(), value.end()), value.end());
		nv = value.size();
		int nx = val.size();
		weight.reserve(nv);
		for (int i = 0; i < nv; i++)
		{
			int count = 0;
			for (double v : val)
			{
				if (v == value[i])
				{
 					count += 1;
				}
			}
			weight.push_back((double) count / nx);
		}
	}

	vector<int> y(value.size());
	std::size_t n(0);
	std::generate(std::begin(y), std::end(y), [&] { return n++; });

	std::sort(std::begin(y),
			  std::end(y),
			  [&](int i1, int i2) { return value[i1] < value[i2]; });

	sortIdx = y;
}

discreteDist::~discreteDist() {}

//==

double discreteDist::getPdf(double x)
{

	//return (1.0 / (sigma * sqrt(2.0 * PI))) * exp(-0.5 * (x - mu)*(x - mu) / (sigma * sigma));
	double pdfval = 0;
	for (int i=0;i<nv;i++)
	{
		if (value[i] == x)
		{
			pdfval = weight[i];
		}
	}
	return pdfval;
}

double discreteDist::getCdf(double x)
{
	double cdfval = 0;
	//return erfc(-x / sqrt(2.0)) / 2.0;
	for (int i = 0; i < nv; i++)
	{
		if (value[i] < x)
		{
			cdfval += weight[i];
		}
	}
	return cdfval;
}

double discreteDist::getMean(void)
{
	double meanval = 0;
	for (int i = 0; i < nv; i++)
	{
		meanval += weight[i] * value[i];
	}
	return meanval;
}

double discreteDist::getStd(void)
{
	double varval = 0;
	for (int i = 0; i < nv; i++)
	{
		varval += weight[i] * value[i] * value[i];
	}
	double meanval = getMean();
	return std::sqrt(varval - meanval * meanval);
}

double discreteDist::getQuantile(double p)
{
	double icdfval= HUGE_VAL;
	double cumwei=0;
	for (int i : this->sortIdx)
	{
		cumwei += weight[i];
		if (cumwei>p)
		{ 
			icdfval = value[i];
			break;
		}
	}
	return icdfval;
}

string discreteDist::getName(void)
{
	return name;
}

vector<double> discreteDist::getParam(void)
{
	vector<double> param;
	param.reserve(nv * 2);
	for (int i=0; i < nv; i++)
	{
		param.push_back(value[i]);
		param.push_back(weight[i]);
	}
	return param;
}


/*
double nnldiscrete(unsigned n, const double* x, double* grad, void* my_func_data)
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


