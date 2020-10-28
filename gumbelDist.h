#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/extreme_value.hpp" // for normal_distribution


using boost::math::extreme_value_distribution;
using std::string;
using std::vector;

class gumbelDist : public RVDist
{
public:
	gumbelDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~gumbelDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	boost::math::extreme_value_distribution<> gumbDist;

private:
	double alp,bet; // NOTE: alpha=1/an, beta=bn in ERA
	//double EC = 0.57721566490153; // euler constant
	double EC = 0.57721566490153286060651209008240243104215933593992;
};
