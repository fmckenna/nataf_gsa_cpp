#pragma once

#include "RVDist.h"
#include "boost/math/distributions/normal.hpp" // for normal_distribution

#include <string>
#include <vector>

using boost::math::normal;
using std::string;
using std::vector;

class normalDist : public RVDist
{
public:
	normalDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~normalDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	normal normDist;

private:
	double mu;
	double sigma;

};
