#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/weibull.hpp" // for normal_distribution


using boost::math::weibull;
using std::string;
using std::vector;

class weibullDist : public RVDist
{
public:
	weibullDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~weibullDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);
	weibull weibDist;

private:
	double an,k;
	void checkParams(void);
	
};
