#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/chi_squared.hpp" // for normal_distribution


using boost::math::chi_squared;
using std::string;
using std::vector;

class chiSquaredDist : public RVDist
{
public:
	chiSquaredDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~chiSquaredDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	chi_squared chiSqDist;

private:
	double k;
	
};
