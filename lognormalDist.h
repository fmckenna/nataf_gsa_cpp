#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/lognormal.hpp" // for normal_distribution


using boost::math::lognormal_distribution;
using std::string;
using std::vector;

class lognormalDist : public RVDist
{
public:
	lognormalDist(string opt = "", vector<double> val = {}, vector<double> add = {0,1});
	virtual ~lognormalDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	boost::math::lognormal_distribution<> lognDist;

private:
	double lambda, zeta;

};
