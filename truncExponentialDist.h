#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/exponential.hpp" // for normal_distribution

using boost::math::exponential;
using std::string;
using std::vector;

class truncExponentialDist : public RVDist
{
public:
	truncExponentialDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~truncExponentialDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	exponential expDist;

private:
	double lambda,a,b;
	double normConst;
	void checkParams();
	
};
