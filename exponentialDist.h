#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/exponential.hpp" // for normal_distribution


using boost::math::exponential;
using std::string;
using std::vector;

class exponentialDist : public RVDist
{
public:
	exponentialDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~exponentialDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	exponential expDist;

private:
	double lambda;
	void checkParams();
};
