#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/uniform.hpp" // for normal_distribution


using boost::math::uniform_distribution;
using std::string;
using std::vector;

class uniformDist : public RVDist
{
public:
	uniformDist(string opt = "", vector<double> val = {}, vector<double> add = {0,1});
	virtual ~uniformDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	boost::math::uniform_distribution<> unifDist;

private:
	double a,b;


};
