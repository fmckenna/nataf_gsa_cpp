#pragma once

#include "RVDist.h"

#include <string>
#include <vector>


#include "boost/math/distributions/beta.hpp" // for normal_distribution


using boost::math::beta_distribution;
using std::string;
using std::vector;

class betaDist : public RVDist
{
public:
	betaDist(string opt = "", vector<double> val = {}, vector<double> add = {0,1});
	virtual ~betaDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

	boost::math::beta_distribution<> betDist;

private:
	double a,b,alp,bet;

};
