#pragma once

#include "RVDist.h"
#include <string>
#include <vector>
#include "boost/math/distributions/gamma.hpp" 


using boost::math::gamma_distribution;
using std::string;
using std::vector;

class gammaDist : public RVDist
{
public:
	gammaDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~gammaDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	double getQuantile2(double p);

	string getName(void);
	vector<double> getParam(void);

	boost::math::gamma_distribution<> gamDist;

private:
	double k,lambda;
	void checkParams();
};
