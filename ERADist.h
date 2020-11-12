#pragma once

#include "exponentialDist.h"
#include "normalDist.h"
#include "gammaDist.h"
#include "betaDist.h"
#include "lognormalDist.h"
#include "uniformDist.h"
#include "chiSquaredDist.h"
#include "gumbelDist.h"
#include "weibullDist.h"
#include "truncExponentialDist.h"
#include "discreteDist.h"

#include <string>
#include <vector>
#include "RVDist.h"
using std::string;
using std::vector;

class ERADist
{

public:
	ERADist(void);
	ERADist(string name, string opt, vector<double> val, vector<double> add = {});
	~ERADist();
	RVDist* theDist;
	string name;
	string opt;
	vector<double> par;



private:

};
