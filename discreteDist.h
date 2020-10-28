#pragma once

#include "RVDist.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

class discreteDist : public RVDist
{
public:
	discreteDist(string opt = "", vector<double> val = {}, vector<double> add = {});
	virtual ~discreteDist();
	double getCdf(double x);
	double getPdf(double x);
	double getMean(void);
	double getStd(void);
	double getQuantile(double p);
	string getName(void);
	vector<double> getParam(void);

private:
	vector<double> value;
	vector<double> weight;
	int nv;
	vector<int> sortIdx;

};
