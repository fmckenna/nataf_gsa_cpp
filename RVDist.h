#pragma once

#include <string>
#include <vector>
using std::string;
using std::vector;

class RVDist
{

public:
	//explicit RVDist(string opt = "", vector<double> val = {});
	RVDist(void);
	virtual ~RVDist(void);
	
	virtual double getPdf(double x) = 0;
	virtual double getCdf(double x) = 0;
	virtual double getMean(void)=0;
	virtual double getStd(void) = 0;
	virtual double getQuantile(double p) = 0;
	virtual string getName(void) = 0;
	virtual vector<double> getParam(void) = 0;

	string name;
	const double PI = 3.141592653589793238462643383279502884197169399375;
private:


};

// For MLE optimization
typedef struct {
	vector<double> xs;
	vector<double> add = {};
} my_samples;

