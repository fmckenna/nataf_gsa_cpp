
#include "ERADist.h"
#include <algorithm>
#include <iostream>
#include <fstream>

extern std::ofstream theErrorFile; // Error log

ERADist::ERADist() {}

ERADist::ERADist(string getName, string getOpt, vector<double> getVal, vector<double> getAdd)
{
	// NAME
	transform(getName.begin(), getName.end(), getName.begin(), ::tolower);
	this->name = getName;
	// OPT
	transform(getOpt.begin(), getOpt.end(), getOpt.begin(), ::toupper);
	this->opt = getOpt;

	if (getName.compare("exponential") == 0)
	{
		this->theDist = new exponentialDist(getOpt, getVal);
	} 
	else if (getName.compare("normal") == 0)
	{
		this->theDist = new normalDist(getOpt, getVal);
	}
	else if (getName.compare("gamma") == 0)
	{
		this->theDist = new gammaDist(getOpt, getVal);
	}
	else if (getName.compare("beta") == 0)
	{
		this->theDist = new betaDist(getOpt, getVal, getAdd);
	}
	else if (getName.compare("lognormal") == 0)
	{
		this->theDist = new lognormalDist(getOpt, getVal);
	}
	else if (getName.compare("uniform") == 0)
	{
		this->theDist = new uniformDist(getOpt, getVal);
	}
	else if (getName.compare("chisquared") == 0)
	{
		this->theDist = new chiSquaredDist(getOpt, getVal);
	}
	else if (getName.compare("chisquare") == 0)
	{
		this->theDist = new chiSquaredDist(getOpt, getVal);
	}
	else if (getName.compare("gumbel") == 0)
	{
		this->theDist = new gumbelDist(getOpt, getVal);
	}
	else if (getName.compare("weibull") == 0)
	{
		this->theDist = new weibullDist(getOpt, getVal);
	}
	else if (getName.compare("truncatedexponential") == 0)
	{
		this->theDist = new truncExponentialDist(getOpt, getVal, getAdd);
	} 
	else if (getName.compare("discrete") == 0)
	{
		this->theDist = new discreteDist(getOpt, getVal);
	}
	else
	{
		theErrorFile << "Error running UQ engine: " << name << " distribution is not supported" << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	std::cout << getName<< "[" << getOpt << "] =================================" << std::endl;

	printf("Mean: %3.2f, Std: %3.2f \n", theDist->getMean(), theDist->getStd());
	printf("Params: ");
	for (auto par : theDist->getParam())
		printf("%3.2f, ", par);
	printf("\n");
	printf("PDF at mean : %3.2f\n\n", theDist->getPdf(theDist->getMean()));

	/*
	std::cout << "Mean: " << theDist->getMean()  << " Std: " <<  theDist->getStd() << std::endl;
	std::cout << "Params: ";
	for (auto i : theDist->getParam())
		std::cout << i << ', ';
	std::cout << std::endl;
	std::cout << theDist->getPdf(theDist->getMean()) << std::endl;
	*/

}

ERADist::~ERADist(void) {}
