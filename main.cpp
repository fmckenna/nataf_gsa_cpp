#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>

#include "jsonInput.h"
#include "ERADist.h"
#include "ERANataf.h"
#include "runGSA.h"


std::ofstream theErrorFile; // Error log

void writeOutputs(jsonInput inp, vector<vector<double>> xvals, vector<vector<double>> gvals, runGSA GsaResults);

int main(int argc, char** argv)
//int main(void)
{ 
	std::string workDir = argv[1];
	std::string osType = argv[2];
	std::string runType = argv[3];

	//std::string workDir = "C:/Users/yisan/Documents/quoFEM/LocalWorkDir/tmp.SimCenter";
	//std::string osType = "Windows";
	//std::string runType = "runningLocal";

	std::cerr << "WORK " << workDir << "\n";
	std::cerr << "OS " << osType << "\n";
	std::cerr << "RUN " << runType << "\n";

	theErrorFile.open(workDir+"/dakota.err");
	
	//
	//  (1) read JSON file
	//
	
	jsonInput inp(workDir);

	//
	//	(2) Construct Nataf Object
	//

	ERANataf T(inp);

	//
	//	(3-1) Random number generator Gaussian(mean=0,var=1) - (batch samples)
	//
	
	int nmc = inp.nmc;
	int nreg = inp.nreg;

	std::mt19937 generator(inp.rseed);
	std::normal_distribution<double> distribution(0.0, 1.0);

	vector<vector<double>> uvals(nmc, vector<double>(inp.nrv, 0.0));
	for (int ns = 0; ns < nmc; ns++)
	{
		for (int nr = 0; nr < inp.nrv; nr++)
			uvals[ns][nr] = distribution(generator);
	}

	std::cout<<std::to_string(inp.nreg)<<std::endl;
	std::cout<<std::to_string(nreg);

	vector<vector<int>> resampIDs(nmc, vector<int>(nreg, 0.0));
	for (int nr = 0; nr < nreg; nr++)
	{
		std::uniform_int_distribution<int> discrete_dist(0, inp.resamplingSize[nr]-1);
		for (int ns = 0; ns < nmc; ns++)
		{
			resampIDs[ns][nr] = discrete_dist(generator);
			std::cout << resampIDs[ns][nr] << std::endl;
		}
	}

	
	//
	//	(4-1) FE Analysis - (batch samples)
	//
	
	vector<vector<double>> gvals, xvals;
	T.simulateAppBatch(osType, runType, inp, uvals, resampIDs, xvals, gvals);


	//
	//	(3-2) and (4-2) alternative (sequential)
	//
/*
	int nmc = inp.nmc;
	std::mt19937 generator(inp.rseed);
	std::normal_distribution<double> distribution(0.0, 1.0);

	vector<vector<double>> gvals, xvals;
	for (int ns = 0; ns < nmc; ns++)
	{
		vector<double> uval;
		for (int nr = 0; nr < inp.nrv; nr++)
			uval.push_back(distribution(generator));

		T.simulateAppSequential(osType, runType, inp, { uval }, xvals, gvals, ns);
		//std::cout<< gvals[ns][0] << " " << gvals[ns][1] << std::endl;
	}	
*/
	//
	//	(5) Global sensitivity analysis
	//

	int Kos = 25;
	runGSA GsaResults(xvals, gvals, inp.groups, Kos); 

	//
	//	Write files/dakota.out, dakotaTab.out
	//

	GsaResults.writeOutputs(inp);
	
	return 0;

}

