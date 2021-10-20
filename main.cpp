#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <chrono>

#include "jsonInput.h"
#include "ERADist.h"
#include "ERANataf.h"
#include "runGSA.h"
#include "runForward.h"


std::ofstream theErrorFile; // Error log

void writeOutputs(jsonInput inp, vector<vector<double>> xvals, vector<vector<double>> gvals, runGSA GsaResults);

int main(int argc, char** argv)
//int main(void)
{ 
	std::string errMsg;

	if (argc != 4) {

		errMsg = "Number of the addition commend line arguments is " + std::to_string(argc-1) + ", but 3 is required. The arguments should always include the working directory / os type / run type";
		std::cout << errMsg << "\n";
		// theErrorFile << errMsg << std::endl;
		// theErrorFile.close();
		// exit(-1);
	}

	std::string workDir = argv[1];
	std::string osType = argv[2];
	std::string runType = argv[3];




	//std::string workDir = "C:/Users/yisan/Documents/quoFEM/LocalWorkDir/tmp.SimCenter";
	//std::string osType = "Windows";
	//std::string runType = "runningLocal";

	std::cout << "WORKDIR: " << workDir << "\n";
	std::cout << "OS: " << osType << "\n";
	std::cout << "RUN: " << runType << "\n";
	std::cout << "\n";
	
	auto start = std::chrono::high_resolution_clock::now();


	theErrorFile.open(workDir+"/dakota.err",std::ofstream::out);
	if (!theErrorFile.is_open()) {
		errMsg = "Error running UQ engine: Failed to creat Dakota.err at " + workDir;
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}
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

	vector<vector<int>> resampIDs(nmc, vector<int>(nreg, 0.0));
	for (int nr = 0; nr < nreg; nr++) // if {0,1},{2,3}, nrg is 2
	{
		if (nmc > inp.resamplingSize[nr])	{ 
			std::uniform_int_distribution<int> discrete_dist(0, inp.resamplingSize[nr] - 1);
			for (int ns = 0; ns < nmc; ns++)
			{
				resampIDs[ns][nr] = discrete_dist(generator);
				std::cout << resampIDs[ns][nr] << std::endl;
			}
		} else {
			std::vector<int> v(nmc); // vector with nmc ints.
			std::iota(std::begin(v), std::end(v), 0); // Fill with 0, 1, ..., 99.
			std::shuffle(v.begin(), v.end(), generator);
			for (int ns = 0; ns < nmc; ns++)
			{
				resampIDs[ns][nr] = v[ns];
				std::cout << resampIDs[ns][nr] << std::endl;
			}

		}
	}
	
	//
	//	(4-1) FE Analysis - (batch samples)
	//
	vector<vector<double>> gvals;
	vector<vector<double>> xvals;
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
	runForward ForwardResults(xvals, gvals);
	//
	//	Write dakotaTab.out
	//
	ForwardResults.writeTabOutputs(inp);

	if (!inp.uqType.compare("Sensitivity Analysis")) {
		//
		//	(5) Global sensitivity analysis
		//
		int Kos = 25;
		runGSA GsaResults(xvals, gvals, inp.groups, Kos);
		//
		//	Write  dakota.out
		//
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time: " << elapsed.count() << " s\n";

		GsaResults.writeOutputs(inp, elapsed.count());
	}
	else if (!inp.uqType.compare("Forward Propagation")) {

		//
		//	(5) Forward analysis
		//
		//
		//	Write  dakota.out
		//
		ForwardResults.writeOutputs(inp);

	}

	return 0;

}

