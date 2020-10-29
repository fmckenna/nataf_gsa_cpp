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

//int main(int argc, char** argv)
int main(void)
{ 
	//std::string workdir = argv[1];
	//std::string osType = argv[2];
	//std::string runType = argv[3];

	std::string workDir = "C:/Users/yisan/Documents/quoFEM/LocalWorkDir/tmp.SimCenter";
	std::string osType = "Windows";
	std::string runType = "runningLocal";

	std::cerr << "WORK " << workDir << "\n";
	std::cerr << "OS " << osType << "\n";
	std::cerr << "RUN " << runType << "\n";

	theErrorFile.open(workDir+"/dakota.err");
	
	// (1) read JSON file
	jsonInput inp(workDir);

	//
	//	Construct Nataf Object
	//

	ERANataf T(inp);

	//
	//	Random number generator N(0,1)
	//

	int nmc = inp.nmc;
	std::mt19937 generator(inp.rseed);
	std::normal_distribution<double> distribution(0.0, 1.0);

	vector<vector<double>> uvals(nmc, vector<double>(inp.nrv, 0.0));
	for (int ns = 0; ns < nmc; ns++)
	{
		for (int nr = 0; nr < inp.nrv; nr++)
			uvals[ns][nr] = distribution(generator);
	}
	
	//
	//	FE Analysis
	//

	vector<vector<double>> gvals, xvals;
	T.simulateApp(osType, runType, inp, uvals, xvals, gvals);

	//
	//	Global sensitivity analysis
	//

	int Kos = 25;
	runGSA GsaResults(xvals, gvals, inp.groups, Kos); 
	
	vector<vector<double>> Si = GsaResults.Simat;
	vector<vector<double>> St = GsaResults.Stmat;

	//
	//	Write files/dakota.out, dakotaTab.out
	//

	writeOutputs(inp, xvals, gvals, GsaResults);
	
	return 0;

}

void writeOutputs(jsonInput inp, 
					vector<vector<double>> xvals, 
					vector<vector<double>> gvals, 
					runGSA GsaResults)
{
	// dakota.out
	string writingloc = inp.workDir + "/dakota.out";
	std::ofstream outfile(writingloc);

	if (!outfile.is_open()) {
		theErrorFile << "Error running UQ engine: Unable to write dakota.out";
		theErrorFile.close();
		exit(-1);
	}

	outfile.setf(std::ios::fixed, std::ios::floatfield); // set fixed floating format
	outfile.precision(4); // for fixed format

	outfile << "* number of input combinations" << std::endl;
	outfile << inp.ngr << std::endl;

	outfile << "* input names" << std::endl;
	for (int i = 0; i < inp.ngr; i++) {
		for (int j = 0; j < inp.groups[i].size() - 1; j++) {
			outfile << inp.rvNames[inp.groups[i][j]] << ",";
		}
		outfile << inp.rvNames[inp.groups[i][inp.groups[i].size() - 1]] << std::endl;
	}

	outfile << "* number of outputs" << std::endl;
	outfile << inp.nqoi << std::endl;

	outfile << "* output names" << std::endl;
	for (int i = 0; i < inp.nqoi; i++) {
		outfile << inp.qoiNames[i] << std::endl;
	}

	outfile << "* ";
	for (int j = 0; j < inp.ngr; j++) {
		outfile << "Sm(" << std::to_string(j + 1) << ")  ";
	}
	for (int j = 0; j < inp.ngr; j++) {
		outfile << "St(" << std::to_string(j + 1) << ")  ";
	}
	outfile << std::endl;

	for (int i = 0; i < inp.nqoi; i++) {

		for (int j = 0; j < inp.ngr; j++) {
			outfile << GsaResults.Simat[i][j] << " ";
		}
		for (int j = 0; j < inp.ngr; j++) {
			outfile << GsaResults.Stmat[i][j] << " ";
		}
		outfile << std::endl;
	}

	outfile << "* number of samples" << std::endl;
	outfile << inp.nmc << std::endl;
	outfile.close();

	// dakotaTab.out
	std::string writingloc1 = inp.workDir + "/dakotaTab.out";
	std::ofstream Taboutfile(writingloc1);

	if (!Taboutfile.is_open()) {
		theErrorFile << "Error running UQ engine: Unable to write dakotaTab.out";
		theErrorFile.close();
		exit(-1);
	}

	Taboutfile.setf(std::ios::fixed, std::ios::floatfield); // set fixed floating format
	Taboutfile.precision(10); // for fixed format

	Taboutfile << "idx         ";
	for (int j = 0; j < inp.nrv + inp.nco; j++) {
		Taboutfile << inp.rvNames[j] << "           ";
	}
	for (int j = 0; j < inp.nqoi; j++) {
		Taboutfile << inp.qoiNames[j] << "            ";
	}
	Taboutfile << std::endl;


	for (int i = 0; i < inp.nmc; i++) {
		Taboutfile << std::to_string(i + 1) << "    ";
		for (int j = 0; j < inp.nrv + inp.nco; j++) {
			Taboutfile << std::to_string(xvals[i][j]) << "    ";
		}
		for (int j = 0; j < inp.nqoi; j++) {
			Taboutfile << std::to_string(gvals[i][j]) << "    ";
		}
		Taboutfile << std::endl;
	}

	theErrorFile.close();

}