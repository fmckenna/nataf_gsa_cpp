#include "runForward.h"

#include <iterator>

runForward::runForward() {}

runForward::runForward(vector<vector<double>> xval,	vector<vector<double>> gmat)
{
	std::cout << std::string("============================FLAG+=============================") << std::endl;

	this->xval = xval;
	this->gval = gmat;
	nmc = xval.size();
	nrv = xval[0].size();

	for (int nr = 0; nr < nrv; nr++) {
		vector<double> xvec;
		for (int ns = 0; ns < nmc; ns++) {
			xvec.push_back(xval[ns][nr]);
		}

		double mean_val = calMean(xvec);
		double stdDev_val = calStd(xvec, mean_val);
		double skewness_val = calSkewness(xvec, mean_val, stdDev_val);
		double kurtosis_val = calKurtosis(xvec, mean_val, stdDev_val);

		std::cout << mean_val << std::endl;
		std::cout << stdDev_val << std::endl;
		std::cout << skewness_val << std::endl;
		std::cout << kurtosis_val << std::endl;

		mean.push_back(mean_val);
		stdDev.push_back(stdDev_val);
		skewness.push_back(skewness_val);
		kurtosis.push_back(kurtosis_val);
	}	
}

runForward::~runForward() {};

double runForward::calMean(vector<double> x) {
	double sum = std::accumulate(std::begin(x), std::end(x), 0.0);
	return sum / x.size();
}

double runForward::calStd(vector<double> x, double m) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m);
		});
	return std::sqrt(accum / (x.size()));
}

double runForward::calSkewness(vector<double> x, double m, double s) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m) * (d - m);
		});
	return (accum / (x.size()) / (s*s*s));
}

double runForward::calKurtosis(vector<double> x, double m, double s) {
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m) * (d - m) * (d - m);
		});
	return (accum / (x.size()) / (s * s * s * s));
	 
}
double runForward::writeOutputs(jsonInput inp)
{
	// dakota.out
	string writingloc = inp.workDir + "/dakota.out";
	std::ofstream outfile(writingloc);


	if (!outfile.is_open()) {
		theErrorFile << "Error running UQ engine: Unable to write dakota.out";
		theErrorFile.close();
		exit(-1);
	}


	json outJson;

	outJson["rvNames"] = inp.rvNames;
	outJson["mean"] = mean;
	outJson["standardDeviation"] = stdDev;
	outJson["skewness"] = skewness;
	outJson["kurtosis"] = kurtosis;

	//outfile.setf(std::ios::fixed, std::ios::floatfield); // set fixed floating format
	//outfile.precision(4); // for fixed format

	outfile << outJson.dump(4) << std::endl;
	//outfile << outJson << std::endl;

	//outfile << "* number of input combinations" << std::endl;
	//outfile << inp.ngr << std::endl;

	//outfile << "* input names" << std::endl;
	//for (int i = 0; i < inp.ngr; i++) {
	//	for (int j = 0; j < inp.groups[i].size() - 1; j++) {
	//		outfile << inp.rvNames[inp.groups[i][j]] << ",";
	//	}
	//	outfile << inp.rvNames[inp.groups[i][inp.groups[i].size() - 1]] << std::endl;
	//}

	//outfile << "* number of outputs" << std::endl;
	//outfile << inp.nqoi << std::endl;

	//outfile << "* output names" << std::endl;
	//for (int i = 0; i < inp.nqoi; i++) {
	//	outfile << inp.qoiNames[i] << std::endl;
	//}

	//outfile << "* ";
	//for (int j = 0; j < inp.ngr; j++) {
	//	outfile << "Sm(" << std::to_string(j + 1) << ")  ";
	//}
	//for (int j = 0; j < inp.ngr; j++) {
	//	outfile << "St(" << std::to_string(j + 1) << ")  ";
	//}
	//outfile << std::endl;

	//for (int i = 0; i < inp.nqoi; i++) {

	//	for (int j = 0; j < inp.ngr; j++) {
	//		outfile << Simat[i][j] << " ";
	//	}
	//	for (int j = 0; j < inp.ngr; j++) {
	//		outfile << Stmat[i][j] << " ";
	//	}
	//	outfile << std::endl;
	//}

	//outfile << "* number of samples" << std::endl;
	//outfile << inp.nmc << std::endl;
	//outfile.close();

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
	for (int j = 0; j < inp.nrv + inp.nco + inp.nre; j++) {
		Taboutfile << inp.rvNames[j] << "           ";
	}
	for (int j = 0; j < inp.nqoi; j++) {
		Taboutfile << inp.qoiNames[j] << "            ";
	}
	Taboutfile << std::endl;


	for (int i = 0; i < inp.nmc; i++) {
		Taboutfile << std::to_string(i + 1) << "    ";
		for (int j = 0; j < inp.nrv + inp.nco + inp.nre; j++) {
			Taboutfile << std::to_string(xval[i][j]) << "    ";
		}
		for (int j = 0; j < inp.nqoi; j++) {
			Taboutfile << std::to_string(gval[i][j]) << "    ";
		}
		Taboutfile << std::endl;
	}

	theErrorFile.close();

}