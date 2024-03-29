
#include "ERANataf.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <filesystem>
#if defined(_OPENMP)
  #include <omp.h>
#endif
//

//#include <chrono>
#include "nlopt.hpp"


double natafObjec(unsigned n, const double* rho0, double* grad, void* my_func_data);
using boost::math::normal;

ERANataf::ERANataf() {}

ERANataf::ERANataf(jsonInput inp)
{

	//
	//	Define Marginal Distributions & corr
	//
	this->nrv = inp.nrv;

	vector<ERADist> M(nrv);
	for (int i = 0; i < nrv; i++)
	{
		M[i] = *new ERADist(inp.distNames[i], inp.opts[i], inp.vals[i], inp.adds[i]);
	}
	vector<vector<double>> Rhox = inp.corr;

	std::cout << "[RhoxExt]" << std::endl;
	for (int nri = 0; nri < nrv; nri++)
	{
		for (int nrj = 0; nrj < nrv; nrj++)
			printf("%3.2f  ", Rhox[nri][nrj]);
		printf("\n");
	}
	printf("\n");

	this->M = M;


	//M_ = M;
	// initialize Rhoz = eye(nrv)
	for (int nr = 0; nr < nrv; nr++)
	{
		vector<double> Rhoz_row(nrv, 0.0);
		Rhoz_row[nr] = 1.0;
		Rhoz.push_back(Rhoz_row);
	}


	// Grid points
	int ngrid = 256;
	double zmax = 8.;
	double zmin = -zmax;

	vector<double> points;
	vector<double> weights;
	quadGL(ngrid, zmin, zmax, *&points, *&weights);
	vector<vector<double>> fxiTmp(nrv, vector<double>(ngrid, 0.0));

	// Nataf
	for (int nri = 0; nri < nrv; nri++)
	{
		auto mydisti = M[nri].theDist;
		string namei = mydisti->getName();
		double meani = mydisti->getMean();
		double stdi  = mydisti->getStd();

		for (int nrj = 0; nrj < nri; nrj++)
		{
			auto mydistj = M[nrj].theDist;
			string namej = mydistj->getName();
			double meanj = mydistj->getMean();
			double stdj  = mydistj->getStd();

			double Rhoxij = Rhox[nri][nrj];
			if (Rhoxij == 0)
			{
				continue;
			}
			else if ((namei.compare("normal") == 0) && (namej.compare("normal") == 0))
			{
				Rhoz[nri][nrj] = Rhoxij;
				Rhoz[nrj][nri] = Rhoxij;
				continue;
			}
			else if ((namei.compare("normal") == 0) && (namej.compare("lognormal") == 0))
			{
				double Vj = stdj / meanj;
				Rhoz[nri][nrj] = Rhoxij * Vj / sqrt(log(1.0 + Vj * Vj));
				Rhoz[nrj][nri] = Rhoz[nri][nrj];
				continue;
			}
			else if ((namei.compare("lognormal") == 0) && (namej.compare("normal") == 0))
			{
				double Vi = stdi / meani;
				Rhoz[nri][nrj] = Rhoxij * Vi / sqrt(log(1.0 + Vi * Vi));
				Rhoz[nrj][nri] = Rhoz[nri][nrj];
				continue;
			}
			else if ((namei.compare("lognormal") == 0) && (namej.compare("lognormal") == 0))
			{
				double Vj = stdj / meanj;
				double Vi = stdi / meani;
				Rhoz[nri][nrj] = (log(1.0 + Rhoxij * Vi * Vj) / sqrt(log(1.0 + Vi * Vi) * log(1 + Vj * Vj)));
				Rhoz[nrj][nri] = Rhoz[nri][nrj];
				continue;
			}
			else
			{
				// solving Nataf
				double Rhozij[1] = { Rhoxij };
				
				if (fxiTmp[nri][0] ==0)
					for (int ng = 0; ng < ngrid; ng++)
						fxiTmp[nri][ng] = (mydisti->getQuantile(cdf(stdNorm, points[ng])) - meani) / (stdi) *weights[ng];
				
				if (fxiTmp[nrj][0] == 0)
					for (int ng = 0; ng < ngrid; ng++)
						fxiTmp[nrj][ng] = (mydistj->getQuantile(cdf(stdNorm, points[ng])) - meanj) / (stdj) * weights[ng];

				my_NatafInfo addVars;
				addVars.points = points;
				addVars.fxii = fxiTmp[nri];
				addVars.fxij = fxiTmp[nrj];

				addVars.Rhoxij = Rhoxij;
				addVars.ngrid = ngrid;

				// optimization
				double lb[1] = { -1. };
				double ub[1] = { 1. };
				nlopt_opt optim = nlopt_create(NLOPT_LN_BOBYQA, 1); // derivative-free algorithm
				nlopt_set_upper_bounds(optim, ub);

				nlopt_set_min_objective(optim, natafObjec, &addVars);
				nlopt_set_xtol_rel(optim, 1e-6);
				double minf; // `*`the` `minimum` `objective` `value,` `upon` `return`*` 

				//double params;
				
				if (nlopt_optimize(optim, Rhozij, &minf) < 0) {
					std::string errMsg = "Error running UQ engine: Nataf optimization failed (nlopt failed)";
					std::cout << errMsg << "\n";
					theErrorFile << errMsg << std::endl;
					theErrorFile.close();
					exit(-1);
				}
				else if (minf > 1.e-5)
				{
					std::string errMsg = "Error running UQ engine: Nataf optimization did not converge (nlopt failed)";
					std::cout << errMsg << "\n";
					theErrorFile << errMsg << std::endl;
					theErrorFile.close();
					exit(-1);
				}
				else {
					//printf("found minimum at f(%g) = %0.10g\n", Rhozij[0], minf);
				}

				Rhoz[nri][nrj] = Rhozij[0];
				Rhoz[nrj][nri] = Rhozij[0];
			}
		}
	}

	RhozMat = Eigen::MatrixXd(nrv, nrv);
	for (int i = 0; i < nrv; ++i) {
		for (int j = 0; j < nrv; ++j) {
			RhozMat(i, j) = Rhoz[i][j];
		}
	}

	if (nrv>0) {
	    RhozInv = RhozMat.inverse();
	}
	// Cholesky decomposition of correlation matrix
	Eigen::LLT<Eigen::MatrixXd> llt(RhozMat);
	
	if (llt.info() == Eigen::NumericalIssue)
	{
		std::string errMsg = "Error running UQ engine: Nataf transformation is not applicable (not positive definite)";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}
	else {
		A = llt.matrixL(); // 'Lower'
	}

	std::cout << "[RhoOpt]" << std::endl;
	for (int nri = 0; nri < nrv; nri++)
	{
		for (int nrj = 0; nrj < nrv; nrj++)
			printf("%3.2f  ", RhozMat(nri,nrj));
		printf("\n");
	}
	printf("\n");
}
ERANataf::~ERANataf() {}

vector<vector<double>> ERANataf::X2U(int nmc, vector<vector<double>> x)
{
	vector<vector<double>> u(nmc, vector<double>(nrv, 0.0));
	Eigen::MatrixXd zvec(nrv, nmc);
	
	//normal stdNorm(0., 1.);
	for (int nr = 0; nr < nrv; nr++)
	{
		RVDist* theDist = M[nr].theDist;
		for (int ns = 0; ns < nmc; ns++)
			zvec(nr,ns) = quantile(stdNorm, theDist->getCdf(x[ns][nr]));
	}
	Eigen::VectorXd uvec = A.colPivHouseholderQr().solve(zvec);

	for (int nr = 0; nr < nrv; nr++)
	{
		for (int ns = 0; ns < nmc; ns++)
			u[ns][nr]=uvec(ns,nr);
	}
	return u;
}

vector<vector<double>> ERANataf::U2X(int nmc, vector<vector<double>> u)
{
	vector<vector<double>> x(nmc, vector<double>(nrv, 0.0));
	Eigen::MatrixXd uvec(nrv,nmc); // transpose
	for (int ns = 0 ; ns <nmc ; ns++)
	{ 
		for (int nr = 0; nr < nrv; nr++)
			uvec(nr,ns) = u[ns][nr];
	}

	Eigen::MatrixXd zvec = A*uvec;

	for (int nr = 0; nr < nrv; nr++)
	{

		//auto start_ = std::chrono::high_resolution_clock::now();

		RVDist *theDist = M[nr].theDist;
		for (int ns = 0; ns < nmc; ns++)
		{ 
			x[ns][nr] = theDist->getQuantile(normCdf(zvec(nr, ns)));
		}
		
		//auto finish_ = std::chrono::high_resolution_clock::now();
		//std::chrono::duration<double> elapsed_ = finish_ - start_;
		//std::cout << "++++++++++++++++++ Estim - Elapsed time: " << elapsed_.count() << " s\n";
	}

	return x;
}


double ERANataf::getJointPdf(vector<double> x)
{
	Eigen::VectorXd u(nrv);
	Eigen::VectorXd phi(nrv);
	Eigen::VectorXd f(nrv);
	for (int nr = 0; nr < nrv; nr++)
	{
		u(nr)=quantile(stdNorm, M[nr].theDist->getCdf(x[nr]));
		phi(nr)=pdf(stdNorm, u[nr]);
		f(nr)=M[nr].theDist->getPdf(x[nr]);
	}
	double phi_n = mvnpdfR(u);
	double jointPdf = 1;
	for (int nr = 0; nr < nrv; nr++)
	{
		jointPdf *= f(nr) / phi(nr) ;
	}
	jointPdf *= phi_n;
	return jointPdf;
}



double ERANataf::getJointCdf(vector<double> x)
{
	Eigen::VectorXd u(nrv);
	for (int nr = 0; nr < nrv; nr++)
	{
		u(nr) = quantile(stdNorm, M[nr].theDist->getCdf(x[nr]));
	}
	std::cout << "++++++++++++++++" << std::endl;
	std::cout << RhozMat << std::endl;
	std::cout << "++++++++++++++++" << std::endl;
	return mvncdfR(u);
}




void ERANataf::quadGL(int N, double a, double b, vector<double> &x, vector<double> &w)
{
	// Modified from ERANataf.m
	//(Written by Greg von Winckel - 02 / 25 / 2004)

	N = N - 1;
	int N1 = N + 1;
	int N2 = N + 2;

	vector<double> y, y0, Lp;
	double maxEps = 0;
	for (int i = 0; i < N+1; i++)
	{
		// xu = linspace(-1, 1, N1)';
		double xval = 2. / (N) *i - 1.;
		double yi = cos((2. * i + 1.) * PI / (2. * N + 2.)) + (0.27 / N1) * sin(PI * xval * N / N2);
		y.push_back(yi);
		y0.push_back(2.);
		Lp.push_back(0.);
		maxEps = std::max(abs(y[i] - y0[i]), maxEps);
	}


	while (maxEps > 1.e-10)
	{
		//vector<double> L_ = 1;
		//vector<double> L = y;

		vector<double> Lold,L;
		for (int i = 0; i < N1; i++)
		{
			Lold.push_back(1.);
			L.push_back(y[i]);
		}

		for (int k = 2; k < N1 + 1; k++)
		{
			for (int i = 0; i < N1; i++)
			{
				double Lnew = ((2. * k - 1.) * y[i] * L[i] - (k - 1.) * Lold[i]) / k;
				Lold[i] = L[i];
				L[i] = Lnew;
			}
		}
		maxEps = 0;
		for (int i = 0; i < N1; i++)
		{
			Lp[i] = N2*(Lold[i]-y[i]*L[i])/(1-y[i]*y[i]);
			y0[i] = y[i];
			y[i] = y0[i] - L[i] / Lp[i];
			maxEps = std::max(abs(y[i] - y0[i]), maxEps);
		}
	}

	for (int i = 0; i < N1; i++)
	{
		x.push_back((a * (1 - y[i]) + b * (1 + y[i])) / 2);
		w.push_back((b - a) / ((1 - y[i] * y[i]) * Lp[i] * Lp[i]) * (N2 / N1) * (N2 / N1));
	}

}

double ERANataf::mvnpdfR(Eigen::VectorXd u)
{

	Eigen::VectorXd& x = u;
	double sqrt2pi = std::sqrt(2 * PI);
	double quadform = (u).transpose() * RhozInv * (u);
	double norm = std::pow(sqrt2pi, -nrv)*std::pow(RhozMat.determinant(), -0.5);

	return norm * exp(-0.5 * quadform);
}


double ERANataf::mvncdfR(Eigen::VectorXd u)
{

	if (nrv == 1)
	{
		return cdf(stdNorm, u(0));
	}
	else
	{ 
		// Genz, A. and F. Bretz (1999) "Numerical Computation of Multivariat t Probabilities with Application to Power Calculation of Multiple
		// Contrasts", J.Statist.Comput.Simul., 63:361-378. - 

		// A is C in genz
		// a=-inf in genz
		// b=u    in genz

		std::default_random_engine generator;
		std::uniform_real_distribution<double> distribution(0.0, 1.0);

		double errEst = 0;
		int N = 0;
		double intVal = 0.0, varSum = 0.0;
		double gam = 0.997;
		int Nmax = 200000;
		while ((errEst > 1.e-5 || N<5) && (N<Nmax))
		{
			// Initialize
			double d = 0;
			double e = cdf(stdNorm, u(0) / A(0, 0));
			double f = e - d;
			vector<double> y;
			double sumcy;

			// Loop to evaluate fq
			for (int nri = 0 ; nri < nrv-1; nri++)
			{
				double w = distribution(generator);
				y.push_back(quantile(stdNorm, (d + w * (e - d))));
				sumcy = 0.0;
				for (int nrj = 0; nrj < nri+1; nrj++)
					sumcy += A(nri + 1, nrj) * y[nrj];
				d = 0;
				e = cdf(stdNorm, (u(nri + 1) - sumcy) / A(nri + 1, nri + 1));
				f = (e - d) * f;
			}
			N++;
			varSum += (N-1)*(f-intVal) * (f - intVal)/N;
			intVal += (f-intVal)/N;
			errEst = gam * std::sqrt(varSum/N / (N - 1));
			//intsum += f(m - 1);
			//varsum += f(m - 1) * f(m - 1);
		}
		return intVal;
	}

	
}

double ERANataf::normCdf(double x)
{
	// from http://www.johndcook.com/cpp_phi.html
	// constants
	double a1 = 0.254829592;
	double a2 = -0.284496736;
	double a3 = 1.421413741;
	double a4 = -1.453152027;
	double a5 = 1.061405429;
	double p = 0.3275911;

	// Save the sign of x
	int sign = 1;
	if (x < 0)
		sign = -1;
	x = fabs(x) / sqrt(2.0);

	// A&S formula 7.1.26
	double t = 1.0 / (1.0 + p * x);
	double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x);

	return 0.5 * (1.0 + sign * y);
}

void ERANataf::simulateAppBatch(string osType, string runType, jsonInput inp, vector<vector<double>> u, vector<vector<int>> resampIDs, vector<vector<double>> &x, vector<vector<double>> &gval)
{
	//
	// If we find result.out in the templete dir. emit error;
	//

	std::string existingResultsFile = inp.workDir + "/templatedir/results.out";
	if (std::filesystem::exists(existingResultsFile)) {
		//*ERROR*
		std::string errMsg = "Error running SimCenterUQ: your templete directory already contains results.out file. Please clean up the directory where input file is located.";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	std::string existingParamsFile = inp.workDir + "/templatedir/results.out";
	if (std::filesystem::exists(existingParamsFile)) {
		//*ERROR*
		std::string errMsg = "Error running SimCenterUQ: your templete directory already contains params.in file. Please clean up the directory where input file is located.";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	//
	// Change from u to x;
	//

	std::cout << "[First 10 samples]" << std::endl;

	x = U2X(inp.nmc, u);
	std::vector<double> zero_vector(inp.nre, 0);
	for (int ns = 0; ns < inp.nmc; ns++)
	{
		// for resampling
		x[ns].insert(x[ns].end(), zero_vector.begin(), zero_vector.end());
		for (int ng = 0; ng < inp.nreg; ng++)
		{
			for (int nr : inp.resamplingGroups[ng])
			{
				x[ns][nr] = inp.vals[nr][resampIDs[ns][ng]];
			}
		}
		// for constants
		for (int j = 0; j < inp.nco; j++)
			x[ns].push_back(inp.constants[j]);

	}
	for (int ns = 0; ns < std::min(inp.nmc,10) ; ns++)
	{
		for (int nr = 0; nr < inp.nrv; nr++)
		{
			std::cout << x[ns][nr] << "  ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	//
	// Run Apps
	//

	std::string copyDir = inp.workDir + "/templatedir";

	std::cerr << "workdir:" << inp.workDir << "\n";
	std::cerr << "copyDir:" << copyDir << "\n";
	std::cerr << "runningFEM analysis.." << "\n\n";

	const size_t COL = 1;
	const size_t  ROW = 1;
	//gval.resize(inp.nmc, vector<int>(inp.nqoi, 0));

	gval.resize(inp.nmc);
	int i;
        #pragma omp parallel for shared(gval) private(i)
	for (i = 0; i < inp.nmc; i++)
	{
		gval[i] = simulateAppOnce(i, inp.workDir, copyDir, inp.nrv + inp.nco + inp.nre, inp.nqoi, inp.rvNames, x[i], osType, runType);
	}

}

vector<double> ERANataf::simulateAppOnce(int i, string workingDirs, string copyDir, int nrvcore, int nqoi, vector<string> rvNames, vector<double> xs, string osType, string runType)
{

	//
	// (1) create "workdir.i " folder :need C++17 to use the files system namespace 
	//

	string workDir = workingDirs + "/workdir." + std::to_string(i + 1);

	//std::cerr << "workDir:" + workDir + "\n";

	//
	// (2) copy files from templatedir to workdir.i
	//

	const auto copyOptions =
		std::filesystem::copy_options::update_existing
		| std::filesystem::copy_options::recursive;

	/*
  const auto copyOptions = std::filesystem::copy_options::overwrite_existing;
	*/

	try
	{
		std::filesystem::copy(copyDir, workDir, copyOptions);
	}
	catch (std::exception & e)
	{
		std::cout << e.what() << "\n";
	}


	//std::filesystem::current_path(workDir); //======= Not good for parallel

	/*
	if (ok != true) {
	  std::cerr << "my_nataf - could not copy files to " << workDir << "\n";
	}
	*/

	//
	// (3) write param.in file
	//

	string params = workDir + "/params.in";
	std::ofstream writeFile(params.data());
	if (writeFile.is_open()) {
		writeFile << std::to_string(nrvcore) + "\n";
		for (int j = 0; j < nrvcore; j++) {
			writeFile << rvNames[j] + " ";
			writeFile << std::to_string(xs[j]) + "\n";
		}
		writeFile.close();
	}

	//
	// (4) run workflow_driver.bat(e.g. It will make "SimCenterInput.tcl" and run OpenSees)
	//

	std::string workflowDriver = "workflow_driver";
	if ((osType.compare("Windows") == 0) && (runType.compare("runningLocal") == 0)) {
	    workflowDriver = "workflow_driver.bat";
        //workflowDriver = "workflow_driver.bat >nul 2>nul";
	}
	string workflowDriver_string = "cd " + workDir + " && " + workDir + "/" + workflowDriver;

	const char* workflowDriver_char = workflowDriver_string.c_str();
	system(workflowDriver_char);
	if (i == 1) {
		std::cout << workflowDriver_char << "\n\n";
	}
	//
	// (5) get the values in "results.out"
	//

	string results = workDir + "/results.out";
	std::ifstream readFile(results.data());

	if (!readFile.is_open()) {
		//*ERROR*
		std::string errMsg = "Error running FEM: results.out missing in workdir." + std::to_string(i + 1) + ". Check your FEM inputs.";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	vector<double> g_tmp;
	if (readFile.is_open()) {
		int j = 0;
		double g;
		while (readFile >> g) {
			g_tmp.push_back(g);
			j++;
		}
		readFile.close();

		if (j == 0) {
			std::string errMsg = "Error running FEM: results.out file at workdir." + std::to_string(i + 1) + " is empty.";
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		if (j != nqoi) {
			//*ERROR*
			std::string errMsg = "Error reading FEM results: the number of outputs in results.out (" + std::to_string(j) + ") does not match the number of QoIs specified (" + std::to_string(nqoi) + ")";
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);
		}
	}

	return g_tmp;
}


void ERANataf::simulateAppSequential(string osType, string runType, jsonInput inp, vector<vector<double>> u, vector<vector<double>>& xval, vector<vector<double>>& gval, int idx)
{
	//
	// Change from u to x;
	//

	vector<vector<double>> x = U2X(1, u);
	for (int j = 0; j < inp.nco; j++)
		x[0].push_back(inp.constants[j]);

	//
	// Run Apps
	//

	std::string copyDir = inp.workDir + "/templatedir";

	//std::cerr << "workdir:" + inp.workDir + "\n";
	std::cerr << "copyDir:" + copyDir + "\n";
	
	//
	// (1) create "workdir.idx " folder :need C++17 to use the files system namespace 
	//

	string workDir = inp.workDir + "/workdir." + std::to_string(idx + 1);

	std::cerr << "workDir:" + workDir + "\n";

	//
	// (2) copy files from templatedir to workdir.idx
	//

	const auto copyOptions =
		std::filesystem::copy_options::update_existing
		| std::filesystem::copy_options::recursive;

	/*
	const auto copyOptions = std::filesystem::copy_options::overwrite_existing;
	*/

	try
	{
		std::filesystem::copy(copyDir, workDir, copyOptions);
	}
	catch (std::exception & e)
	{
		std::cout << e.what();
	}


	std::filesystem::current_path(workDir);

	/*
	if (ok != true) {
		std::cerr << "my_nataf - could not copy files to " << workDir << "\n";
	}
	*/

	//
	// (3) write param.in file
	//

	string params = workDir + "/params.in";
	std::ofstream writeFile(params.data());
	if (writeFile.is_open()) {
		writeFile << std::to_string(inp.nrv + inp.nco) + "\n";
		for (int j = 0; j < inp.nrv+inp.nco; j++) {
			writeFile << inp.rvNames[j] + " ";
			writeFile << std::to_string(x[0][j]) + "\n";
		}
		writeFile.close();
	}
	
	//
	// (4) run workflow_driver.bat(e.g. It will make "SimCenterInput.tcl" and run OpenSees)
	//

	std::string workflowDriver = "workflow_driver";
	if ((osType.compare("Windows") == 0) && (runType.compare("runningLocal") == 0))
		workflowDriver = "workflow_driver.bat";

	string workflowDriver_string = workDir + "/" + workflowDriver;

	const char* workflowDriver_char = workflowDriver_string.c_str();
	system(workflowDriver_char);

	//
	// (5) get the values in "results.out"
	//

	string results = workDir + "/results.out";
	std::ifstream readFile(results.data());

	if (!readFile.is_open()) {
		//*ERROR*

		std::string errMsg = "Error reading FEM results: check your inputs ";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}



	vector<double> g_tmp(inp.nqoi);
	if (readFile.is_open()) {
		int j = 0;
		double g;
		while (readFile >> g) {
			g_tmp[j] = g;
			j++;
		}
		readFile.close();
	}

	gval.push_back(g_tmp);
	xval.push_back(x[0]);
}




double natafObjec(unsigned n, const double* rho0, double* grad, void* my_func_data)
{

	my_NatafInfo* info = (my_NatafInfo*)my_func_data;
	vector<double> points = info->points;
	vector<double> fxii = info->fxii;
	vector<double> fxij = info->fxij;
	double Rhoxij = info->Rhoxij;
	double ngrid = info->ngrid;
	const double PI = 3.1415926535897932384626433;

	//normal stdNorm(0., 1.);
	double integVal = 0;
	//double rho_2xi, xi, xi2, eta, eta2, feta, fxi, bNormPdf;
	double xi, eta, feta;

	int nn = 0;
	double rho = rho0[0];
	double coef1 = 1. / (2. * (1. - rho0[0] * rho0[0]));
	double coef2 = 1. / (2. * PI * sqrt(1. - rho0[0] * rho0[0]));
	for (int ngi = 0; ngi < ngrid; ngi++)
	{
		xi = points[ngi];
		feta = fxii[ngi];
		for (int ngj = 0; ngj < ngrid; ngj++)
		{
			eta = points[ngj];
			//fxi = fxij[ngj];
			//bNormPdf = (1. / (2. * PI * sqrt(1. - rho0[0] * rho0[0])) * exp(-1. / (2. * (1. - rho0[0] * rho0[0])) * (xi * xi - 2. * rho0[0] * xi * eta + eta * eta)));
			integVal += feta * fxij[ngj] * exp(- coef1 * ((xi - eta)*(xi - eta) + 2*(1-rho)*xi*eta));
		}
	}
	integVal = coef2 * integVal;
	//printf("*");
	return abs(integVal - Rhoxij);
}


