#include "runGSA.h"

#include <iterator>

runGSA::runGSA() {}

runGSA::runGSA(vector<vector<double>> xval,
	vector<vector<double>> gmat,
	vector<vector<int>> combs_tmp,
	int Kos)
{
	this->xval = xval;
	this->gval = gmat;
	this->combs_tmp = combs_tmp;
	nmc = xval.size();

	nrv = xval[0].size();
	ncombs = combs_tmp.size(); 
	int nqoi = gmat[0].size();
	Kos = std::min(Kos, int(nmc / 5));

	// for each QoI
	for (int j = 0; j < nqoi; j++) {

		vector<double> gvec;
		double sqDiff = 0;
		gvec.reserve(nmc);
		for (int i = 0; i < nmc; i++) {
			gvec.push_back(gmat[i][j]);
			
		}
		

		// check if the variance is zero
		double mean = 0;
		for (int i = 0; i < nmc; i++)
			mean += gvec[i];

		mean = mean / double(nmc);
		for (int i = 0; i < nmc; i++)
			sqDiff += (gmat[i][j] - mean) * (gmat[i][j] - mean);

		//double var = sqDiff / nmc;
		if (sqDiff < 1.e-10) {
			theErrorFile << "Error running FEM: the variance of output is zero. Output value is " << mean;
			theErrorFile.close();
			exit(1);
		};

		vector<double> Sij, Stj;


		// repeat GSA with different Kos untill success
		double failIdx = -100, i = 1;
		while ((failIdx == -100) || (Kos/i < 0.5))
		{
			Sij = doGSA(gvec, ceil(Kos/i),'M');
			failIdx = Sij[0];
			i *=2;
		}

		failIdx = -100, i = 1;
		while ((failIdx == -100) || (Kos / i < 0.5))
		{
			Stj = doGSA(gvec, ceil(Kos/i), 'T');
			failIdx = Stj[0];
			i *= 2;
		}

		vector<double> Si_temp, Kos,St_temp;


		Simat.push_back(Stj);
		Stmat.push_back(Sij);
	}
}

vector<double> runGSA::doGSA(vector<double> gval,int Kos,char Opt)
{
	vector<vector<int>> combs;

	if (Opt == 'T')
	{
		vector<int> allSet(nrv);
		std::iota(allSet.begin(), allSet.end(), 0);

		for (auto comb : combs_tmp)
		{
			vector<int> cnc;
			//std::set_difference(allSet.begin(), allSet.end(), comb.begin(), comb.end(), cnc.begin());
			std::set_difference(allSet.begin(), allSet.end(), comb.begin(), comb.end(), std::inserter(cnc, cnc.begin()));
			combs.push_back(cnc);
		}
	}
	else
	{
		combs = combs_tmp;
	}

	double V = calVar(gval);
	double Vi;
	vector<double> Si;
	Si.reserve(ncombs);

	for (int nc = 0; nc < ncombs; nc++)
	{

		const int endm = combs[nc].size(); // (nx+ng)-1
		const int endx = endm - 1;			// (nx)-1

		if (endm == 0)
		{
			if (Opt == 'T')
			{
				Si.push_back(1.); // total
			}
			else
			{
				Si.push_back(0.);   // main
			}
			printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);
			continue;
		}
		else if (endm == nrv)
		{
			if (Opt == 'T')
			{
				Si.push_back(0.); // total
			}
			else
			{
				Si.push_back(1.);   // main
			}
			printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);
			continue;
		}

		mat data(endm + 1, nmc);

		for (int ne = 0; ne < endm; ne++)
		{
			int idx = combs[nc][ne];

			if (idx > nrv - 1) {
				theErrorFile << "Error running UQ engine: combination index exceeds the bound" << std::endl;
				theErrorFile.close();
				exit(-1);
			}

			for (int ns = 0; ns < nmc; ns++)
			{
				data(ne, ns) = xval[ns][idx];
			}
		}


		for (int ns = 0; ns < nmc; ns++)
		{
			data(endm, ns) = gval[ns];
		}

		gmm_full model;
		bool status = model.learn(data, Kos, maha_dist, random_subset, 30, 100, V *1.e-3, false);

		if (status == false)
		{
			theErrorFile << "Error running UQ engine: GSA learning failed" << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		mat mu = model.means;   //nrv x Ko
		cube cov = model.fcovs; //nrv x nrv x Ko
		rowvec pi = model.hefts;   //1 x Ko 
		rowvec mug = mu.row(endm);    //1 x Ko

		vector<double> mui;
		mui.reserve(nmc);
		// Used to calculate conditional mean and covariance
		cube SiginvSig(1, endm, Kos);
		mat muk(endm, Kos);
		for (int k = 0; k < Kos; k++)
		{
			mat Sig12 = cov.subcube(0, endm, k, endx, endm, k);
			mat Sig11 = cov.subcube(0, 0, k, endx, endx, k);
			muk.col(k) = mu.submat(0, k, endx, k);
			SiginvSig.slice(k) = solve(Sig11, Sig12).t();
		}

		//model.means.print("means:");
		//model.fcovs.print("fcovs:");

		for (int i = 0; i < nmc; i++)
		{
			rowvec pik_tmp(Kos, fill::zeros);
			colvec muki(Kos);
			mat xi = data.submat(0, i, endx, i);

			for (int k = 0; k < Kos; k++)
			{
				mat tmp = SiginvSig.slice(k);
				mat muval = muk.col(k);
				muki.subvec(k, k) = mug(k) + SiginvSig.slice(k) * (xi - muval);
				pik_tmp(k) = pi(k) * mvnPdf(xi, muval, cov.subcube(0, 0, k, endx, endx, k));

			}

			rowvec piki = pik_tmp / sum(pik_tmp);
			mat tmp = piki * muki;
			mui.push_back(tmp(0, 0));
		}

		Vi = calVar(mui);
		//Si.push_back(Vi / V);

		if (Opt == 'T')
		{
			Si.push_back(1 - Vi / V); // total
		}
		else
		{
			Si.push_back(Vi / V);   // main
		}

		printf("GSA i=%i, Si=%.2f, K=%i \n", nc + 1, Si[nc], Kos);

		if (isinf(Si[nc]) || isnan(Si[nc]))
		{
			return { -100 };
		}
	}

	return Si;
}


double runGSA::mvnPdf(mat x, mat mu, mat cov) 
{
	
	double n = size(x)(1);
	double sqrt2pi = std::sqrt(2 * PI);
	mat xmu = x - mu;
	mat quadform = xmu.t() * inv(cov) * xmu;
	double norm = std::pow(sqrt2pi, -n)*std::pow(abs(det(cov)), -0.5);
	//std::cout << norm << std::endl;

	return norm * std::exp(-0.5 * quadform(0,0));
}

double runGSA::calMean(vector<double> x) {
	double sum = std::accumulate(std::begin(x), std::end(x), 0.0);
	return sum / x.size();

}

runGSA::~runGSA() {};

double runGSA::calVar(vector<double> x) {
	double m = calMean(x);
	double accum = 0.0;
	std::for_each(std::begin(x), std::end(x), [&](const double d) {
		accum += (d - m) * (d - m);
		});
	//std::cout << (accum / (x.size())) << std::endl;
	return (accum / (x.size()));
}



double runGSA::writeOutputs(jsonInput inp)
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
			outfile << Simat[i][j] << " ";
		}
		for (int j = 0; j < inp.ngr; j++) {
			outfile << Stmat[i][j] << " ";
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