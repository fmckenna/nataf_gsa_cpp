#include "runGSA.h"

#include <iterator>

runGSA::runGSA() {}

runGSA::runGSA(vector<vector<double>> xval,
	vector<vector<double>> gmat,
	vector<vector<int>> combs_tmp,
	int Kos)
{
	this->xval = xval;
	this->combs_tmp = combs_tmp;
	nmc = xval.size();
	this->Kos = std::min(Kos, int(nmc / 5));

	nrv = xval[0].size();
	ncombs = combs_tmp.size(); 
	int nqoi = gmat[0].size();

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

		vector<double> Stj = doGSA(gvec, 'M');
		vector<double> Sij = doGSA(gvec, 'T');

		vector<double> Si_temp, St_temp;


		Simat.push_back(Stj);
		Stmat.push_back(Sij);
	}
}

vector<double> runGSA::doGSA(vector<double> gval,char Opt)
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
				
				if ((isinf(pik_tmp(k))) || (isnan(pik_tmp(k))))
				{
					pik_tmp(k) = 1.e10;
					continue;
				}
				//double a = mvnPdf(xi, muval, cov.subcube(0, 0, k, endx, endx, k));
				//mat v = cov.subcube(0, 0, k, endx, endx, k);
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