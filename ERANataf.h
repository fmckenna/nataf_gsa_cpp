#pragma once
#include "ERADist.h"
#include "jsonInput.h"

using std::string;
using std::vector;
<<<<<<< HEAD
#include "Eigen/Dense"
=======
#include "lib_eigen/Eigen/Dense"
>>>>>>> upstream/master


extern std::ofstream theErrorFile; // Error log

class ERANataf
{

public:
	ERANataf(void);
	ERANataf(jsonInput inp);
	~ERANataf();

	int nrv;
	vector<vector<double>> Rhox;
	vector<vector<double>> Rhoz;
	Eigen::MatrixXd RhozMat, RhozInv;

	vector<vector<double>> X2U(int nmc, vector<vector<double>> x);
	vector<vector<double>> U2X(int nmc, vector<vector<double>> u);
	//ERADist **M_;
	vector<ERADist> M;
	void simulateAppBatch(string osType, 
						 string runType, 
						 jsonInput inp, 
					 	 vector<vector<double>> u,
						 vector<vector<int>> idx,
						 vector<vector<double>> &xvals, 
						 vector<vector<double>> &gvals);
	void simulateAppSequential(string osType,
						string runType,
						jsonInput inp,
						vector<vector<double>> u,
						vector<vector<double>>& xvals,
						vector<vector<double>>& gvals,
						int idx);

private:
	const double PI = 3.1415926535897932384626433;
	void quadGL(int N, double a, double b, vector<double>& x, vector<double>& w);
	Eigen::MatrixXd A;
	normal stdNorm;
	double mvnpdfR(Eigen::VectorXd u);
	double mvncdfR(Eigen::VectorXd u);
	double getJointPdf(vector<double> x);
	double getJointCdf(vector<double> x);
	double normCdf(double x);

};

// For MLE optimization
typedef struct NatafStr {
	vector<double> points = {};
	vector<double> fxii = {};
	vector<double> fxij = {};
	double Rhoxij=0;
	int ngrid=0;
} my_NatafInfo;
