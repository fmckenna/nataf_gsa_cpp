#pragma once

#include "lib_json/json.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <regex>
extern std::ofstream theErrorFile; // Error log

using json = nlohmann::json;
using std::string;
using std::vector;

class jsonInput
{
public:
	jsonInput(string workDir);
	virtual ~jsonInput(void);

	string workDir;
	int nmc;
	int nrv;
	int nco;
	int nre;
	int nqoi;
	int rseed;
	int ngr;
	string UQmethod;
	vector<string> distNames;
	vector<vector<double>> vals;
	vector<vector<double>> resampleCandidates;
	vector<double> constants;
	vector<string> opts;
	vector<string> rvNames;
	vector<string> qoiNames;
	vector<vector<double>> corr;
	vector<vector<double>> adds;
	vector<vector<int>> groups;
	vector<vector<int>> resamplingGroups;

private:
	void getPnames(string distname, string optname, vector<std::string>& par_char);
	void fromTextToVects(string groupTxt, vector<string> groupDomain, vector<vector<int>>& groupVect, vector<int>& flattenVect);

};

