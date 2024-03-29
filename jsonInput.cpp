
#include "jsonInput.h"

jsonInput::jsonInput(string workDir)
{
	this->workDir = workDir;
	std::ifstream myfile(workDir + "/templatedir/dakota.json");
	if (!myfile.is_open()) {
		std::string errMsg = "Error running UQ engine: Unable to open JSON";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	json UQjson = json::parse(myfile, nullptr, false);
	if (UQjson.is_discarded())
	{
		std::string errMsg = "Error reading json: JSON syntax is broken";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	//json UQjson = json::parse(myfile);

	uqType = UQjson["UQ_Method"]["uqType"];

	if ((uqType.compare("Forward Propagation") == 0) || (uqType.compare("Sensitivity Analysis") == 0)) {
		// pass
	} else
	{
		//*ERROR*
	  std::string errMsg = "Error reading json: 'Forward Analysis' or 'Sensitivity Analysis' backend is called, but the user requested " + UQjson["UQ_Method"]["uqType"].get<std::string>();
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}

	nmc = UQjson["UQ_Method"]["samplingMethodData"]["samples"];
	rseed = UQjson["UQ_Method"]["samplingMethodData"]["seed"];
	UQmethod = UQjson["UQ_Method"]["samplingMethodData"]["method"];
	//
	// Specify parameters in each distributions.
	//

	//std::vector<int> corrIdx;
	std::vector<int> randIdx, constIdx, resampIdx;
	int count = 0;
	nrv = 0;
	nco = 0;
	nre = 0;

	std::string resampGroupTxt;
	if (UQjson["UQ_Method"].find("RVdataGroup") != UQjson["UQ_Method"].end()) {
		// if the key "sensitivityGroups" exists
		resampGroupTxt = UQjson["UQ_Method"]["RVdataGroup"];
		resampGroupTxt.erase(remove(resampGroupTxt.begin(), resampGroupTxt.end(), ' '), resampGroupTxt.end());
	} else {
		resampGroupTxt = "";
	}
	vector<vector<string>> resamplingGroupsString;
	vector<string> flattenResamplingGroups;

	fromTextToStr(resampGroupTxt, resamplingGroupsString, flattenResamplingGroups);
	nreg = resamplingGroupsString.size();

	auto it = std::unique(flattenResamplingGroups.begin(), flattenResamplingGroups.end());
	bool isUnique = (it == flattenResamplingGroups.end());
	if (!isUnique) {
		//*ERROR*
		std::string errMsg = "Error reading input: groups of random variables should be mutually exclusive";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);

	}

	for (auto& elem : UQjson["randomVariables"])
	{
		if (elem.find("inputType") == elem.end())
		{
			//*ERROR*
			std::string errMsg = "Error reading json: input file does not have the key 'inputType'";
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);
		}
		// if key "correlationMatrix" exists
		
		// name of distribution
		std::string distName = elem["distribution"];
		std::transform(distName.begin(), distName.end(), distName.begin(), ::tolower); // make it lower case
		distName.erase(remove_if(distName.begin(), distName.end(), isspace), distName.end()); // remove space

		// type of input (PAR, MOM, or DAT)
		std::string inpType = elem["inputType"];
		std::string inpTypeSub = inpType.substr(0, 3);
		for (int i = 0; i < 3; i++) {
			inpTypeSub[i] = toupper(inpTypeSub[i]);
		}

		// get parameter names for each dist
		std::vector<std::string> pnames;
		getPnames(distName, inpTypeSub, pnames);

		if (std::find(flattenResamplingGroups.begin(), flattenResamplingGroups.end(), elem["name"]) != flattenResamplingGroups.end()) {
			//is_resamplingGroup
			if (!(distName.compare("discrete") == 0) || !(inpTypeSub.compare("DAT")) == 0) {
				//*ERROR*
				string InputType;
				if (!inpTypeSub.compare("DAT")) {
					InputType = "Dataset";
				}
				else if (!inpTypeSub.compare("PAR")) {
					InputType = "Parameters";
				}
				else {
					InputType = "Moments";
				}
				std::string errMsg = "Error reading input: RVs specified in UQ tab should have the option <Dataset-Discrete>. Your input is <" + InputType + "-" + distName + ">";
				std::cout << errMsg << "\n";
				theErrorFile << errMsg << std::endl;
				theErrorFile.close();
				exit(-1);
			}
			
			resampIdx.push_back(count);
			nre++;
			count++;
			continue;
		}


		if (distName.compare("constant") == 0) {
			constIdx.push_back(count);
			nco++;
			count++;
			continue;
		}
		if ((distName.compare("discrete") == 0) && (inpTypeSub.compare("PAR")) == 0) {
			if (elem[pnames[0]].size() == 1) {
				// discrete distribution with only one quantity = constant
				constIdx.push_back(count);
				nco++;
				count++;
				continue;
			}
		}

		// save name of random variable etc
		rvNames.push_back(elem["name"]);
		distNames.push_back(distName);
		opts.push_back(inpTypeSub);

		// if "DAT" 
		if (opts[nrv].compare("DAT") == 0) {

			// Sample set inside vals
			//std::string directory = elem["dataDir"];
			//std::string tmpName = elem["name"];
			//std::string directory = ".\\templatedir\\"+ tmpName +".in";

			std::string tmpName = elem["name"];
			std::filesystem::path dir = workDir;
			//std::filesystem::path relPath = "templatedir\\" + tmpName + ".in";
			std::filesystem::path relPath = std::filesystem::path("templatedir") / std::filesystem::path(tmpName+ ".in") ;
			std::filesystem::path directory = dir / relPath;
			std::ifstream data_table(directory);
			if (!data_table.is_open()) {
				//*ERROR*
				std::string errMsg = "Error reading json: cannot open data file " + directory.u8string();
				std::cout << errMsg << "\n";
				theErrorFile << errMsg << std::endl;
				theErrorFile.close();
				exit(-1);
			}

			std::vector<double> vals_tmp;
			double samps = 0.0;
			while (data_table >> samps)
			{ 
				vals_tmp.push_back(samps);
				if (data_table.peek() == ',')
					data_table.ignore();
			}
			data_table.close();
			resampleCandidates.push_back(vals_tmp);
			vals.push_back(vals_tmp);

			if (vals_tmp.size() < 1) { //*ERROR*
				int a = vals_tmp.size();
				std::cout << a << std::endl;

				std::string errMsg = "Error reading json: data file of " + rvNames[nrv] + " has less then one sample.";
				std::cout << errMsg << "\n";
				theErrorFile << errMsg << std::endl;
				theErrorFile.close();
				exit(-1);
			}

			// Save boundary informaions
			if (distNames[nrv].compare("binomial") == 0) {
				adds.push_back({ elem["n"],0.0 }); // not used
			}
			else if (distNames[nrv].compare("beta") == 0) {
				adds.push_back({ elem["lowerbound"],elem["upperbound"] });
			}
			else if (distNames[nrv].compare("truncatedexponential") == 0) {
				adds.push_back({ elem["a"],elem["b"] });
			}
			else
			{
				adds.push_back({}); // default
			}
		}
		else // if "PAR" or "MOM" 
		{
			// Parameter (moment) values inside vals
			if (distNames[nrv].compare("discrete") == 0) {
				std::vector<double> vals_tmp;
				int numdisc = elem[pnames[0]].size();
				for (int i = 0; i < numdisc; i++)
				{
					vals_tmp.push_back(elem[pnames[0]][i]);
					vals_tmp.push_back(elem[pnames[1]][i]);
				}
				vals.push_back(vals_tmp);
			}
			else
			{
				std::vector<double> vals_temp;
				for (auto& pn : pnames)
				{
					if (elem.find(pn) != elem.end())
					{ 
						vals_temp.push_back(elem[pn]); // get parameter values
					}
					else
					{
						std::string errMsg = "Error reading json: cannot find <" + pn + "> in " + distName + " from input json.";
						std::cout << errMsg << "\n";
						theErrorFile << errMsg << std::endl;
						theErrorFile.close();
						exit(-1);
					}

				}
				vals.push_back(vals_temp);
			}
			adds.push_back({});
		}
		randIdx.push_back(count);
		nrv++;
		count++;
	}

	//
	// get resamples
	//


	for (int i : resampIdx)
	{

		auto elem = UQjson["randomVariables"][i];


		std::string tmpName = elem["name"];
		std::filesystem::path dir = workDir;
		//std::filesystem::path relPath = "templatedir\\" + tmpName + ".in";
		std::filesystem::path relPath = std::filesystem::path("templatedir") / std::filesystem::path(tmpName+ ".in") ;
		std::filesystem::path directory = dir / relPath;

		std::ifstream data_table(directory);
		if (!data_table.is_open()) {
			//*ERROR*
			std::string errMsg = "Error reading json: cannot open data file: " + directory.u8string();
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);
		}

		std::vector<double> vals_tmp;
		double samps = 0.0;
		while (data_table >> samps)
		{
			vals_tmp.push_back(samps);
			if (data_table.peek() == ',')
				data_table.ignore();
		}
		data_table.close();

		if (vals_tmp.size() < 1) { //*ERROR*
			int a = vals_tmp.size();
			std::cout << a << std::endl;

			//*ERROR*
			std::string errMsg = "Error reading json: data file of " + rvNames[nrv] + " has less then one sample.";
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);

		}

		vals.push_back(vals_tmp);
		resampleCandidates.push_back(vals_tmp);
		rvNames.push_back(elem["name"]);
	}

	//
	// get constants
	//
	
	//for (auto& elem : UQjson["randomVariables"])
	for (int i : constIdx)
	{
		// name of distribution
		auto elem = UQjson["randomVariables"][i];
		string distname = elem["distribution"];
		std::transform(distname.begin(), distname.end(), distname.begin(), ::tolower); // make lower case

		// input type
		std::string inpType = elem["inputType"];
		std::string inpTypeSub = inpType.substr(0, 3);
		std::transform(inpTypeSub.begin(), inpTypeSub.end(), inpTypeSub.begin(), ::toupper); // make upper case

		// parameter name
		std::vector<std::string> pnames;
		getPnames(distname, inpTypeSub, pnames);

		// If constant
		if (distname.compare("constant") == 0) 
		{
			// *name of random variable
			rvNames.push_back(elem["name"]);
			constants.push_back(elem[pnames[0]]);
		}
		// If constant (discrete)
		else if ((distname.compare("discrete") == 0) && (inpTypeSub.compare("PAR")) == 0) 
		{ 
			if (elem["value"].size() == 1) {
				// discrete distribution with only one quantity = constant
				rvNames.push_back(elem["name"]);
				constants.push_back(elem[pnames[0]][0]);
			}
		}
	}



	//
	// get edp names
	//

	nqoi = 0;
	for (auto& elem : UQjson["EDP"]) {
		// *name of distribution
		if (elem["length"] == 1) {
			qoiNames.push_back(elem["name"]);
			nqoi++;
		} else if (elem["length"] > 1) {
			std::string name = elem["name"];
			for (int j=0; j < elem["length"]; j++) {
				qoiNames.push_back(name + "_" + std::to_string(j+1));
				nqoi++;
			}
		}
	}

	//
	// get correlation matrix
	//

	//vector<vector<double>> corr;
	//for (int i = 0; i < nrv; i++)
	//{
	//	vector<double> corr_row(nrv, 0.2);
	//	corr_row[i] = 1.0; // overwrite diagonal elements
	//	corr.push_back(corr_row);
	//}

	//corr.reserve(nrv*nrv);
	
	if (UQjson.find("correlationMatrix") != UQjson.end()) {
		corr = *new vector<vector<double>>(nrv, vector<double>(nrv, 0.0));
		// if key "correlationMatrix" exists
		for (int i = 0; i < nrv; i++) {
			for (int j = 0; j < nrv; j++) {
				corr[i][j] = UQjson["correlationMatrix"][randIdx[i] + randIdx[j] * (nrv + nco)];
			}
		}
	}
	else
	{
		//corr.assign(nrv* nrv, 0);
		for (int i = 0; i < nrv; i++) {
			vector<double> corr_row(nrv, 0.0);
			corr_row[i] = 1.0; // overwrite diagonal elements
			corr.push_back(corr_row);
		}

//}
	}


	//
	// get resampling group index matrix
	//

	fromTextToId(resampGroupTxt, rvNames, resamplingGroups);

	//
	// get group index matrix
	//

	if (UQjson["UQ_Method"].find("sensitivityGroups") != UQjson["UQ_Method"].end()) {
		// if the key "sensitivityGroups" exists

		std::string groupTxt = UQjson["UQ_Method"]["sensitivityGroups"];
		groupTxt.erase(remove(groupTxt.begin(), groupTxt.end(), ' '), groupTxt.end()); // remove any white spaces
		fromTextToId(groupTxt, rvNames, groups);
	}
	else {
		for (int i = 0; i < nrv; i++) {
			groups.push_back({i});
		}
		std::cout << nreg <<std::endl;
		for (int i = 0; i < nreg; i++) {
			for (int j = 0; j < size(resamplingGroups[i]); j++) {
				groups.push_back({ resamplingGroups[i][j] });
			}
		}

	}
	ngr = groups.size();



	for (int i = 0; i < nreg; i++) {
		int length_old = size(vals[resamplingGroups[i][0]]);
		int length_data;
		for (int j = 1; j < size(resamplingGroups[i]); j++) {
			length_data = size(vals[resamplingGroups[i][j]]);
			if (length_data != length_old)
			{
				std::string errMsg = "Error reading json: RVs in the same group do not have the same number of samples";
				std::cout << errMsg << "\n";
				theErrorFile << errMsg << std::endl;
				theErrorFile.close();
				exit(-1);
			}
			length_old = length_data;
		}
		resamplingSize.push_back(length_data);
	}

}

void
jsonInput::fromTextToId(string groupTxt, vector<string>& groupPool, vector<vector<int>>& groupIdVect)
{
	int nrv = size(groupPool);
	std::regex re(R"(\{([^}]+)\})"); // will get string inside {}
	std::sregex_token_iterator it(groupTxt.begin(), groupTxt.end(), re, 1);
	std::sregex_token_iterator end;
	while (it != end) {
		std::stringstream ss(*it++);
		std::vector<int> groupID;
		std::vector<string> groupString;
		while (ss.good()) {
			std::string substr;
			getline(ss, substr, ',');  // incase we have multiple strings inside {}
			groupString.push_back(substr);

			std::vector<std::string>::iterator itr = std::find(groupPool.begin(), groupPool.end(), substr);
			if (itr != groupPool.cend()) { // from names (a,b,{a,b}) to idx's (1,2,{1,2})		
				int index_rvn = std::distance(groupPool.begin(), itr); // start from 0
				groupID.push_back((int)index_rvn);
				std::string errMsg;
				if (index_rvn >= nrv) {
					std::string errMsg = "Error reading json: RV group (for Sobol) cannot contain constant variable: " + groupPool[index_rvn]  ;
					std::cout << errMsg << "\n";
					theErrorFile << errMsg << std::endl;
					theErrorFile.close();
					exit(-1);
				}
			}
			else {
				// *ERROR*
				std::string errMsg = "Error reading json: element <" + substr + "> inside the variable groups not found.";
				std::cout << errMsg << "\n";
				theErrorFile << errMsg << std::endl;
				theErrorFile.close();
				exit(-1);

			}
		}
		groupIdVect.push_back(groupID);
	}
}

void
jsonInput::fromTextToStr(string groupTxt, vector<vector<string>>& groupStringVector, vector<string>& flattenStringVect)
{
	std::regex re(R"(\{([^}]+)\})"); // will get string inside {}
	std::sregex_token_iterator it(groupTxt.begin(), groupTxt.end(), re, 1);
	std::sregex_token_iterator end;
	while (it != end) {
		std::stringstream ss(*it++);
		std::vector<string> groupString;
		while (ss.good()) {
			std::string substr;
			getline(ss, substr, ',');  // incase we have multiple strings inside {}
			groupString.push_back(substr);
			flattenStringVect.push_back(substr);
		}
		groupStringVector.push_back(groupString);
	}
}

void 
jsonInput::getPnames(string distname, string optname, vector<std::string>& par_char)
{
	if (optname.compare("PAR") == 0) { // Get parameters

		if (distname.compare("binomial") == 0) {  // Not used
			par_char.push_back("n");
			par_char.push_back("p");
		}
		else if (distname.compare("geometric") == 0) {  // Not used
			par_char.push_back("p");
		}
		else if (distname.compare("negativebinomial") == 0) {  // Not used
			par_char.push_back("k");
			par_char.push_back("p");
		}
		else if (distname.compare("poisson") == 0) {
			par_char.push_back("lambda");
			par_char.push_back("t");
		}
		else if (distname.compare("uniform") == 0) {
			par_char.push_back("lowerbound");
			par_char.push_back("upperbound");
		}
		else if (distname.compare("normal") == 0) {
			par_char.push_back("mean");
			par_char.push_back("stdDev");
		}
		else if (distname.compare("lognormal") == 0) {
			par_char.push_back("lambda");
			par_char.push_back("zeta");
		}
		else if (distname.compare("exponential") == 0) {
			par_char.push_back("lambda");
		}
		else if (distname.compare("gamma") == 0) {
			par_char.push_back("k");
			par_char.push_back("lambda");
		}
		else if (distname.compare("beta") == 0) {
			par_char.push_back("alphas");
			par_char.push_back("betas");
			par_char.push_back("lowerbound");
			par_char.push_back("upperbound");
		}
		else if (distname.compare("gumbelMin") == 0) {  // Not used
			par_char.push_back("an");
			par_char.push_back("bn");
		}
		else if (distname.compare("gumbel") == 0) {
			par_char.push_back("alphaparam");
			par_char.push_back("betaparam");
		}
		else if (distname.compare("frechet") == 0) {  // Not used
			par_char.push_back("an");
			par_char.push_back("k");
		}
		else if (distname.compare("weibull") == 0) {
			par_char.push_back("scaleparam"); //an
			par_char.push_back("shapeparam"); //k
		}
		else if (distname.compare("gev") == 0) {  // Not used
			par_char.push_back("beta");
			par_char.push_back("alpha");
			par_char.push_back("epsilon");
		}
		else if (distname.compare("gevmin") == 0) {  // Not used
			par_char.push_back("beta");
			par_char.push_back("alpha");
			par_char.push_back("epsilon");
		}
		else if (distname.compare("pareto") == 0) {  // Not used
			par_char.push_back("x_m");
			par_char.push_back("alpha");
		}
		else if (distname.compare("rayleigh") == 0) {  // Not used
			par_char.push_back("alpha");
		}
		else if (distname.compare("chisquare") == 0) {
			par_char.push_back("k");
		}
		else if (distname.compare("discrete") == 0) {
			par_char.push_back("Values");
			par_char.push_back("Weights");
		}
		else if (distname.compare("truncatedexponential") == 0) {
			par_char.push_back("lambda");
			par_char.push_back("a");
			par_char.push_back("b");
		}
		else if (distname.compare("constant") == 0) {
			par_char.push_back("value");
		}
		else {
			std::string errMsg = "Error reading json: cannot interpret distribution name: " + distname;
			std::cout << errMsg << "\n";
			theErrorFile << errMsg << std::endl;
			theErrorFile.close();
			exit(-1);
		}
	}
	else if (optname.compare("MOM") == 0) { // Get Moments
		if (distname.compare("normal") == 0) {  // Not used
			par_char.push_back("mean");
			par_char.push_back("stdDev");  // 
		}
		else if (distname.compare("lognormal") == 0) {  // Not used
			par_char.push_back("mean");
			par_char.push_back("stdDev");  // 
		}
		else if (distname.compare("geometric") == 0) {  // Not used
			par_char.push_back("mean");
		}
		else if (distname.compare("poisson") == 0) {
			par_char.push_back("mean");
		}
		else if (distname.compare("exponential") == 0) {
			par_char.push_back("mean");
		}
		else if (distname.compare("beta") == 0) {
			par_char.push_back("mean");
			par_char.push_back("standardDev");
			par_char.push_back("lowerbound");
			par_char.push_back("upperbound");
		}
		else if (distname.compare("gev") == 0) {
			par_char.push_back("mean");
			par_char.push_back("standardDev");
			par_char.push_back("epsilon");
		}
		else if (distname.compare("gevmin") == 0) {  // Not used
			par_char.push_back("mean");
			par_char.push_back("standardDev");
			par_char.push_back("epsilon");
		}
		else if (distname.compare("rayleigh") == 0) {  // Not used
			par_char.push_back("mean");
		}
		else if (distname.compare("chisquare") == 0) {
			par_char.push_back("mean");
		}
		else if (distname.compare("constant") == 0) {
			par_char.push_back("value");
		}
		else if (distname.compare("truncatedexponential") == 0) {
			par_char.push_back("mean");
			par_char.push_back("a");
			par_char.push_back("b");
		}
		else
		{
			par_char.push_back("mean");
			par_char.push_back("standardDev");
		}
	}
	else if (optname.compare("DAT") == 0) { // Get DATA	
		if (distname.compare("binomial") == 0) {
			par_char.push_back("n");
		}
		else if (distname.compare("beta") == 0) {
			par_char.push_back("lowerbound");
			par_char.push_back("upperbound");
		}
		else if (distname.compare("truncatedexponential") == 0) {
			par_char.push_back("a");
			par_char.push_back("b");
		}
		else if (distname.compare("constant") == 0) {
			par_char.push_back("value");
		}
	}
	else {
		std::string errMsg = "Error reading json: input type should be one of PAR/MOM/DAT";
		std::cout << errMsg << "\n";
		theErrorFile << errMsg << std::endl;
		theErrorFile.close();
		exit(-1);
	}
}

jsonInput::~jsonInput(void) {}
