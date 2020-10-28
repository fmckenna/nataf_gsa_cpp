
#include "jsonInput.h"

jsonInput::jsonInput(string workDir)
{
	this->workDir = workDir;
	std::ifstream myfile(workDir + "/templatedir/dakota.json");
	if (!myfile.is_open()) {
		theErrorFile << "Error running UQ engine: Unable to open dakota.json";
		theErrorFile.close();
		exit(-1);
	}

	json UQjson = json::parse(myfile);

	// 
	// Get variables
	// 

	nmc = UQjson["UQ_Method"]["samplingMethodData"]["samples"];
	rseed = UQjson["UQ_Method"]["samplingMethodData"]["seed"];
	UQmethod = UQjson["UQ_Method"]["samplingMethodData"]["method"];

	//
	// Specify parameters in each distributions.
	//

	//std::vector<int> corrIdx;
	std::vector<int> randIdx, constIdx;
	int count = 0;
	nrv = 0;
	nco = 0;
	for (auto& elem : UQjson["randomVariables"])
	{
		if (elem.find("inputType") == elem.end())
		{
			//*ERROR*
			theErrorFile << "Error reading json: input file does not have the key 'inputType'" << std::endl;
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
			std::string directory = elem["dataDir"];
			std::ifstream data_table(directory);
			if (!data_table.is_open()) {
				//*ERROR*
				theErrorFile << "Error reading json: cannot open data file at " << directory << std::endl;
				theErrorFile.close();
				exit(-1);
			}

			std::vector<double> vals_tmp;
			double samps = 0.0;
			while (data_table >> samps)
			{ 
				vals_tmp.push_back(samps);
			}
			data_table.close();
			vals.push_back(vals_tmp);

			if (vals_tmp.size() < 3) { //*ERROR*
				theErrorFile << "Error reading json: data file of " << rvNames[nrv] << " has less then three samples." << std::endl;
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
					vals_temp.push_back(elem[pn]); // get parameter values
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
		qoiNames.push_back(elem["name"]);
		nqoi++;
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
	// get group index matrix
	//

	if (UQjson["UQ_Method"].find("sensitivityGroups") != UQjson["UQ_Method"].end()) {
		// if the key "sensitivityGroups" exists

		std::string groupTxt = UQjson["UQ_Method"]["sensitivityGroups"];
		std::regex re(R"(\{([^}]+)\})"); // will get string inside {}
		std::sregex_token_iterator it(groupTxt.begin(), groupTxt.end(), re, 1);
		std::sregex_token_iterator end;
		while (it != end) {
			std::stringstream ss(*it++);
			std::vector<int> aGroup; 
			while (ss.good()) {
				std::string substr;
				getline(ss, substr, ',');  // incase we have multiple strings inside {}
				std::vector<std::string>::iterator itr = std::find(rvNames.begin(), rvNames.end(), substr);
				if (itr != rvNames.cend()) { // from names (a,b,{a,b}) to idx's (1,2,{1,2})		
					int index_rvn = std::distance(rvNames.begin(), itr); // start from 0
					aGroup.push_back((int)index_rvn);

					if (index_rvn > nrv) {
						// If it is a constant variable
						theErrorFile << "Error reading json: RV group (for Sobol) cannot contain constant variable: ";
						theErrorFile << rvNames[index_rvn - 1] << std::endl;
						theErrorFile.close();
						exit(-1);
					}
				}
				else {
					// *ERROR*
					theErrorFile << "Error reading json: element <" << substr << "> inside the sensitivity groups not found." << std::endl;
					theErrorFile.close();
					exit(-1);
				}
			}
			groups.push_back(aGroup);
		}
	}
	else {
		for (int i = 0; i < nrv; i++) {
			groups.push_back({i});
		}
	}
	ngr = groups.size();

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
			theErrorFile << "Error reading json: cannot interpret distribution name: " << distname;
			theErrorFile.close();
			exit(-1);
			// NA
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
		theErrorFile << "Error reading json: input type should be one of PAR/MOM/DAT";
		theErrorFile.close();
		exit(-1);
	}
}
jsonInput::~jsonInput(void) {}