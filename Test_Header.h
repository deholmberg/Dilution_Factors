// The purpose of this header file is to contain extra information and functions to make my 
// "First_Pass_Asymmetry.C" file less cluttered.

using namespace std;

// This struct holds information about each of the kinematic bins used in the analysis
struct Bin{

	// This holds the bin boundary information, including midpoints
	double Q2_Min = 0.0; double Q2_Max = 0.0; double Q2_Mid = 0.0;
	double X_Min = 0.0;  double X_Max = 0.0;  double X_Mid = 0.0;

	void SetQ2Bins( double Q2min, double Q2max ){
		Q2_Min = Q2min; Q2_Max = Q2max; Q2_Mid = (Q2min+Q2max)/2.0;
	}	
	void SetXBins( double Xmin, double Xmax ){
		X_Min = Xmin; X_Max = Xmax; X_Mid = (Xmin+Xmax)/2.0;
	}

	string Target_Type = "NH3"; // Kept as this for now

	// This holds the dilution factor information
	double Df = 0.0; double Err_Df = 0.0;
  // This holds the DF information using the PF in the calculation
  double Df2= 0.0; double Err_Df2= 0.0;

	// This holds the packing fraction information
	double Pf = 0.0; double Err_Pf = 0.0;

	void SetDf( double df, double errdf ){
		Df = df; Err_Df = errdf;
	}
 
  void SetDf2( double df2, double errdf2 ){
    Df2 = df2; Err_Df2 = errdf2;
  }

	void SetPf( double pf, double errpf ){
		Pf = pf; Err_Pf = errpf;
	}

	// This holds the value of the structure function ratio R for this bin
	double R = 0.0; double Err_R = 0.0;

	void SetR( double r, double errr){
		R = r; Err_R = errr;
	}
	// This holds kinematic variables that exist for each bin. These come from the
	// CLAS6 parameterizations and models
	double Eta = 0.0; double Depol_Factor = 0.0; double A2 = 0.0;
	void SetKinematics( double eta, double depol_factor, double a2 ){
		Eta = eta; Depol_Factor = depol_factor; A2 = a2;
	}

	// This holds the number of counts for the Faraday Cup-normalized counts.
	// NOTE: you MUST already include the FC normalization BEFORE you plug
	// values into here!!!
	// Assume electron helicity is + for pointing downstream, - points upstream
	// Assume proton helicity is + for pointing upstream, - points downstream
	// Helicity states given as (electron, proton)
	// n- = spin anti-aligned states (+,+), (-,-)
	// n+ = spin parallel states (+,-) ,(-,+)
	double n_m = 0.0; double n_p = 0.0;
	double Nm = 0; double Np = 0; // These hold the total counts for each helicity state;
				      // used for error propagation
	double N_total = 0.0;
	double FC_p = 0.0; double FC_m = 0.0; // Tracks total FC charge for each state
	double FC_total = 0.0;
	void AddNormalizedCounts( double nm, double np ){
		n_m += nm; n_p += np;
		//FC_p += 1.0/nm; FC_m += 1.0/np;
		//Nm++; Np++;
	}
	void AddCounts(double nm, double np){
		Nm += nm; Np += np;
	}
	void SetNormalizedCounts(){
		n_m = Nm / FC_m; n_p = Np / FC_p;
	}
 
  // This stuff holds the total number of counts for different target types.
  // For each target type, the FC charge and DIS electron count is stored.
  double nA = 0.0; double nCH = 0.0; double nC = 0.0; double nET = 0.0; double nF = 0.0; // raw counts
  double fA = 0.0; double fCH = 0.0; double fC = 0.0; double fET = 0.0; double fF = 0.0; // FC charges
  void SetAllCounts(double NA, double NCH, double NC, double NET, double NF){
    nA = NA; nCH = NCH; nC = NC; nET = NET; nF = NF;
  }
  void SetAllFCCharges(double FA, double FCH, double FC, double FET, double FF){
    fA = FA; fCH = FCH; fC = FC; fET = FET; fF = FF;
  }
  void AccumulateCounts(double& NA, double& NCH, double& NC, double& NET, double& NF){
    NA += nA; NCH += nCH; NC += nC; NET += nET; NF += nF;
  }
  void AccumulateCharges(double& FA, double& FCH, double& FC, double& FET, double& FF){
    FA += fA; FCH += fCH; FC += fC; FET += fET; FF += fF;
  }
	
	// This seems like a messy way to do things, and it is. BUT, here's my apologetic:
	// There's a FC charge for the beam's +/- states. However, these will be binned according
	// to the alignment of the beam helicity state relative to that of the target. By accumulating
	// the FC charge in this way, I can account for this alignment.
	//void SetThisRunFC(double thisFCp, double thisFCm ){
	//	FC_p += thisFCp; FC_m += thisFCm;
	//}
	double A_ll = 0.0; double Err_A_ll = 0.0; double A1_val = 0.0; double Err_A1_val = 0.0;
	double Asymmetry(double pbpt){
		if( n_m != 0.0 || n_p != 0.0 ){
			double asym = ((n_m - n_p)/(n_m + n_p)) / (pbpt * Df);
			A_ll = asym;
			return asym;
		}
		else return -999; // Error case
	}
	double Asymmetry_Error(double pbpt, double pbpt_err){
		if( n_m != 0.0 || n_p != 0.0 || FC_p != 0.0 || FC_m != 0.0 ){
			// This calculates the propagation of the dilution factor uncertainty
			double df_numerator = pow((n_p*n_p - n_m*n_m),2)*Err_Df*Err_Df;
			double pbpt_num = pow((n_p*n_p - n_m*n_m),2)*pbpt_err*pbpt_err;
			double asym_err = abs(2.0/(pow((n_m+n_p),2)*pbpt*Df )) * sqrt( (n_m*n_p*n_p)/FC_m + (n_p*n_m*n_m)/FC_p + df_numerator/(Df*Df) + pbpt_num/(pbpt*pbpt) );
			Err_A_ll = asym_err;
			return asym_err;
		}
		else return -999; // Error case
	}
	double RawAsymmetry(){
    SetNormalizedCounts();
		if( n_m != 0.0 || n_p != 0.0 ){
			double asym = ((n_m - n_p)/(n_m + n_p));
			//A_ll = asym;
			return asym;
		}
		else return -999; // Error case
	}
	double RawAsymmetryError(){
    SetNormalizedCounts();
		if( n_m != 0.0 || n_p != 0.0 || FC_p != 0.0 || FC_m != 0.0 ){
			// This calculates the propagation of the dilution factor uncertainty
			double asym_err = abs(2.0/(pow((n_m+n_p),2) )) * sqrt( (n_m*n_p*n_p)/FC_m + (n_p*n_m*n_m)/FC_p );
			//Err_A_ll = asym_err;
			return asym_err;
		}
		else return -999; // Error case
	}


	// Returns the value for A1
	double A1(double asym){
		double a1 = (asym / Depol_Factor) - (Eta*A2);
		double a1_err = Err_A_ll / Depol_Factor;
		A1_val = a1; Err_A1_val = a1_err;
		return a1;
	}

	string WriteOutput(){
		stringstream s;
		s << Q2_Min <<"	"<< Q2_Max<<"	"<<X_Min<<"	"<<X_Max<<"	"<<Eta<<"	"<<Depol_Factor<<"	"<<Df<<"	"<<Err_Df<<"	"<<A2<<"	";
		s << n_m <<"	"<< n_p <<"	"<< A_ll <<"	"<< Err_A_ll<<"	"<<A1_val<<"	"<<Err_A1_val<<"	";
		string out = s.str();
		return out;
	}
	string WriteBins(){
		stringstream s;
		s << Q2_Min <<"	"<< Q2_Max<<"	"<<X_Min<<"	"<<X_Max<<"	"<< (Nm+Np) <<"	"<< (FC_p+FC_m);
		string out = s.str();
		return out;
	}
	void PrintBinInfo(){
		cout << "=========================\n";
		cout << "Information for this bin:\n";
		cout << "Q^2 = "<< Q2_Mid << " GeV^2\n";
		cout << "x = " << X_Mid << endl;
		cout << "A1 = "<< A1_val << " +- "<< Err_A1_val << endl;
		cout << "=========================\n";
	}

};


// This struct holds information about each of the runs in the dataset, mainly
// to keep track of the Faraday Cup information for each run
struct Run{

	int Run_Number = 0;
	string Target_Type = "NONE";

	// This contains the gated Faraday Cup charge for the run
	double FC_Gated_Pos = 0.0; double FC_Gated_Neg = 0.0;

	void SetGatedFC( double fcpos, double fcneg ){
		FC_Gated_Pos = fcpos; FC_Gated_Neg = fcneg;
	}

	// This has the information on NMR target polarization and HWP status
	double Target_Pol_NMR = 0.0; int insertedHWP = 0;


	void SetSignInfo( double nmr_tpol, int hwp ){
		Target_Pol_NMR = nmr_tpol;
		if( hwp == 1 ) insertedHWP = -1;
		else if( hwp == 0 ) insertedHWP = 1;
		else cout << "Couldn't set HWP position; check info; HWP status is set to 0\n";
	}

};

// This function returns a vector and returns "Run" structs for each
// of the runs that have FCup information
vector<Run> ReadFCupInfo(){
    string line;
    vector<Run> All_Runs;
    ifstream fcin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/FCup_16000_17000.txt");
    while( getline( fcin, line ) ){
	stringstream sin(line);
	int runnum; double fcpos, fcneg, fczero, fctotal;
	sin >> runnum >> fcpos >> fcneg >> fczero >> fctotal;
	Run thisRun; thisRun.Run_Number = runnum;
	thisRun.SetGatedFC( fcpos, fcneg );
	All_Runs.push_back( thisRun ); // Adds the run with FC info to the list	
    }
    fcin.close();
    return All_Runs;
}

// This reads in the RCDB information for each run and assigns it to each run in All_Runs vector
// Run_Number  Target  Beam_Cur_Req  Beam_Energy  Num_Events  HWP_Status  Target_Pol  Solenoid_Scale  Torus_Scale
void SetRCDBInfo( vector<Run>& All_Runs ){
    string line;
    ifstream rcdbin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/RCDB_Run_Data.txt");
    getline(rcdbin,line); // Throw away header row
    while( getline( rcdbin, line ) ){
	int runnum, HWPstat; string target; double beamcur, beamE, numEvents, tpol, sol_scale, tor_scale;
	stringstream sin(line);
	sin >> runnum >> target >> beamcur >> beamE >> numEvents >> HWPstat >> tpol >> sol_scale >> tor_scale;
	// This is super inefficient, but I'll fix it later
	for(size_t i=0; i<All_Runs.size(); i++){
	    if( All_Runs[i].Run_Number == runnum ){
		All_Runs[i].SetSignInfo( tpol, HWPstat );
		All_Runs[i].Target_Type = target;
		i = All_Runs.size(); // kills loop
	    }
	}	
    }
    rcdbin.close();
}


// Ok, this is really inefficient and dumb, but I loop through "updated_info.txt" twice.
// First, I go through to get the boundaries of the x, Q2 bins so I can make the array of bins easier to use.
vector< vector<Bin>> SetKinematicBins(){
    string line;  
    //ifstream sortin("Input_Text_Files/updated_info.txt");
    ifstream sortin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/myDFinfo.txt");
    vector<double> xbins, q2bins; // Holds the bin boundaries given in the rgc file
    // This loop reads in from the file, but only saves all of the unique values of Q2 and x so I 
    // can sort them afterwards.
    double xBinWidth = 0.0125; // The total size of an x bin is x_mid +/- xBinWidth
    getline(sortin,line); // Throw out header line
	while(getline(sortin,line)){
        //xbin, Q2bin, W2, Q2, x, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2
		double junk1, junk2, w2, q2_mid, x_mid, q2min, q2max;
		stringstream sin(line); sin >> junk1 >> junk2 >> w2 >> q2_mid >> x_mid >> q2min >> q2max;
		double xmin = x_mid - xBinWidth; double xmax = x_mid + xBinWidth;
		if( xbins.size() != 0 && q2bins.size() != 0){
			bool minIsInsideQ2Bin = false;
			bool minIsInsideXBin  = false;
			bool maxIsInsideQ2Bin = false;
			bool maxIsInsideXBin  = false;
			for(size_t i=0; i<xbins.size(); i++){
				if(xbins[i] == xmin) minIsInsideXBin = true;
				if(xbins[i] == xmax) maxIsInsideXBin = true;
			}
			for(size_t i=0; i<q2bins.size(); i++){
				if(q2bins[i] == q2min) minIsInsideQ2Bin = true;
				if(q2bins[i] == q2max) maxIsInsideQ2Bin = true;
			}
			if(!minIsInsideXBin) xbins.push_back(xmin);
			if(!maxIsInsideXBin) xbins.push_back(xmax);
			if(!minIsInsideQ2Bin)q2bins.push_back(q2min);
			if(!maxIsInsideQ2Bin)q2bins.push_back(q2max);
		}
		else{
			xbins.push_back(xmin); xbins.push_back(xmax);
			q2bins.push_back(q2min); q2bins.push_back(q2max);		
		}
	}
    sortin.close();
    // Now sort the values:
    sort( xbins.begin(), xbins.end() );
    sort( q2bins.begin(), q2bins.end() );
    // Define the number of bins in x, Q2 that are to be used
    int num_XBins = xbins.size()-1;
    int num_Q2Bins = q2bins.size()-1;
    // Create a list of all the mid values for each of the bins
    vector<string> q2mid; // implemented later in the loop

    // Now, I read in the values for the NH3 dilution factors
    vector< vector<Bin> > All_Bins( abs(num_Q2Bins), vector<Bin>(abs(num_XBins)) );
    //ifstream dfin("Input_Text_Files/updated_info.txt");
    ifstream dfin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/myDFinfo.txt");
    getline(dfin, line); // Once again, throw away the header line
    while( getline( dfin, line ) ){
	stringstream sin(line);
        //xbin, Q2bin, W2, Q2, x, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2
	double junk1, junk2, w2, q2_mid, x_mid, q2min, q2max, df_nh3, err_df_nh3;
	double q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2;
	sin >> junk1 >> junk2 >> w2 >> q2_mid >> x_mid >> q2min >> q2max >> df_nh3 >> err_df_nh3 >> q2_avg >> W >> nu >> y >> Eprime >> theta >> eps >> D >> eta >> F1 >> F2 >> R >> A1 >> A2 >> g1 >> g2;
	double xmin = x_mid - xBinWidth; double xmax = x_mid + xBinWidth;
	Bin thisBin;
	thisBin.SetQ2Bins( q2min, q2max ); thisBin.SetXBins( xmin, xmax );
	thisBin.SetDf( df_nh3, err_df_nh3 );
	thisBin.SetKinematics( eta, D, A2 );
	// Now take this bin and slot it into the appropriate place in All_Bins
	for(int i=0; i<num_Q2Bins; i++){
	    for(int j=0; j<num_XBins; j++){
		if( thisBin.X_Mid > xbins[j] && thisBin.X_Mid < xbins[j+1] &&
		    thisBin.Q2_Mid > q2bins[i] && thisBin.Q2_Mid < q2bins[i+1])
			All_Bins[i][j] = thisBin;
	    }
	}
    }
    dfin.close();
    return All_Bins;
}

// This goes through the same as above, but it also sets the A1 value from what's in "updated_info.txt"
vector< vector<Bin>> SetA1Bins(){
    string line;  
    //ifstream sortin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/updated_info.txt");
    //ifstream sortin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/old_data.txt");
    ifstream sortin("Input_Text_Files/test.txt");
    vector<double> xbins, q2bins; // Holds the bin boundaries given in the rgc file
    // This loop reads in from the file, but only saves all of the unique values of Q2 and x so I 
    // can sort them afterwards.
    double xBinWidth = 0.0125; // The total size of an x bin is x_mid +/- xBinWidth
    getline(sortin,line); // Throw out header line
	while(getline(sortin,line)){
        //xbin, Q2bin, W2, Q2, x, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2
		double junk1, junk2, w2, q2_mid, x_mid, q2min, q2max;
		stringstream sin(line); sin >> junk1 >> junk2 >> w2 >> q2_mid >> x_mid >> q2min >> q2max;
		double xmin = x_mid - xBinWidth; double xmax = x_mid + xBinWidth;
		if( xbins.size() != 0 && q2bins.size() != 0){
			bool minIsInsideQ2Bin = false;
			bool minIsInsideXBin  = false;
			bool maxIsInsideQ2Bin = false;
			bool maxIsInsideXBin  = false;
			for(size_t i=0; i<xbins.size(); i++){
				if(xbins[i] == xmin) minIsInsideXBin = true;
				if(xbins[i] == xmax) maxIsInsideXBin = true;
			}
			for(size_t i=0; i<q2bins.size(); i++){
				if(q2bins[i] == q2min) minIsInsideQ2Bin = true;
				if(q2bins[i] == q2max) maxIsInsideQ2Bin = true;
			}
			if(!minIsInsideXBin) xbins.push_back(xmin);
			if(!maxIsInsideXBin) xbins.push_back(xmax);
			if(!minIsInsideQ2Bin)q2bins.push_back(q2min);
			if(!maxIsInsideQ2Bin)q2bins.push_back(q2max);
		}
		else{
			xbins.push_back(xmin); xbins.push_back(xmax);
			q2bins.push_back(q2min); q2bins.push_back(q2max);		
		}
	}
    sortin.close();
    // Now sort the values:
    sort( xbins.begin(), xbins.end() );
    sort( q2bins.begin(), q2bins.end() );
    // Define the number of bins in x, Q2 that are to be used
    int num_XBins = xbins.size()-1;
    int num_Q2Bins = q2bins.size()-1;
    // Create a list of all the mid values for each of the bins
    vector<string> q2mid; // implemented later in the loop

    // Now, I read in the values for the NH3 dilution factors
    vector< vector<Bin> > All_Bins( abs(num_Q2Bins), vector<Bin>(abs(num_XBins)) );
    //ifstream dfin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/updated_info.txt");
    //ifstream dfin("/w/hallb-scshelf2102/clas12/holmberg/DF_Calculation/Input_Text_Files/old_data.txt");
    ifstream dfin("Input_Text_Files/test.txt");
    getline(dfin, line); // Once again, throw away the header line
    while( getline( dfin, line ) ){
	stringstream sin(line);
        //xbin, Q2bin, W2, Q2, x, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2
	double junk1, junk2, w2, q2_mid, x_mid, q2min, q2max, df_nh3, err_df_nh3;
	double q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2;
	sin >> junk1 >> junk2 >> w2 >> q2_mid >> x_mid >> q2min >> q2max >> df_nh3 >> err_df_nh3 >> q2_avg >> W >> nu >> y >> Eprime >> theta >> eps >> D >> eta >> F1 >> F2 >> R >> A1 >> A2 >> g1 >> g2;
	double xmin = x_mid - xBinWidth; double xmax = x_mid + xBinWidth;
	Bin thisBin;
	thisBin.SetQ2Bins( q2min, q2max ); thisBin.SetXBins( xmin, xmax );
	thisBin.A1_val = A1;
	thisBin.SetDf( df_nh3, err_df_nh3 );
	thisBin.SetKinematics( eta, D, A2 );
	// Now take this bin and slot it into the appropriate place in All_Bins
	for(int i=0; i<num_Q2Bins; i++){
	    for(int j=0; j<num_XBins; j++){
		if( thisBin.X_Mid > xbins[j] && thisBin.X_Mid < xbins[j+1] &&
		    thisBin.Q2_Mid > q2bins[i] && thisBin.Q2_Mid < q2bins[i+1])
			All_Bins[i][j] = thisBin;
	    }
	}
    }
    dfin.close();
    return All_Bins;
}

