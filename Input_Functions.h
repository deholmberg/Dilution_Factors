// This header file contains the information required to make the bins for
// the dilution factors in x, Q2.

using namespace std;

string THIS_DIR = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Input_Text_Files/";

// NOTE: NIST values are given to the precision of their listed uncertainty on the NIST website
const double p_mass = 0.93827208816; // proton mass in GeV/c^2 from NIST
const double n_mass = 0.93956542052; // neutron mass in GeV/c^2 from NIST
const double ppm_mass=0.13957039; // +-pion mass in GeV/c^2 from Particle Data Group
const double p0_mass= 4.5936/1000.0; // neutral pion mass in GeV/c^2 from Particle Data Group
const double k_mass = 0.493677; // K+- mass in GeV/c^2 from Particle Data Group
const double e_mass = 0.51099895/1000.0; // electron mass in GeV/c^2 from NIST
const double d_mass = 1.87561294257; // deuteron mass in GeV/c^2 from NIST
const double beam_energy = 10.5473; // electron beam energy in GeV

// These are parameters used in the calculations of the dilution factors and packing fraction
/*
const double lC = 1.678; // cm (length of carbon target)
const double LHe = 5.86; // cm (length of LHe bath)
const double lCH = 3.18; // cm (length of CH2 target)
const double lCD = 2.686;// cm (length of CD2 target)
const double pC = 1.7926/(12.0); // mol/cm^3 (carbon target density)
const double pCH = 0.9425/(14.027); // mol/cm^3 (CH2 target density)
const double pCD = 1.0979/(16.039); // mol/cm^3 (CD2 target density)
const double pA = 0.867/(17.031); // mol/cm^3 (NH3 target density)
const double pD = 1.007/(20.048); // mol/cm^3 (ND3 target density)
*/
// These values are the scaled lengths and densities

const double dLc= 0.005; // Contraction for carbon of 0.5%
const double dL = 0.021; // The percentage of the length for CH2 and CD2 contraction, at 2.1%
const double lC = 1.678*(1.0-dLc); // cm (length of carbon target)
const double LHe = 5.86*(1.0-0.018); // cm (length of LHe bath)
const double lCH = 3.18*(1.0-dL); // cm (length of CH2 target)
const double lCD = 2.686*(1.0-dL);// cm (length of CD2 target)
const double pC = 1.7926/(12.0*pow((1.0-dLc),3)); // mol/cm^3 (carbon target density)
const double pCH = 0.9425/(14.027*pow((1.0-dL),3)); // mol/cm^3 (CH2 target density)
const double pCD = 1.0979/(16.039*pow((1.0-dL),3)); // mol/cm^3 (CD2 target density)
const double pA = 0.867/(17.031); // mol/cm^3 (NH3 target density)
const double pD = 1.007/(20.048); // mol/cm^3 (ND3 target density)


// This first set of terms holds the constant info for the NH3 dilution factors
// Terms to hold some stuff so it's not so cluttered
double a = pC*lC*LHe + (pCH-pC)*lCH*lC - pCH*lCH*LHe;
double b = 9.0*pA*pC - 2.0*pCH*(pA+3.0*pC);
double c = 9.0*lC*pC*LHe*pA - 2.0*lCH*pCH*LHe*pA - lC*lCH*b;
// These hold the same values as above, only for the ND3 dilution factors
double aa = pC*lC*LHe + (pCD-pC)*lCD*lC - pCD*lCD*LHe;
double bb = 9.0*pD*pC - 2.0*pCD*(pD+3.0*pC);
double cc = 9.0*lC*pC*LHe*pD - 2.0*lCD*pCD*LHe*pD - lC*lCD*bb;

// Numerator terms for NH3 calculation
double uCH = abs( -pC*lC*LHe/a );
double uCB = abs( pCH*lCH*LHe/a );
double uFO = abs( lC*lCH*(pC-pCH)/a );
// Denominator terms for NH3 calculation
double dCH = abs( 9.0*lC*pC*LHe*pA/c );
double dCB = abs( -2.0*lCH*pCH*LHe*pA/c );
double dFO = abs( -b*lC*lCH/c );
// Overall coefficient
double Coeff = abs( -9.0*pA*a/c );

// Numerator terms for ND3 calculation
double uCHnd3 = abs( -pC*lC*LHe/aa );
double uCBnd3 = abs( pCD*lCD*LHe/aa );
double uFOnd3 = abs( lC*lCD*(pC-pCD)/aa );
// Denominator terms for ND3 calculation
double dCHnd3 = abs( 9.0*lC*pC*LHe*pD/cc );
double dCBnd3 = abs( -2.0*lCD*pCD*LHe*pD/cc );
double dFOnd3 = abs( -bb*lC*lCD/cc );
// Overall coefficient
double Coeffnd3 = abs( -9.0*pD*aa/cc );

// These variables are used to calculate the packing fraction (PF) for the NH3 target cell
double A = lC*lCH*( 9*pA*pC - 2*pCH*(pA+3*pC));
double B = (6.0)*lC*pC*lCH*pCH;
double C = 2*lCH*pCH*LHe*pA - 9*lC*pC*LHe*pA + A;
double D = -2*lCH*pCH*LHe*pA;
double E = 9*lC*pC*LHe*pA;
// These variables are used to calculate the packing fraction (PF) for the ND3 target cell
double AA = lC*lCD*( 9*pD*pC - 2*pCD*(pD+3*pC));
double BB = (6.0)*lC*pC*lCD*pCD;
double CC = 2*lCD*pCD*LHe*pD - 9*lC*pC*LHe*pD + AA;
double DD = -2*lCD*pCD*LHe*pD;
double EE = 9*lC*pC*LHe*pD;

// Calculates the NH3 packing fraction
double NH3PF(double nA, double nCH, double nC, double nET, double nF){
	double numerator = B*(nA - nET);
	double denominator = C*nET + D*nC + E*nCH - A*nF;
	if( denominator != 0.0 ) return numerator / denominator;
	else return 0.0; // error case
}
// Calculates the error in the NH3 PF
double NH3PFError(double nA, double nCH, double nC, double nET, double nF, double dnA, double dnCH, double dnC, double dnET, double dnF){
	double denominator = C*nET + D*nC + E*nCH - A*nF;
	if( denominator != 0.0 ){
		double DPA = B / denominator;
		double DPET = -DPA - (B*C*(nA - nET))/pow(denominator,2);
		double DPC = -D*B*(nA-nET)/pow(denominator,2);
		double DPCH = -E*B*(nA-nET)/pow(denominator,2);
		double DPF = A*B*(nA-nET)/pow(denominator,2);
		return sqrt( dnA*dnA*DPA*DPA + dnET*dnET*DPET*DPET + dnC*dnC*DPC*DPC + dnCH*dnCH*DPCH*DPCH + dnF*dnF*DPF*DPF );
	}
	else return 0.0; // error case

}
// Calculates the ND3 packing fraction
double ND3PF(double nA, double nCD, double nC, double nET, double nF){
	double numerator = BB*(nA - nET);
	double denominator = CC*nET + DD*nC + EE*nCD - AA*nF;
	if( denominator != 0.0 ) return numerator / denominator;
	else return 0.0; // error case
}
// Calculates the error in the ND3 PF
double ND3PFError(double nA, double nCD, double nC, double nET, double nF, double dnA, double dnCD, double dnC, double dnET, double dnF){
	double denominator = CC*nET + DD*nC + EE*nCD - AA*nF;
	if( denominator != 0.0 ){
		double DPA = BB / denominator;
		double DPET = -DPA - (BB*CC*(nA - nET))/pow(denominator,2);
		double DPC = -DD*BB*(nA-nET)/pow(denominator,2);
		double DPCH = -EE*BB*(nA-nET)/pow(denominator,2);
		double DPF = AA*BB*(nA-nET)/pow(denominator,2);
		return sqrt( dnA*dnA*DPA*DPA + dnET*dnET*DPET*DPET + dnC*dnC*DPC*DPC + dnCD*dnCD*DPCH*DPCH + dnF*dnF*DPF*DPF );
	}
	else return 0.0; // error case
}

// Calculates the dilution factor for NH3 targets using raw statistics only
double NH3DF(double nA, double nCH, double nC, double nET, double nF){
	double numerator = Coeff*(nA-nET)*(uCH*nCH - uCB*nC - uFO*nF + nET);
	double denominator = nA*(dCH*nCH - dCB*nC - dFO*nF - nET);
	if( denominator != 0.0 ) return numerator / denominator;
	else return 0; // error case
}

double NH3DFError(double nA, double nCH, double nC, double nET, double nF, double dnA, double dnCH, double dnC, double dnET, double dnF){
  double term1 = dnF*dnF*nA*nA*pow((nA-nET),2)*pow( ( (dCB*uFO-uCB*dFO)*nC+(uCH*dFO-dCH*uFO)*nCH+(uFO+dFO)*nET), 2 ); 
  double term2 = dnA*dnA*nET*nET*pow( (-uCB*nC+uCH*nCH+nET-uFO*nF),2)*pow( (dCB*nC-dCH*nCH+nET+dFO*nF),2); 
  double term3 = dnCH*dnCH*nA*nA*pow( (nA-nET),2)*pow( ( (uCH*dCB-uCB*dCH)*nC+(uCH+dCH)*nET+(uCH*dFO-dCH*uFO)*nF),2 );
  double term4 = dnC*dnC*nA*nA*pow( (nA-nET),2 )*pow( ( (uCH*dCB-uCB*dCH)*nCH+(uCB+dCB)*nET+(uCB*dFO-dCB*uFO)*nF ) ,2);
  double subt1 = uFO*dFO*nF*nF - ((dCH*uFO+uCH*dFO)*nCH+2*dFO*nET)*nF + nA*((uCB+dCB)*nC-(uCH+dCH)*nCH+(uFO+dFO)*nF);
  double subt2 = nC*( (dCB*uFO+uCB*dFO)*nF -((uCB*dCH+uCH*dCB)*nCH) -2*dCB*nET );
  double term5 = dnET*dnET*nA*nA*pow( (uCB*dCB*nC*nC+uCH*dCH*nCH*nCH+2*dCH*nCH*nET-nET*nET +subt1 +subt2 ), 2);
  
  double denom = nA*nA*pow( (dCB*nC -dCH*nCH +nET +dFO*nF), 2);
  double num = term1 + term2 + term3 + term4 +term5;
  
  if( denom != 0.0 && num > 0.0 ) return Coeff*sqrt(num)/denom;
  else return 0.0; // error case
  
}

// Calculates the dilution factor for ND3 targets using raw statistics only
double ND3DF(double nA, double nCD, double nC, double nET, double nF){
	//double numerator = Coeffnd3*(nA-nET)*(uCHnd3*nCD - uCBnd3*nC - uFOnd3*nF + nET);
	double numerator = Coeffnd3*(nA-nET)*(uCHnd3*nCD - uCBnd3*nC - uFOnd3*nF - nET);
	double denominator = nA*(dCHnd3*nCD - dCBnd3*nC + dFOnd3*nF - nET);
	if( denominator != 0.0 ) return numerator / denominator;
	else return 0; // error case
}

double ND3DFError(double nA, double nCD, double nC, double nET, double nF, double dnA, double dnCD, double dnC, double dnET, double dnF){
  double term1 = dnF*dnF*nA*nA*pow((nA-nET),2)*pow( ( (dCBnd3*uFOnd3-uCBnd3*dFOnd3)*nC+(uCHnd3*dFOnd3-dCHnd3*uFOnd3)*nCD+(uFOnd3+dFOnd3)*nET), 2 ); 
  double term2 = dnA*dnA*nET*nET*pow( (-uCBnd3*nC+uCHnd3*nCD+nET-uFOnd3*nF),2)*pow( (dCBnd3*nC-dCHnd3*nCD+nET+dFOnd3*nF),2); 
  double term3 = dnCD*dnCD*nA*nA*pow( (nA-nET),2)*pow( ( (uCHnd3*dCBnd3-uCBnd3*dCHnd3)*nC+(uCHnd3+dCHnd3)*nET+(uCHnd3*dFOnd3-dCHnd3*uFOnd3)*nF),2 );
  double term4 = dnC*dnC*nA*nA*pow( (nA-nET),2 )*pow( ( (uCHnd3*dCBnd3-uCBnd3*dCHnd3)*nCD+(uCBnd3+dCBnd3)*nET+(uCBnd3*dFOnd3-dCBnd3*uFOnd3)*nF ) ,2);
  double subt1 = uFOnd3*dFOnd3*nF*nF - ((dCHnd3*uFOnd3+uCHnd3*dFOnd3)*nCD+2*dFOnd3*nET)*nF + nA*((uCBnd3+dCBnd3)*nC-(uCHnd3+dCHnd3)*nCD+(uFOnd3+dFOnd3)*nF);
  double subt2 = nC*( (dCBnd3*uFOnd3+uCBnd3*dFOnd3)*nF -((uCBnd3*dCHnd3+uCHnd3*dCBnd3)*nCD) -2*dCBnd3*nET );
  double term5 = dnET*dnET*nA*nA*pow( (uCBnd3*dCBnd3*nC*nC+uCHnd3*dCHnd3*nCD*nCD+2*dCHnd3*nCD*nET-nET*nET +subt1 +subt2 ), 2);
  
  double denom = nA*nA*pow( (dCBnd3*nC -dCHnd3*nCD +nET +dFOnd3*nF), 2);
  double num = term1 + term2 + term3 + term4 +term5;
  
  if( denom != 0.0 && num > 0.0 ) return Coeffnd3*sqrt(num)/denom;
  else return 0.0; // error case
  
}

// These constants are used to make the calulation of the DF_NH3 using the PF a little easier
const double q = 3*pA / (2*pC*pCH); const double r = pC - pCH + LHe*(pCH/lC - pC/lCH);
const double s = -LHe*pCH/lC; const double t = LHe*pC/lCH; const double u = pCH - pC;

// Calculates the dilution factor using a previously calculated value of the packing fraction (PF)
double NH3DF_Thru_PF(double pf, double nA, double nCH, double nC, double nET, double nF){
  if( nA != 0.0 ) return (q*pf/nA)*(r*nET + s*nC + t*nCH + u*nF );
  else return 0.0; // error case
}

// Calculates the error in the dilution factor calculated with the PF as input; does NOT include error
// propagation from the PF right now
double NH3DF_Thru_PF_Error(double pf, double nA, double nCH, double nC, double nET, double nF, double dpf, double dnA, double dnCH, double dnC, double dnET, double dnF){
  // These terms collect items and makes things more collected
  double DFA = (-q*pf/(nA*nA))*( r*nET + s*nC + t*nCH + u*nF );
  double DFET = q*pf*r/nA;
  double DFNC = q*pf*s/nA;
  double DFCH = q*pf*t/nA;
  double DFNF = q*pf*u/nA;
  if( nA != 0.0 ){
    return sqrt( dnA*dnA*DFA*DFA + dnET*dnET*DFET*DFET + dnC*dnC*DFNC*DFNC + dnCH*dnCH*DFCH*DFCH + dnF*dnF*DFNF*DFNF );
  }
  else return 0.0;

}

// These functions calculate various pieces of the DF calculated through the PF for NH3
double NH3DF_Thru_PF_Numerator(double pf, double nA, double nCH, double nC, double nET, double nF){
  return (r*nET + s*nC + t*nCH + u*nF );
}
double NH3DF_Thru_PF_Numerator_Error(double pf, double nA, double nCH, double nC, double nET, double nF, double dpf, double dnA, double dnCH, double dnC, double dnET, double dnF){
  return sqrt( dnET*dnET*r*r + dnC*dnC*s*s + dnCH*dnCH*t*t + dnF*dnF*u*u );
}
// Plots the ammonia target counts for NH3 targets
double NH3DF_Thru_PF_NA(double pf, double nA, double nCH, double nC, double nET, double nF){
  if( nA != 0 )  return q/nA;
  else return 0;
}
double NH3DF_Thru_PF_NA_Error(double pf, double nA, double nCH, double nC, double nET, double nF, double dpf, double dnA, double dnCH, double dnC, double dnET, double dnF){
  if( nA != 0 )  return dnA*q/(nA*nA);
  else return 0;
}

// These constants are used to make the calulation of the DF_ND3 using the PF a little easier
const double qq = 3*pD / (2*pC*pCD); const double rr = pC - pCD + LHe*(pCD/lC - pC/lCD);
const double ss = -LHe*pCD/lC; const double tt = LHe*pC/lCD; const double uu = pCD - pC;

// These functions calculate various pieces of the DF calculated through the PF for ND3
double ND3DF_Thru_PF_Numerator(double pf, double nA, double nCD, double nC, double nET, double nF){
  return (rr*nET + ss*nC + tt*nCD + uu*nF );
}
double ND3DF_Thru_PF_Numerator_Error(double pf, double nA, double nCD, double nC, double nET, double nF, double dpf, double dnA, double dnCD, double dnC, double dnET, double dnF){
  return sqrt( dnET*dnET*rr*rr + dnC*dnC*ss*ss + dnCD*dnCD*tt*tt + dnF*dnF*uu*uu );
}
// Plots the ammonia target counts for ND3 targets
double ND3DF_Thru_PF_NA(double pf, double nA, double nCD, double nC, double nET, double nF){
  if( nA != 0 )  return qq/nA;
  else return 0;
}
double ND3DF_Thru_PF_NA_Error(double pf, double nA, double nCD, double nC, double nET, double nF, double dpf, double dnA, double dnCD, double dnC, double dnET, double dnF){
  if( nA != 0 )  return dnA*qq/(nA*nA);
  else return 0;
}


// Calculates the dilution factor using a previously calculated value of the packing fraction (PF)
double ND3DF_Thru_PF(double pf, double nA, double nCD, double nC, double nET, double nF){
  if( nA != 0.0 ) return (qq*pf/nA)*(rr*nET + ss*nC + tt*nCD + uu*nF );
  else return 0.0; // error case
}
// Calculates the error in the dilution factor calculated with the PF as input; does NOT include error
// propagation from the PF right now
double ND3DF_Thru_PF_Error(double pf, double nD, double nCD, double nC, double nET, double nF, double dpf, double dnD, double dnCD, double dnC, double dnET, double dnF){
  // These terms collect items and makes things more collected
  double DFD = (-qq*pf/(nD*nD))*( rr*nET + ss*nC + tt*nCD + uu*nF );
  double DFET = qq*pf*rr/nD;
  double DFNC = qq*pf*ss/nD;
  double DFCD = qq*pf*tt/nD;
  double DFNF = qq*pf*uu/nD;
  if( nD != 0.0 ){
    return sqrt( dnD*dnD*DFD*DFD + dnET*dnET*DFET*DFET + dnC*dnC*DFNC*DFNC + dnCD*dnCD*DFCD*DFCD + dnF*dnF*DFNF*DFNF );
  }
  else return 0.0;

}

//Run_Number  Target  Beam_Cur_Req  Beam_Energy  Num_Events  HWP_Status  Target_Pol  Solenoid_Scale  Torus_Scale
// This class holds information about each of the runs from RCDB
class Run{
    private:
	int RunNumber = 0; 
	string TargetType;
	string BeamCurrent = "N/A"; // in nC
	double BeamEnergy = 0; // in GeV
        int NumberEvents = 0;
	int HWPStatus = -1; // 1=in, 0=out, -1=uninitialized
	double TargetPolarization = 0; // NMR Tpol for NH3 and ND3; zero for other targets
	double OfflineTPol = 0; // From Ishara's analysis
	double OfflineTPolErr = 0;
	double SolenoidScale = 0; // Current direction in solenoid
	double TorusScale = 0; // Current direction in torus
	string Epoch; // Run range 'epoch' that this run is a part of; used for DF/PF calculation
	int EpochNum; // Number of the epoch saved as an integer
    public:
	// Constructor
	Run() = default;
	// Mutators
	void SetRunInfo(int rn, string tt, string bc, double be, int ne, int hwp, double tp, double ss, double ts, string epoch){
	    RunNumber = rn; TargetType = tt; BeamCurrent = bc; BeamEnergy = be; NumberEvents = ne;
	    HWPStatus = hwp; TargetPolarization = tp; SolenoidScale = ss; TorusScale = ts; Epoch = epoch;
	    int length = epoch.length();
	    string epochnum; epochnum += epoch[length-2]; epochnum += epoch[length-1];
	    //cout << epoch <<" "<< epochnum <<" "<< rn <<" ";
	    EpochNum = stoi(epochnum); //cout << EpochNum << endl;
	}
	void SetOffileTPol(double tpol, double tpolerr){
	    OfflineTPol = tpol; OfflineTPolErr = tpolerr;
	}
	// Accessors
	int getRunNumber() const{ return RunNumber;}
	string getTargetType() const{ return TargetType;}
	string getBeamCurrent() const{ return BeamCurrent;}
	double getBeamEnergy() const{ return BeamEnergy;}
	int getNumberEvents() const{ return NumberEvents;}
	int getHWPStatus() const{ return HWPStatus;}
	double getTargetPolarization() const{ return TargetPolarization;}
	double getOfflineTPol() const{ return OfflineTPol;}
	double getOfflineTPolErr() const{ return OfflineTPolErr;}
	double getSolenoidScale() const{ return SolenoidScale;}
	double getTorusScale() const{ return TorusScale;}
	string getEpoch() const{ return Epoch;}
	int getEpochNum() const{ return EpochNum;}
	void Print(){
	  cout <<RunNumber<<" "<<TargetType<<" "<<BeamCurrent<<" "<<BeamEnergy<<" "<<NumberEvents<<" ";
	  cout <<HWPStatus<<" "<<TargetPolarization<<" "<<SolenoidScale<<" "<<TorusScale<<" "<<Epoch<<endl;
	}
	// Destructor
	~Run() =  default;
};

// This class holds the run for a given run period
class RunPeriod{
    private:
	vector<Run> Runs;
    public:
	// Constructor
	RunPeriod(){ SetRunPeriod(); }
	// Mutator
	void SetRunPeriod(){
	    string fPath = THIS_DIR + "RGC_Run_Info.txt";
	    ifstream fin( fPath.c_str()  );
	    if( !fin.fail() ){
	        string line;
	        getline(fin,line); // throw away header row
	        while(getline(fin,line)){
	            stringstream sin(line);
		    int rn; string tt; string bc; double be; int ne; int hwp; double tp; double ss; double ts; string epoch;
		    sin >> rn >> tt >> bc >> be >> ne >> hwp >> tp >> ss >> ts >> epoch;
		    Run thisRun;
		    thisRun.SetRunInfo(rn, tt, bc, be, ne, hwp, tp, ss, ts, epoch);
		    Runs.push_back(thisRun);
		}
	    }
	    else cout <<"Couldn't find input file; run period not set.\n";
	    fin.close();
	    string fPath2 = THIS_DIR + "Offline_NMR_Values.txt";
	    ifstream fin2( fPath2.c_str() );
	    if( !fin2.fail() ){
		string line; getline( fin2, line ); // throw away header row
		while(getline(fin2,line)){
		    stringstream sin(line);
		    int run, cell; string species; double avgOnline, avgOffline, avgOfflineErr, runDone;
		    sin >> run >> species >> cell >> avgOnline >> avgOffline >> avgOfflineErr >> runDone;
		    // Now slot into the appropriate run
		    for(size_t i=0; i<Runs.size(); i++){
			if(Runs[i].getRunNumber() == run){
			    Runs[i].SetOffileTPol( avgOffline, avgOfflineErr );
			    i = Runs.size();
			}
		    }
		}
	    }
	    fin2.close();
	}
	// Accessors
	double getTargetPolarization(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getTargetPolarization();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}
	double getOfflineTPol(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getOfflineTPol();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}
	double getOfflineTPolErr(int runnum) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == runnum ){
		    if( Runs[i].getTargetType() == "NH3" || Runs[i].getTargetType() == "ND3" )
		       	return Runs[i].getOfflineTPolErr();
		    else return 0;
		}
	    }
	    cout <<"Couldn't find run "<< runnum <<". Check inputs.\n";
	    return 0;
	}

	int getEpochNum(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run ) return Runs[i].getEpochNum();
	    }
	    return -1;
	}
	// If the target is empty or foils, then ignore the "misc" error bit. Else, include the bit
	bool applyMiscBit(int run) const{
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getRunNumber() == run && Runs[i].getTargetType() == "Empty" ) return false;
	    }
	    return true;
	}
	vector<int> getRunEpoch(string epoch) const{
	    vector<int> RunsInEpoch;
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getEpoch() == epoch ) RunsInEpoch.push_back( Runs[i].getRunNumber() );
	    }
	    if( RunsInEpoch.size() == 0 ) cout << "Couldn't find runs for epoch "<<epoch<<endl;
	    return RunsInEpoch;
	}
	vector<int> getSubEpoch(string epoch, int targetPol) const{
	    // By selecting the sign of the target polarization, you can select the subset of runs from the
	    // epoch that correspond to the sign of the target polarization
	    vector<int> RunsInSubEpoch;
	    for( auto run : Runs ){
		if( run.getEpoch() == epoch && targetPol < 0 ) RunsInSubEpoch.push_back( run.getRunNumber() );
		else if( run.getEpoch() == epoch && targetPol > 0) RunsInSubEpoch.push_back( run.getRunNumber() );
	    }
	    if( RunsInSubEpoch.size() == 0 ) cout << "Couldn't find runs for epoch "<< epoch <<" with P_target ~ "<< targetPol << endl;
	    return RunsInSubEpoch;
	}
	vector<int> getDataSet(string Dataset, string TorPol, string SolPol) const{
	    // Returns the subset of runs consisting of this dataset, torus polarity, solenoid polarity
	    vector<int> RunsInDataSet;
	    for( auto run : Runs ){
		if( Dataset == "Su22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16137 && run.getRunNumber() <= 16772 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Fa22" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Pos" &&
		    run.getSolenoidScale() > 0 && run.getRunNumber() >= 16843 && run.getRunNumber() <= 17408 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Neg" && run.getTorusScale() < 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17811 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
		else if( Dataset == "Sp23" && TorPol == "Pos" && run.getTorusScale() > 0 && SolPol == "Neg" &&
		    run.getSolenoidScale() < 0 && run.getRunNumber() >= 17482 && run.getRunNumber() <= 17811 ){
		    RunsInDataSet.push_back( run.getRunNumber() );
		}
	    }
	    return RunsInDataSet;
	}
	vector<int> getAllGoodRuns() const{
	    vector<int> AllGoodRuns;
	    for(size_t i=0; i<Runs.size(); i++){
		AllGoodRuns.push_back( Runs[i].getRunNumber() );
	    }
	    if( AllGoodRuns.size() == 0 ) cout << "Couldn't find any good runs\n";
	    return AllGoodRuns;
    
	}
	void Print(){
	    for(size_t i=0; i<Runs.size(); i++){
		Runs[i].Print();
	    }
	}
	// Destructor
	~RunPeriod() = default;
};

// This vector contains the bin boundaries for the Q2 bins
const vector<double> Q2_Bin_Bounds = { /*0.9188, 1.0969,*/ 1.3094, 1.5632, 1.8661, 2.2277,
	2.6594, 3.1747, 3.7899, 4.5243, 5.4009, 6.4475, 7.6969, 9.1884, 10.9689};
// This vector contains the bin boundaries for the x bins
//const vector<double> X_Bin_Bounds = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
//	0.5, 0.55, 0.6, 0.65, 0.7, 0.75};
//const vector<double> X_Bin_Bounds = {   0.0875, 0.1125, 0.1375, 0.1625, 0.1875, 0.2125, 0.2375, 0.2625,
//	0.2875, 0.3125, 0.3375, 0.3625, 0.3875, 0.4125, 0.4375, 0.4625, 0.4875, 0.5125, 0.5375, 0.5625, 0.5875, 
//	0.6125, 0.6375, 0.6625, 0.6875, 0.7125, 0.7375, 0.7625, 0.7875 };

const vector<double> X_Bin_Bounds = {0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
	0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,
	0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675,
	0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875 };

// This struct holds information about the counts of different target types
class Count{

    private:
	string TargetType = "N/A";
	double Nm = 0.0; double Np = 0.0;
	double FCm = 0.0; double FCp = 0.0;
	double NormNm = 0.0; double NormNp = 0.0;
    public:
	// Constructor
	Count() = default;
	Count(string targtype) : TargetType(targtype){}
	// Mutators
	void SetTarget(string targtype){
		TargetType = targtype;
	}
	void AddCounts(double nm, double np){
		Nm += nm; Np += np;
	}
	void AddFCCharge(double fcm, double fcp){
		FCm += fcm; FCp += fcp;
	}
	void SetNormCounts(){
	    if( FCm != 0.0 && FCp != 0.0 ){
		NormNm = Nm/FCm; NormNp = Np/FCp;
	    }
	    // Don't print anything if it's an empty bin
	    else if(Nm > 0.0 && Np > 0.0){
		cout << "ERROR: Invalid FC charge. Unable to set FC values.\n";
		cout << "FCm = "<< FCm <<", FCp = "<< FCp << endl;
	    }
	}
	// Accessors
	double getNm() const{return Nm;}
	double getNp() const{return Np;}
	double getNt() const{return Nm+Np;}
	double getFCm() const{return FCm;}
	double getFCp() const{return FCp;}
	double getFCt() const{return FCm+FCp;}
	double getNormNm() const{
	    if( FCm > 0.0 ) return Nm / FCm;
	    else return 0.0;
	}
	double getNormNp() const{
	    if( FCp > 0.0 ) return Np / FCp;
	    else return 0.0;
	}
	//double getNormNt() const{return NormNm+NormNp;}
	double getNormNt() const{
	    if( FCm > 0.0 && FCp > 0.0 )
		return (Nm+Np)/(FCm+FCp);
	    else
		return 0.0;
	}

	void Print() const{
		cout << "Target = "<< TargetType;
		cout << ",	Nm = "<< Nm <<",	Np = "<< Np;
	       	cout << ",	FCm = "<< FCm <<",	FCp = "<< FCp;
		cout << ",	NormNm = "<<NormNm<<",	NormNp = "<<NormNp << endl;
	}
	// Destructor
	~Count() = default;
};

// This struct holds the information for each of the x, Q2 bins
class Bin{

    private:
	// Minimum, maximum, and middle value for the bin (mid = (max+min)/2)
	double Q2_Min = 0.0; double Q2_Max = 0.0; double Q2_Mid = 0.0;
	double X_Min  = 0.0; double X_Max  = 0.0; double X_Mid  = 0.0;
	// Holds the calculated dilution factor and its associated error for this bin
	double DF_NH3 = 0.0; double ErrDF_NH3 = 0.0; // This calculates the DF using unique PF_bath for each bin
	double DF_FixedPF_NH3 = 0.0;   // This DF value is calculated using a fixed value for the PF_bath across all bins
	double ErrDF_FixedPF_NH3 = 0.0;
	double DF_ND3 = 0.0; double ErrDF_ND3 = 0.0;
	double DF_FixedPF_ND3 = 0.0;   // This DF value is calculated using a fixed value for the PF_bath across all bins
	double ErrDF_FixedPF_ND3 = 0.0;
	// Holds the calculated packing fraction and the associated error for this bin
	double PF_bath_NH3 = 0.0; double ErrPF_bath_NH3 = 0.0;
	double PF_cell_NH3 = 0.0; double ErrPF_cell_NH3 = 0.0;
	double PF_bath_ND3 = 0.0; double ErrPF_bath_ND3 = 0.0;
	double PF_cell_ND3 = 0.0; double ErrPF_cell_ND3 = 0.0;
	// Holds information for various pieces of the DF calculation; for debugging
	double Numerator_NH3_DF_Thru_PF = 0.0; double Numerator_NH3_DF_Thru_PF_Error = 0.0;
	double NA_Counts_NH3_DF_Thru_PF = 0.0; double NA_Counts_NH3_DF_Thru_PF_Error = 0.0;
	double Numerator_ND3_DF_Thru_PF = 0.0; double Numerator_ND3_DF_Thru_PF_Error = 0.0;
	double NA_Counts_ND3_DF_Thru_PF = 0.0; double NA_Counts_ND3_DF_Thru_PF_Error = 0.0;
	// Holds the non-FC corrected raw double-spin asymmetry for NH3 and ND3
	double AllRawNoFCNH3 = 0.0; double AllRawNoFCND3 = 0.0;
	double AllRawNoFCNH3Err = 0.0; double AllRawNoFCND3Err = 0.0;
	// Holds the raw double-spin asymmetry for ND3
	double AllRawNH3 = 0.0; double AllRawND3 = 0.0;
	double AllRawNH3Err = 0.0; double AllRawND3Err = 0.0;
	// Holds the physical double-spin asymmetry
	double AllPhysNH3 = 0.0; double AllPhysND3 = 0.0;
	double AllPhysNH3Err = 0.0; double AllPhysND3Err = 0.0;
	// Holds the theoretical value of the double-spin asymmetry for NH3 in this bin
	double AllTheoryNH3 = 0.0; double AllTheoryND3 = 0.0;
        // Holds various other kinematic factors (A2 is from Sebastian's theoretical values)
        double DepolFactorNH3= 0.0; double EtaNH3 = 0.0; double A2NH3 = 0.0;
        double DepolFactorND3= 0.0; double EtaND3 = 0.0; double A2ND3 = 0.0;
	// Holds the value of A1 from data and theory
	double A1_NH3 = 0.0; double ErrA1_NH3 = 0.0;
	double A1_ND3 = 0.0; double ErrA1_ND3 = 0.0;
	double A1_Theory_NH3 = 0.0; double A1_Theory_ND3 = 0.0;

	// Hold the count types for every kind of target
	Count NH3_Counts = Count("NH3"); // Ammonia
	Count ND3_Counts = Count("ND3"); // Deuterated ammonia
	Count C_Counts   = Count("C");   // Carbon foils
	Count CH2_Counts = Count("CH2"); // CH2
	Count CD2_Counts = Count("CD2"); // CD2
	Count ET_Counts  = Count("ET");  // Empty target, with LHe
	Count F_Counts   = Count("F");   // Empty target, no LHe (Aluminum foils only)
	Count Input_Loop = Count("Input");//USED ONLY IN "Get_DF_By_Sector.C" FOR WRITING TO INPUT TEXT FILES
    public:
	// Constructors
	Bin() = default;

	// Mutators
	void SetQ2Bins(double qmin, double qmax){
		Q2_Min = qmin; Q2_Max = qmax; Q2_Mid = (qmin+qmax)/2.0;
	}
	void SetXBins(double xmin, double xmax){
		X_Min = xmin; X_Max = xmax; X_Mid = (xmin+xmax)/2.0;
	}
	void SetDF_NH3(double df, double dferr){
		DF_NH3 = df; ErrDF_NH3 = dferr;
	}
	void SetDF_ND3(double df, double dferr){
		DF_ND3 = df; ErrDF_ND3 = dferr;
	}
	void SetDF_FixedPF_NH3(double df, double dferr){
		DF_FixedPF_NH3 = df; ErrDF_FixedPF_NH3 = dferr;
	}
	void SetDF_FixedPF_ND3(double df, double dferr){
		DF_FixedPF_ND3 = df; ErrDF_FixedPF_ND3 = dferr;
	}
	void SetAllRaw(double allraw, string target){
		if( target == "NH3" ) AllRawNH3 = allraw;
		else if( target == "ND3" ) AllRawND3 = allraw;
	}
	void SetAllRawErr(double allrawerr, string target){
		if( target == "NH3" ) AllRawNH3Err = allrawerr;
		else if( target == "ND3" ) AllRawND3Err = allrawerr;
	}
	void SetAllRawNoFC(double allraw, string target){
		if( target == "NH3" ) AllRawNoFCNH3 = allraw;
		else if( target == "ND3" ) AllRawNoFCND3 = allraw;
	}
	void SetAllRawNoFCErr(double allrawerr, string target){
		if( target == "NH3" ) AllRawNoFCNH3Err = allrawerr;
		else if( target == "ND3" ) AllRawNoFCND3Err = allrawerr;
	}
	void SetAllTheory(double allth, string target){
		if( target == "NH3" ) AllTheoryNH3 = allth;
		else if( target == "ND3" ) AllTheoryND3 = allth;
	}
	void SetAllPhys(double allphys, double errallphys, string target){
		if( target == "NH3" ){ AllPhysNH3 = allphys; AllPhysNH3Err = errallphys;}
		else if( target == "ND3" ){ AllPhysND3 = allphys; AllPhysND3Err = errallphys;}
	}
	void SetKinematicFactors(double depol, double eta, double a1, double a2, string target){
		if( target == "NH3" ){ 
		    DepolFactorNH3 = depol; EtaNH3 = eta; A1_Theory_NH3 = a1; A2NH3 = a2;
		}
		else if( target == "ND3" ){ 
		    DepolFactorND3 = depol; EtaND3 = eta; A1_Theory_ND3 = a1; A2ND3 = a2;
		}
	}
	void SetA1(double a1, double err_a1, string target ){
		if( target == "NH3" ){
		    A1_NH3 = a1; ErrA1_NH3 = err_a1;
		}
		else if( target == "ND3" ){
		    A1_ND3 = a1; ErrA1_ND3 = err_a1;
		}
	}
	void CalculateA1(double PBPT, double pbpt_err){
		//cout << "For bin "<< Q2_Mid <<" "<< X_Mid <<": DepolNH3 = "<< DepolFactorNH3 <<", EtaNH3 = "<< EtaNH3 <<", A2_NH3 = "<< A2NH3;
	        //cout << ", DF_FixedPF_NH3 = "<< DF_FixedPF_NH3 <<", PbPt = "<< PBPT << endl;
		double pbpt = abs(PBPT);
		if( DepolFactorNH3 > 0 && EtaNH3 > 0 && A2NH3 > 0 && pbpt > 0 && DF_FixedPF_NH3 > 0 ){
		    AllPhysNH3 = (AllRawNH3 / (pbpt*DF_FixedPF_NH3*DepolFactorNH3)) - (EtaNH3*A2NH3);
		    AllPhysNH3Err = (AllRawNH3 / (pbpt*DF_FixedPF_NH3*DepolFactorNH3))*sqrt( pbpt_err*pbpt_err/(pbpt*pbpt) + ErrDF_FixedPF_NH3*ErrDF_FixedPF_NH3/(DF_FixedPF_NH3*DF_FixedPF_NH3) );
		}
		else{
		    cout << "ERROR: Invalid values for bin "<< Q2_Mid <<" "<< X_Mid << endl;
		}
	}
	void CalculateDF(){
		NH3_Counts.SetNormCounts(); double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		ND3_Counts.SetNormCounts(); double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		CH2_Counts.SetNormCounts(); double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		CD2_Counts.SetNormCounts(); double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		C_Counts.SetNormCounts();   double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		ET_Counts.SetNormCounts();  double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		F_Counts.SetNormCounts();   double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();
		// Calculate NH3 dilution factor for this bin
		if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 ){
			double dfnh3 = NH3DF( nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			double errdfnh3 = NH3DFError( nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_NH3 = dfnh3; ErrDF_NH3 = errdfnh3;
		}
		else if(nA > 0.0){
			cout <<"ERROR in CalculateDF(): Some FC charges are zero for NH3 calculation; check inputs.\n";
		}
		// Calculate ND3 dilution factor for this bin
		// Using CH2 in place of CD2 for the time being
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 ){
			double dfnd3 = ND3DF( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errdfnd3 = ND3DFError( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_ND3 = dfnd3; ErrDF_ND3 = errdfnd3;
		}
		else if(nD > 0.0){
			cout <<"ERROR in CalculateDF(): Some FC charges are zero for ND3 calculation; check inputs.\n";
		}
		
	}
	void CalculatePF(){
		NH3_Counts.SetNormCounts(); double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		ND3_Counts.SetNormCounts(); double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		CH2_Counts.SetNormCounts(); double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		CD2_Counts.SetNormCounts(); double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		C_Counts.SetNormCounts();   double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		ET_Counts.SetNormCounts();  double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		F_Counts.SetNormCounts();   double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();
		// Calculate NH3 packing fraction for this bin
		if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 ){
			double pfnh3 = NH3PF( nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			double errpfnh3 = NH3PFError( nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			PF_bath_NH3 = pfnh3; ErrPF_bath_NH3 = errpfnh3;
			PF_cell_NH3 = (LHe / 5.0)*pfnh3; ErrPF_cell_NH3 = (LHe / 5.0)*errpfnh3;
		}
		else if(nA > 0.0){
			cout <<"***********************************************************************************\n";
			cout <<"ERROR in CalculatePF(): Some FC charges are zero for NH3 calculation; check inputs.\n";
			cout <<"For Bin Q^2 = "<< Q2_Mid <<" (GeV^2), X = "<< X_Mid << endl;
			cout <<"  -> fA = " << fA << endl;
			cout <<"  -> fCH = " << fCH << endl;
			cout <<"  -> fC = " << fC << endl;
			cout <<"  -> fET = " << fET << endl;
			cout <<"  -> fF = " << fF << endl;
			cout <<"***********************************************************************************\n";
		}
		// Calculate ND3 packing fraction for this bin
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 ){
			double pfnd3 = ND3PF( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errpfnd3 = ND3PFError( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			PF_bath_ND3 = pfnd3; ErrPF_bath_ND3 = errpfnd3;
			PF_cell_ND3 = (LHe / 5.0)*pfnd3; ErrPF_cell_ND3 = (LHe / 5.0)*errpfnd3;
		}
		else if(nD > 0.0){
			cout <<"***********************************************************************************\n";
			cout <<"ERROR in CalculatePF(): Some FC charges are zero for ND3 calculation; check inputs.\n";
			cout <<"For Bin Q^2 = "<< Q2_Mid <<" (GeV^2), X = "<< X_Mid << endl;
			cout <<"  -> fD = " << fD << endl;
			cout <<"  -> fCD = " << fCD << endl;
			cout <<"  -> fC = " << fC << endl;
			cout <<"  -> fET = " << fET << endl;
			cout <<"  -> fF = " << fF << endl;
			cout <<"***********************************************************************************\n";
		}
	}
	// This function is called in a member function of the DataSet object, using the same name. Don't call this
	// object on its own; it should only be used by the member function in the DataSet class!!!
	void CalculateDFThruPF( double thisPF, string targetType ){
		NH3_Counts.SetNormCounts(); double nA = NH3_Counts.getNt(); double fA = NH3_Counts.getFCt();
		ND3_Counts.SetNormCounts(); double nD = ND3_Counts.getNt(); double fD = ND3_Counts.getFCt();
		CH2_Counts.SetNormCounts(); double nCH= CH2_Counts.getNt(); double fCH= CH2_Counts.getFCt();
		CD2_Counts.SetNormCounts(); double nCD= CD2_Counts.getNt(); double fCD= CD2_Counts.getFCt();
		C_Counts.SetNormCounts();   double nC = C_Counts.getNt();   double fC = C_Counts.getFCt();
		ET_Counts.SetNormCounts();  double nET = ET_Counts.getNt(); double fET = ET_Counts.getFCt();
		F_Counts.SetNormCounts();   double nF = F_Counts.getNt();   double fF = F_Counts.getFCt();        
		// Calculate NH3 packing fraction for this bin
		if( targetType == "NH3" ){
		    if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 && thisPF>0 ){
			double dfnh3 = NH3DF_Thru_PF( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			double errdfnh3 = NH3DF_Thru_PF_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_FixedPF_NH3 = dfnh3; ErrDF_FixedPF_NH3 = errdfnh3;
			Numerator_NH3_DF_Thru_PF = NH3DF_Thru_PF_Numerator( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			Numerator_NH3_DF_Thru_PF_Error = NH3DF_Thru_PF_Numerator_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			NA_Counts_NH3_DF_Thru_PF = NH3DF_Thru_PF_NA( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			NA_Counts_NH3_DF_Thru_PF_Error = NH3DF_Thru_PF_NA_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
		    }
		    else{
			cout <<"ERROR in CalculateDFThruPF(): Some FC charges are zero for NH3 calculation; check inputs.\n";
		    }
		}
		// Calculate ND3 packing fraction for this bin
		else if( targetType == "ND3" ){
		    if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 && thisPF>0 ){
			double dfnd3 = ND3DF_Thru_PF( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errdfnd3 = ND3DF_Thru_PF_Error(thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_FixedPF_ND3 = dfnd3; ErrDF_FixedPF_ND3 = errdfnd3;
			Numerator_ND3_DF_Thru_PF = ND3DF_Thru_PF_Numerator( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			Numerator_ND3_DF_Thru_PF_Error = ND3DF_Thru_PF_Numerator_Error(thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			NA_Counts_ND3_DF_Thru_PF = ND3DF_Thru_PF_NA( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			NA_Counts_ND3_DF_Thru_PF_Error = ND3DF_Thru_PF_NA_Error(thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			//cout << dfnd3 <<" +- "<< errdfnd3 << endl;
		    }
		    else{
			cout <<"ERROR in CalculateDFThruPF(): Some FC charges are zero for ND3 calculation; check inputs.\n";
		    }
		}
		else cout <<"ERROR in CalculateDFThruPF(): Invalid target type. Can only calculate DF(PF) for NH3 and ND3 targets.\n";
	}
	// Calculates the raw double-spin asymmetry for this bin
	void CalculateAllRaw(){
	    NH3_Counts.SetNormCounts();
	    double normNm = NH3_Counts.getNormNm(); double normNp = NH3_Counts.getNormNp();
	    double nm = NH3_Counts.getNm(); double np = NH3_Counts.getNp();
	    //double nm = NH3_Counts.getNormNm(); double np = NH3_Counts.getNormNp();
	    //double fcm= NH3_Counts.getFCm(); double fcp = NH3_Counts.getFCp();
	    if( normNm > 0.0 && normNp > 0.0 ){ 
		AllRawNH3 = (normNm - normNp) / (normNm + normNp);
		//AllRawNH3Err = (2.0/pow((nm + np),2))*sqrt( (np*np*nm)/(fcm*fcm) + (nm*nm*np)/(fcp*fcp));
		AllRawNH3Err = 0.5 * sqrt( (nm + np)/(nm*np) );
		AllRawNoFCNH3 = (nm - np)/(nm + np);
		AllRawNoFCNH3Err = 0.5 * sqrt( (nm + np)/(nm*np) );
	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawNH3 = 0\n";
		AllRawNH3 = 0.0;
		AllRawNH3Err = 0.0;
	    }
	    double normND3Nm = ND3_Counts.getNormNm(); double normND3Np = ND3_Counts.getNormNp();
	    double nmND3 = ND3_Counts.getNm(); double npND3 = ND3_Counts.getNp();
	    if( normND3Nm > 0.0 && normND3Np > 0.0 ){
		AllRawND3 = (normND3Nm - normND3Np) / (normND3Nm + normND3Np);
		AllRawND3Err = (2.0/pow((nmND3 + npND3),2))*sqrt(npND3*npND3*nmND3 + nmND3*nmND3*npND3);
		AllRawNoFCND3 = (nmND3 - npND3)/(nmND3 + npND3);
		AllRawNoFCND3Err = (2.0/pow((nmND3 + npND3),2))*sqrt(npND3*npND3*nmND3 + nmND3*nmND3*npND3);
	    }
	    else{
		//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawND3 = 0\n";
		AllRawND3 = 0.0;
		AllRawND3Err = 0.0;
	    }
	}
	// Calculates the physical double-spin asymmetry for this bin
	void CalculateAllPhys( double PbPt, double PbPtErr ){
    	    CalculateAllRaw(); // Make sure the raw asymmetry is set first
	    if( AllRawNH3 != 0 && DF_FixedPF_NH3 != 0){
		AllPhysNH3 = AllRawNH3 / (PbPt*DF_FixedPF_NH3);
		AllPhysNH3Err = (1.0/(PbPt*DF_FixedPF_NH3))*sqrt( pow(AllRawNH3Err,2) + pow(AllRawNH3*ErrDF_FixedPF_NH3,2)/pow(DF_FixedPF_NH3,2) + pow(AllRawNH3*PbPtErr,2)/(PbPt*PbPt)  );
	    }
	    if( AllRawND3 != 0 ){
		AllPhysND3 = AllRawND3 / (PbPt*DF_FixedPF_ND3);
	    }
	}

	void AddCounts(double nm, double np, string targtype){
		if(targtype == "NH3") NH3_Counts.AddCounts(nm, np);
		else if(targtype == "ND3") ND3_Counts.AddCounts(nm, np);
		else if(targtype == "CH2") CH2_Counts.AddCounts(nm, np);
		else if(targtype == "CD2") CD2_Counts.AddCounts(nm, np);
		else if(targtype == "C") C_Counts.AddCounts(nm, np);
		else if(targtype == "ET") ET_Counts.AddCounts(nm, np);
		else if(targtype == "F") F_Counts.AddCounts(nm, np);
		else if(targtype == "Input") Input_Loop.AddCounts(nm, np);
		else{
		    cout <<"ERROR: Couldn't find targtype type. Added zero counts\n";
		}
	}
	void AddFCCharge(double fcm, double fcp, string targtype){
		if(targtype == "NH3") NH3_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "ND3") ND3_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "CH2") CH2_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "CD2") CD2_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "C") C_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "ET") ET_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "F") F_Counts.AddFCCharge(fcm, fcp);
		else if(targtype == "Input") Input_Loop.AddFCCharge(fcm, fcp);
		else{
		    cout <<"ERROR: Couldn't find targtype type. Added zero counts\n";
		}
	}
	void NormalizeAllCounts(){
		NH3_Counts.SetNormCounts(); ND3_Counts.SetNormCounts();
		CH2_Counts.SetNormCounts(); CD2_Counts.SetNormCounts();
		C_Counts.SetNormCounts(); ET_Counts.SetNormCounts();
		F_Counts.SetNormCounts(); Input_Loop.SetNormCounts();
	}
	// Accessors
	double getBinQ2() const{ return Q2_Mid; }
	double getBinX() const{ return X_Mid; }
	double getDF_NH3() const{ return DF_NH3; }
	double getErrDF_NH3() const{ return ErrDF_NH3; }
	double getDF_FixedPF_NH3() const{ return DF_FixedPF_NH3; }
	double getErrDF_FixedPF_NH3() const{ return ErrDF_FixedPF_NH3; }
	double getDF_FixedPF_ND3() const{ return DF_FixedPF_ND3; }
	double getErrDF_FixedPF_ND3() const{ return ErrDF_FixedPF_ND3; }
	double getDF_ND3() const{ return DF_ND3; }
	double getErrDF_ND3() const{ return ErrDF_ND3; }
	double getPF_bath_NH3() const{ return PF_bath_NH3; }
	double getErrPF_bath_NH3() const{ return ErrPF_bath_NH3; }
	double getPF_cell_NH3() const{ return PF_cell_NH3; }
	double getErrPF_cell_NH3() const{ return ErrPF_cell_NH3; }
	double getPF_bath_ND3() const{ return PF_bath_ND3; }
	double getErrPF_bath_ND3() const{ return ErrPF_bath_ND3; }
	double getPF_cell_ND3() const{ return PF_cell_ND3; }
	double getErrPF_cell_ND3() const{ return ErrPF_cell_ND3; }
	double getDepolFactorNH3() const{ return DepolFactorNH3; }
	double getEtaNH3() const{ return EtaNH3; }
	double getA2NH3() const{ return A2NH3; }
	double getDepolFactorND3() const{ return DepolFactorND3; }
	double getEtaND3() const{ return EtaND3; }
	double getA2ND3() const{ return A2ND3; }
	double getA1NH3() const{ return A1_NH3;}
	double getErrA1NH3() const{ return ErrA1_NH3;}
	double getA1ND3() const{ return A1_ND3;}
	double getErrA1ND3() const{ return ErrA1_ND3;}
	double getAllPhysNH3() const{ return AllPhysNH3; }
	double getAllPhysNH3Err() const{ return AllPhysNH3Err; }
	
	double getNumeratorDFThruPF(string target) const{
	    if( target == "NH3" ) return Numerator_NH3_DF_Thru_PF;
	    else if( target == "ND3" ) return Numerator_ND3_DF_Thru_PF;
	    else return 0.0;
	}
	double getNumeratorDFThruPFError(string target) const{
	    if( target == "NH3" ) return Numerator_NH3_DF_Thru_PF_Error;
	    else if( target == "ND3" ) return Numerator_ND3_DF_Thru_PF_Error;
	    else return 0.0;
	}
	double getNACountsDFThruPF(string target) const{
	    if( target == "NH3" ) return NA_Counts_NH3_DF_Thru_PF;
	    else if( target == "ND3" ) return NA_Counts_ND3_DF_Thru_PF;
	    else return 0.0;
	}
	double getNACountsDFThruPFError(string target) const{
	    if( target == "NH3" ) return NA_Counts_NH3_DF_Thru_PF_Error;
	    else if( target == "ND3" ) return NA_Counts_ND3_DF_Thru_PF_Error;
	    else return 0.0;
	}
	double getAllTheory(string target) const{
	    if( target == "NH3" ) return AllTheoryNH3;
	    else if( target == "ND3" ) return AllTheoryND3;
	    else return 0.0;
	}
	double getAllRaw(string target) const{
	    if( target == "NH3" ) return AllRawNH3;
	    else if( target == "ND3" ) return AllRawND3;
	    else return 0.0;
	}
	double getAllRawErr(string target) const{
	    if( target == "NH3" ) return AllRawNH3Err;
	    else if( target == "ND3" ) return AllRawND3Err;
	    else return 0.0;
	}
	double getAllRawNoFC(string target) const{
	    if( target == "NH3" ) return AllRawNoFCNH3;
	    else if( target == "ND3" ) return AllRawNoFCND3;
	    else return 0.0;
	}
	double getAllRawNoFCErr(string target) const{
	    if( target == "NH3" ) return AllRawNoFCNH3Err;
	    else if( target == "ND3" ) return AllRawNoFCND3Err;
	    else return 0.0;
	}

	double getA1Theory(string target) const{
	    if( target == "NH3" ) return A1_Theory_NH3; //DepolFactorNH3*(A1_Theory_NH3 + EtaNH3*A2NH3);
	    else if( target == "ND3" ) return A1_Theory_ND3; //DepolFactorND3*(A1_Theory_ND3 + EtaND3*A2ND3);
	    else return 0.0;
	}
	bool IsInBin(double qmid, double xmid) const{
	    if( qmid < Q2_Max && qmid > Q2_Min && xmid < X_Max && xmid > X_Min ) return true;
	    else return false;
	}
	double getNm(string target) const{
		if(target == "NH3") return NH3_Counts.getNm();
		else if(target == "ND3") return ND3_Counts.getNm();
		else if(target == "CH2") return CH2_Counts.getNm();
		else if(target == "CD2") return CD2_Counts.getNm();
		else if(target == "C") return C_Counts.getNm();
		else if(target == "ET") return ET_Counts.getNm();
		else if(target == "F") return F_Counts.getNm();
		else if(target == "Input") return Input_Loop.getNm();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNp(string target) const{
		if(target == "NH3") return NH3_Counts.getNp();
		else if(target == "ND3") return ND3_Counts.getNp();
		else if(target == "CH2") return CH2_Counts.getNp();
		else if(target == "CD2") return CD2_Counts.getNp();
		else if(target == "C") return C_Counts.getNp();
		else if(target == "ET") return ET_Counts.getNp();
		else if(target == "F") return F_Counts.getNp();
		else if(target == "Input") return Input_Loop.getNp();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNt(string target) const{
		if(target == "NH3") return NH3_Counts.getNt();
		else if(target == "ND3") return ND3_Counts.getNt();
		else if(target == "CH2") return CH2_Counts.getNt();
		else if(target == "CD2") return CD2_Counts.getNt();
		else if(target == "C") return C_Counts.getNt();
		else if(target == "ET") return ET_Counts.getNt();
		else if(target == "F") return F_Counts.getNt();
		else if(target == "Input") return Input_Loop.getNt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNm(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNm();
		else if(target == "ND3") return ND3_Counts.getNormNm();
		else if(target == "CH2") return CH2_Counts.getNormNm();
		else if(target == "CD2") return CD2_Counts.getNormNm();
		else if(target == "C") return C_Counts.getNormNm();
		else if(target == "ET") return ET_Counts.getNormNm();
		else if(target == "F") return F_Counts.getNormNm();
		else if(target == "Input") return Input_Loop.getNormNm();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNp(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNp();
		else if(target == "ND3") return ND3_Counts.getNormNp();
		else if(target == "CH2") return CH2_Counts.getNormNp();
		else if(target == "CD2") return CD2_Counts.getNormNp();
		else if(target == "C") return C_Counts.getNormNp();
		else if(target == "ET") return ET_Counts.getNormNp();
		else if(target == "F") return F_Counts.getNormNp();
		else if(target == "Input") return Input_Loop.getNormNp();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getNormNt(string target) const{
		if(target == "NH3") return NH3_Counts.getNormNt();
		else if(target == "ND3") return ND3_Counts.getNormNt();
		else if(target == "CH2") return CH2_Counts.getNormNt();
		else if(target == "CD2") return CD2_Counts.getNormNt();
		else if(target == "C") return C_Counts.getNormNt();
		else if(target == "ET") return ET_Counts.getNormNt();
		else if(target == "F") return F_Counts.getNormNt();
		else if(target == "Input") return Input_Loop.getNormNt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getAllNt() const{
		double nh3nt = NH3_Counts.getNt(); double nd3nt = ND3_Counts.getNt();
		double ch2nt = CH2_Counts.getNt(); double cd2nt = CD2_Counts.getNt();
		double carnt = C_Counts.getNt(); double empnt = ET_Counts.getNt();
		double foint = F_Counts.getNt();
		return nh3nt + nd3nt + ch2nt + cd2nt + carnt + empnt + foint;
	}
	double getFCm(string target) const{
		if(target == "NH3") return NH3_Counts.getFCm();
		else if(target == "ND3") return ND3_Counts.getFCm();
		else if(target == "CH2") return CH2_Counts.getFCm();
		else if(target == "CD2") return CD2_Counts.getFCm();
		else if(target == "C") return C_Counts.getFCm();
		else if(target == "ET") return ET_Counts.getFCm();
		else if(target == "F") return F_Counts.getFCm();
		else if(target == "Input") return Input_Loop.getFCm();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getFCp(string target) const{
		if(target == "NH3") return NH3_Counts.getFCp();
		else if(target == "ND3") return ND3_Counts.getFCp();
		else if(target == "CH2") return CH2_Counts.getFCp();
		else if(target == "CD2") return CD2_Counts.getFCp();
		else if(target == "C") return C_Counts.getFCp();
		else if(target == "ET") return ET_Counts.getFCp();
		else if(target == "F") return F_Counts.getFCp();
		else if(target == "Input") return Input_Loop.getFCp();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}
	double getFCt(string target) const{
		if(target == "NH3") return NH3_Counts.getFCt();
		else if(target == "ND3") return ND3_Counts.getFCt();
		else if(target == "CH2") return CH2_Counts.getFCt();
		else if(target == "CD2") return CD2_Counts.getFCt();
		else if(target == "C") return C_Counts.getFCt();
		else if(target == "ET") return ET_Counts.getFCt();
		else if(target == "F") return F_Counts.getFCt();
		else if(target == "Input") return Input_Loop.getFCt();
		else{
		    cout <<"ERROR: Couldn't find target type. Return zero counts\n";
		    return 0.0;
		}
	}

	double getAllFCt() const{
		double nh3fct = NH3_Counts.getFCt(); double nd3fct = ND3_Counts.getFCt();
		double ch2fct = CH2_Counts.getFCt(); double cd2fct = CD2_Counts.getFCt();
		double carfct = C_Counts.getFCt(); double empfct = ET_Counts.getFCt();
		double foifct = F_Counts.getFCt();
		return nh3fct + nd3fct + ch2fct + cd2fct + carfct + empfct + foifct;
	}
	// Used for making count ratio plots
	double getCountRatio(string target1, string target2){
		double targ1NormCounts = 0; double targ2NormCounts = 0;

		if(target1 == "NH3") targ1NormCounts = NH3_Counts.getNormNt();
		else if(target1 == "ND3") targ1NormCounts = ND3_Counts.getNormNt();
		else if(target1 == "CH2") targ1NormCounts = CH2_Counts.getNormNt();
		else if(target1 == "CD2") targ1NormCounts = CD2_Counts.getNormNt();
		else if(target1 == "C") targ1NormCounts = C_Counts.getNormNt();
		else if(target1 == "ET") targ1NormCounts = ET_Counts.getNormNt();
		else if(target1 == "F") targ1NormCounts = F_Counts.getNormNt();
		if(target2 == "NH3") targ2NormCounts = NH3_Counts.getNormNt();
		else if(target2 == "ND3") targ2NormCounts = ND3_Counts.getNormNt();
		else if(target2 == "CH2") targ2NormCounts = CH2_Counts.getNormNt();
		else if(target2 == "CD2") targ2NormCounts = CD2_Counts.getNormNt();
		else if(target2 == "C") targ2NormCounts = C_Counts.getNormNt();
		else if(target2 == "ET") targ2NormCounts = ET_Counts.getNormNt();
		else if(target2 == "F") targ2NormCounts = F_Counts.getNormNt();

		if( targ1NormCounts != 0 && targ2NormCounts != 0 ){
			return targ1NormCounts / targ2NormCounts;
		}

		return 0;
	}
	double getCountRatioErr(string target1, string target2){
		double targ1Counts = 0; double targ2Counts = 0;
		double targ1FCt = 0; double targ2FCt = 0;

		if(target1 == "NH3") targ1Counts = NH3_Counts.getNt();
		else if(target1 == "ND3") targ1Counts = ND3_Counts.getNt();
		else if(target1 == "CH2") targ1Counts = CH2_Counts.getNt();
		else if(target1 == "CD2") targ1Counts = CD2_Counts.getNt();
		else if(target1 == "C") targ1Counts = C_Counts.getNt();
		else if(target1 == "ET") targ1Counts = ET_Counts.getNt();
		else if(target1 == "F") targ1Counts = F_Counts.getNt();
		if(target2 == "NH3") targ2Counts = NH3_Counts.getNt();
		else if(target2 == "ND3") targ2Counts = ND3_Counts.getNt();
		else if(target2 == "CH2") targ2Counts = CH2_Counts.getNt();
		else if(target2 == "CD2") targ2Counts = CD2_Counts.getNt();
		else if(target2 == "C") targ2Counts = C_Counts.getNt();
		else if(target2 == "ET") targ2Counts = ET_Counts.getNt();
		else if(target2 == "F") targ2Counts = F_Counts.getNt();

		if(target1 == "NH3") targ1FCt = NH3_Counts.getFCt();
		else if(target1 == "ND3") targ1FCt = ND3_Counts.getFCt();
		else if(target1 == "CH2") targ1FCt = CH2_Counts.getFCt();
		else if(target1 == "CD2") targ1FCt = CD2_Counts.getFCt();
		else if(target1 == "C") targ1FCt = C_Counts.getFCt();
		else if(target1 == "ET") targ1FCt = ET_Counts.getFCt();
		else if(target1 == "F") targ1FCt = F_Counts.getFCt();
		if(target2 == "NH3") targ2FCt = NH3_Counts.getFCt();
		else if(target2 == "ND3") targ2FCt = ND3_Counts.getFCt();
		else if(target2 == "CH2") targ2FCt = CH2_Counts.getFCt();
		else if(target2 == "CD2") targ2FCt = CD2_Counts.getFCt();
		else if(target2 == "C") targ2FCt = C_Counts.getFCt();
		else if(target2 == "ET") targ2FCt = ET_Counts.getFCt();
		else if(target2 == "F") targ2FCt = F_Counts.getFCt();

		if( targ1Counts != 0 && targ2Counts != 0 && targ1FCt != 0 && targ2FCt != 0 ){
			double targ1Err = sqrt( targ1Counts )/targ1FCt;
			double targ2Err = sqrt( targ2Counts )/targ2FCt;
			double targ1Norm = targ1Counts / targ1FCt; double targ2Norm = targ2Counts / targ2FCt;
			//return sqrt( targ1Err*targ1Err/(targ2Counts*targ2Counts) + targ2Err*targ2Err*targ1Counts*targ1Counts/pow(targ2Counts,4));
			return sqrt( (targ1Norm/(targ2Norm*targ2Norm*targ1FCt)) + ((targ1Norm*targ1Norm)/(pow(targ2Norm,3)*targ2FCt)) );
		}
		
		return 0;

	}

	string ReturnBinContents() const{
		string object = to_string(Q2_Min)+"	"+to_string(Q2_Max)+"	";
		object += to_string(X_Min)+"	"+to_string(X_Max)+"	";
		object += to_string(Input_Loop.getNm()) +"	"+to_string(Input_Loop.getNp())+"	";
		object += to_string(Input_Loop.getFCm())+"	"+to_string(Input_Loop.getFCp())+"\n";
		return object;
	}
	void Print() const{
		cout << "=================================== Bin: Q^2 = "<<Q2_Mid <<" X = "<< X_Mid <<" ===================================\n";
		NH3_Counts.Print(); ND3_Counts.Print();
		CH2_Counts.Print(); CD2_Counts.Print();
		C_Counts.Print(); ET_Counts.Print();
		F_Counts.Print();
		cout << "DF_NH3 = "<<DF_NH3 <<" +/- "<<ErrDF_NH3<< endl;
		cout << "DF_ND3 = "<<DF_ND3 <<" +/- "<<ErrDF_ND3<< endl;
		cout << "=========================================================================================================\n";
	}
	// Destructor
	~Bin() = default;

};

// This function is used in initializing the DataSet class to create a 2x2 matrix of bins
vector<vector<Bin>> SetBinWidths(){
	// In this 2x2 matrix, each column corresponds to one bin in Q2, and
	// each row is a bin in Bjorken x.
	vector<vector<Bin>> AllBins;
	//cout << "Now making the bin widths\n";
	for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){ // Loop over Q2 bins
		double qmin = Q2_Bin_Bounds[i];
		double qmax = Q2_Bin_Bounds[i+1];
		vector<Bin> theseQ2Bins;
		//cout << qmin <<" "<< qmax;
		for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){ // Loop over X bins
			Bin thisBin;
			double xmin = X_Bin_Bounds[j];
			double xmax = X_Bin_Bounds[j+1];
			//cout <<" "<< xmin <<" "<< xmax << endl;
			thisBin.SetQ2Bins(qmin, qmax);
			thisBin.SetXBins(xmin, xmax);
			theseQ2Bins.push_back(thisBin);
		}
		AllBins.push_back(theseQ2Bins);
	}
	// Now initialize the theoretical values for the double-spin asymmetry
	string fPath = THIS_DIR + "Asymmetry_Parameterizations.txt";
	ifstream allin( fPath.c_str() );
	if( !allin.fail() ){
	    string line;
	    getline(allin,line); // Throw away header row
	    while( getline(allin, line) ){
	        stringstream sin(line);
		double xbin, Q2bin, W2, Q2, xmid, Q2min, Q2max, Df_NH3, Err_Df_NH3, Q2_avg, W, nu, y, Eprime, theta, eps, D, eta, F1, F2, R, A1, A2, g1, g2;
		sin >> xbin>>Q2bin>>W2>>Q2>>xmid>>Q2min>>Q2max>>Df_NH3>>Err_Df_NH3>>Q2_avg>>W>>nu>>y>>Eprime>>theta>>eps>>D>>eta>>F1>>F2>>R>>A1>>A2>>g1>>g2;
		double qmid = (Q2min+Q2max) / 2.0;
		double ALL_Theory = D*(A1+eta*A2); // Theoretical value of double spin asymmetry for this bin
		//cout << qmid <<" "<< xmid <<" "<< ALL_Theory << endl;
		for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		    for(size_t j=0; j<AllBins[i].size(); j++){
		        if( AllBins[i][j].IsInBin( qmid, xmid )  ){
		            AllBins[i][j].SetAllTheory( ALL_Theory, "NH3" );
			    AllBins[i][j].SetKinematicFactors( D, eta, A1, A2, "NH3" );
			    //AllBins[i][j].SetAllTheory( ALL_Theory_ND3, "ND3" );
			    // Placeholder for any future ND3 incorporation
		        }
		    }
		}
	    }
	}
	else{
	    cout << "Failed to read in A_ll theory values. Check path to 'Asymmetry_Parameterizations.txt' and try again.\n";
	}
	allin.close();

	return AllBins;
}


class DataSet{

    private:
	vector<vector<Bin>> AllBins;
	Bin NullBin; // Used to pass as an error parameter
    public:
	// Constructors
	DataSet() : AllBins( SetBinWidths() ){}

	// Mutators
	// The following function is used to construct an array of bins from input data
	void AddToBins( string filePath, string targtype ){
	    ifstream fin(filePath.c_str());
	    if( !fin.fail() ){
		string line;
		getline(fin,line); // throw away header row
		bool FCchargeIsSet = false; // Tracks whether or not the FC charge has been appended to the data set
		//Q2_Min   Q2_Max   X_Min   X_Max   Nm_Counts   Np_Counts   FCm_Charge   FCp_Charge
		while(getline(fin,line)){
		    stringstream sin(line);
		    double q2_min, q2_max, x_min, x_max, nm_counts, np_counts, fcm_charge, fcp_charge;
		    sin >> q2_min >> q2_max >> x_min >> x_max >> nm_counts >> np_counts >> fcm_charge >> fcp_charge;
		    if( !FCchargeIsSet ){ 
			for(size_t i=0; i<AllBins.size(); i++){
			    for(size_t j=0; j<AllBins[i].size(); j++){
				AllBins[i][j].AddFCCharge(fcm_charge, fcp_charge, targtype);
			    }
			}
			FCchargeIsSet = true;
		    }
		    double qmid = (q2_min+q2_max)/2.0; double xmid = (x_min+x_max)/2.0;
		    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		        for(size_t j=0; j<AllBins[i].size(); j++){
			    if( AllBins[i][j].IsInBin( qmid, xmid )  ){
			        AllBins[i][j].AddCounts( nm_counts, np_counts, targtype );
			    }
			}
		    }
		}
	    }
	    else{
		//cout << "Couldn't open file path " << filePath <<".\n";
		//cout << "Please check path and try again.\n";
	    }
	    fin.close();
	}
	void SlotCountsIntoBin(double qval, double xval, double hel_state){ // Used in Get_DF script
	    // Helicity state corresponds to either Nm or Np, corrected for HWP and target polarization sign
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qval, xval ) ){
		      if( hel_state == 1 ) AllBins[i][j].AddCounts(1.0, 0.0, "Input");
		      else if( hel_state == -1 ) AllBins[i][j].AddCounts(0.0, 1.0, "Input");
		    }
		}
	    }
	}
	void SlotChargeIntoBin(double fcm, double fcp){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].AddFCCharge( fcm, fcp, "Input" );
		}
	    }    
	}
	void NormalizeAllCounts(){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].NormalizeAllCounts();
		}
	    }
	}
	void CalculateDF(){ // This also normalizes the counts when called
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculateDF();
		}
	    }
	}
	void CalculatePF(){ // This also normalizes the counts when called
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculatePF();
		}
	    }
	}
	void CalculateDFThruPF( string targetType ){ // Calculates the DF using the average of all the measured PF values
	    // Make sure that the PF is already calculated first:
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculatePF();
		}
	    }
	    // Now get the weighted average of the PF for the dataset
	    double totalPF = 0; double denominator = 0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    double thisPF = 0; double thisPFerr = 0;
		    if( targetType == "NH3" ){
		        thisPF = AllBins[i][j].getPF_bath_NH3();
		        thisPFerr = AllBins[i][j].getErrPF_bath_NH3();
		        if( thisPF != 0.0 && thisPFerr != 0.0 ){ 
			    totalPF += thisPF / pow(thisPFerr,2) ;
			    denominator += 1.0 / pow(thisPFerr,2);
		        }
		    }
		    else if( targetType == "ND3" ){
		        thisPF = AllBins[i][j].getPF_bath_ND3();
		        thisPFerr = AllBins[i][j].getErrPF_bath_ND3();
		        if( thisPF != 0.0 && thisPFerr != 0.0 ){ 
			    totalPF += thisPF / pow(thisPFerr,2) ;
			    denominator += 1.0 / pow(thisPFerr,2);
			}
		    }
		}
	    }
	    // Now take the new weighted average PF and write it to each of the bins
	    if( denominator != 0.0 ){
	        double avgPF = totalPF / denominator;
		double avgPFerr = 1.0 / sqrt( denominator );
	        cout << "Average PF for "<< targetType<<" data set is: "<< avgPF <<" "<< avgPFerr << endl;
	        // Now use this value to set the values of DF_Thru_PF with the new value
	        for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		        AllBins[i][j].CalculateDFThruPF( avgPF, targetType );
		    }
	        }
	    }
	    else{
		cout <<"ERROR: Zero passed into weighted avg. PF calculation; new PF values are not set.\n";
	    }
	}
	// Set the raw double-spin asymmetry for NH3
	void CalculateAllRaw(){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculateAllRaw();
		}
	    }
	}
	void SetDFs( double qval, double xval, double df, double dferr ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    double Qmid = AllBins[i][j].getBinQ2(); double Xmid = AllBins[i][j].getBinX();
		    //cout << Epoch.getDF_FixdPF_NH3(qmid,xmid) <<" "<< Epoch.getErrDF_FixedPF_NH3(qmid,xmid) << endl;
		    if( AllBins[i][j].IsInBin( qval, xval ) ){
		        AllBins[i][j].SetDF_FixedPF_NH3( df, dferr  ); 
		    }
		}
	    }
	}

	// Accessors
	Bin getThisBin(double qmid, double xmid) const{
	    //Bin& thisBin;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j];
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    // In the fail case, return an empty bin
	    return NullBin;
	}
	double getAvgPF() const{
	    double total = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		  if( AllBins[i][j].getPF_bath_NH3() != 0.0 ){
		    total += AllBins[i][j].getPF_bath_NH3();
		    counts++;
		  }
		}
	    }
	    if(counts != 0.0 ) return total/counts;
	    else return 0.0;
	}
	double getNumeratorDFThruPF(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNumeratorDFThruPF( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNumeratorDFThruPFError(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNumeratorDFThruPFError( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNACountsDFThruPF(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNACountsDFThruPF( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getNACountsDFThruPFError(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getNACountsDFThruPFError( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	// Accessors for the dilution factors
	double getDF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getDF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getDF_FixedPF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_FixedPF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_FixedPF_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_FixedPF_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getDF_FixedPF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getDF_FixedPF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrDF_FixedPF_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrDF_FixedPF_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

	// Accessors for the bath and cell packing fractions for NH3
	double getPF_bath_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_bath_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_bath_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_bath_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getPF_cell_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_cell_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_cell_NH3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_cell_NH3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	// Accessors for the bath and cell packing fractions for ND3
	double getPF_bath_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_bath_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_bath_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_bath_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getPF_cell_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getPF_cell_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getErrPF_cell_ND3(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getErrPF_cell_ND3();
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllTheory(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllTheory(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRaw(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRaw(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRawErr(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawErr(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRawNoFC(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawNoFC(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	double getAllRawNoFCErr(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRawNoFCErr(target);
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}

        // Getting the average PF for the NH3 targets
	double getAveragePF_bath_NH3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_bath_NH3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_bath_NH3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getAveragePF_cell_NH3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_cell_NH3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_cell_NH3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	// Getting the average PF for the ND3 targets
	double getAveragePF_bath_ND3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_bath_ND3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_bath_ND3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getAveragePF_cell_ND3() const{
	    double avgPF = 0.0; double counts = 0.0;
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].getPF_cell_ND3() > 0.0 ){
			avgPF += AllBins[i][j].getPF_cell_ND3();
			counts++;
		    }
		}
	    }
	    if( counts > 0.0 ) return avgPF / counts;
	    else return 0.0;
	}
	double getCountRatio( double qmid, double xmid, string targ1, string targ2 ){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin(qmid,xmid) ){
			return AllBins[i][j].getCountRatio(targ1, targ2);
		    }
		}
	    }
	    return 0;
	}
	double getCountRatioErr( double qmid, double xmid, string targ1, string targ2 ){
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin(qmid,xmid) ){
			return AllBins[i][j].getCountRatioErr(targ1, targ2);
		    }
		}
	    }
	    return 0;
	}

	// Prints the total counts, FC charge, and normalized counts for the entire data set
	// for every available target type.
	void PrintTotalCounts() const{
	    double NA = 0; double FCA = 0; // NH3 counts and FC charge across all bins
	    double ND = 0; double FCD = 0; // ND3
	    double NCH= 0; double FCCH= 0; // CH2
	    double NCD= 0; double FCCD= 0; // CD2
	    double NC = 0; double FCC = 0; // Carbon
	    double NET= 0; double FCET= 0; // Empty target w/ LHe bath
	    double NF = 0; double FCF = 0; // Empty target w/o LHe (foils only)
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    NA += AllBins[i][j].getNt("NH3"); FCA = AllBins[i][j].getFCt("NH3");
		    ND += AllBins[i][j].getNt("ND3"); FCD = AllBins[i][j].getFCt("ND3");
		    NCH+= AllBins[i][j].getNt("CH2"); FCCH= AllBins[i][j].getFCt("CH2");
		    NCD+= AllBins[i][j].getNt("CD2"); FCCD= AllBins[i][j].getFCt("CD2");
		    NC += AllBins[i][j].getNt("C");   FCC = AllBins[i][j].getFCt("C");
		    NET+= AllBins[i][j].getNt("ET");  FCET= AllBins[i][j].getFCt("ET");
		    NF += AllBins[i][j].getNt("F");   FCF = AllBins[i][j].getFCt("F");
		}
	    }
	    cout << "All Counts for this data set:\n";
	    if( FCA != 0 ) cout << "NH3 Counts: NA = "<< NA <<", FCA = "<< FCA <<", NormNA = " << NA/FCA << endl;
	    else cout << "NH3 Counts: NA = "<< NA <<", FCA = "<< FCA << endl;
	    if( FCD != 0 ) cout << "ND3 Counts: ND = "<< ND <<", FCD = "<< FCD <<", NormND = " << ND/FCD << endl;
	    else cout << "ND3 Counts: ND = "<< ND <<", FCD = "<< FCD << endl;
	    if( FCCH != 0 ) cout << "CH2 Counts: NCH = "<< NCH <<", FCCH = "<< FCCH <<", NormNCH = " << NCH/FCCH << endl;
	    else cout << "CH2 Counts: NCH = "<< NCH <<", FCCH = "<< FCCH << endl;
	    if( FCCD != 0 ) cout << "CD2 Counts: NCD = "<< NCD <<", FCCD = "<< FCCD <<", NormNCD = " << NCD/FCCD << endl;
	    else cout << "CD2 Counts: NCD = "<< NCD <<", FCCD = "<< FCCD << endl;
	    if( FCC != 0 ) cout << "Carbon Counts: NC = "<< NC <<", FCC = "<< FCC <<", NormNC = " << NC/FCC << endl;
	    else cout << "Carbon Counts: NC = "<< NC <<", FCC = "<< FCC << endl;
	    if( FCET != 0 ) cout << "ET Counts: NET = "<< NET <<", FCET = "<< FCET <<", NormNET = " << NET/FCET << endl;
	    else cout << "ET Counts: NET = "<< NET <<", FCET = "<< FCET << endl;
	    if( FCF != 0 ) cout << "Foil Counts: NF = "<< NF <<", FCF = "<< FCF <<", NormNF = " << NF/FCF << endl;
	    else cout << "Foil Counts: NF = "<< NF <<", FCF = "<< FCF << endl;

	}
	void Print() const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].Print();
		}
	    }    
	}

	void WriteToCSV(string outFilePath){ // Creates comma delimited file
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		// Write the titles
		fout<<"Bin_Q2,Bin_X,";
		fout<<"DF_NH3,Err_DF_NH3,DF_Fixed_PF_NH3,ErrDF_FixedPF_NH3,PF_bath_NH3,ErrPF_bath_NH3,PF_cell_NH3,ErrPF_cell_NH3,";
		fout<<"DF_ND3,Err_DF_ND3,DF_Fixed_PF_ND3,ErrDF_FixedPF_ND3,PF_bath_ND3,ErrPF_bath_ND3,PF_cell_ND3,ErrPF_cell_ND3,";
		fout<<"FC_total,N_NH3,N_C,N_CH2,N_ET,N_F,N_CD2,FC_NH3,FC_C,FC_CH2,FC_ET,FC_F,FC_CD2,";
		fout<<"normN_normNH3,normN_C,normN_CH2,normN_ET,normN_F,normN_CD2\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0 || AllBins[i][j].getDF_FixedPF_NH3() != 0.0 || AllBins[i][j].getDF_FixedPF_ND3() ){
			fout<<AllBins[i][j].getBinQ2()<<","<<AllBins[i][j].getBinX()<<",";
			fout<<AllBins[i][j].getDF_NH3()<<","<<AllBins[i][j].getErrDF_NH3()<<",";
			fout<<AllBins[i][j].getDF_FixedPF_NH3()<<","<<AllBins[i][j].getErrDF_FixedPF_NH3()<<",";
			fout<<AllBins[i][j].getPF_bath_NH3()<<","<<AllBins[i][j].getErrPF_bath_NH3()<<",";
			fout<<AllBins[i][j].getPF_cell_NH3()<<","<<AllBins[i][j].getErrPF_cell_NH3()<<",";

			fout<<AllBins[i][j].getDF_ND3()<<","<<AllBins[i][j].getErrDF_ND3()<<",";
			fout<<AllBins[i][j].getDF_FixedPF_ND3()<<","<<AllBins[i][j].getErrDF_FixedPF_ND3()<<",";
			fout<<AllBins[i][j].getPF_bath_ND3()<<","<<AllBins[i][j].getErrPF_bath_ND3()<<",";
			fout<<AllBins[i][j].getPF_cell_ND3()<<","<<AllBins[i][j].getErrPF_cell_ND3()<<",";

			fout<<AllBins[i][j].getAllFCt()<<",";
			fout<<AllBins[i][j].getNt("NH3")<<","<<AllBins[i][j].getNt("C")<<","<<AllBins[i][j].getNt("CH2")<<","<<AllBins[i][j].getNt("ET")<<",";
			fout<<AllBins[i][j].getNt("F")<<","<<AllBins[i][j].getNt("CD2")<<",";
			fout<<AllBins[i][j].getFCt("NH3")<<","<<AllBins[i][j].getFCt("C")<<","<<AllBins[i][j].getFCt("CH2")<<","<<AllBins[i][j].getFCt("ET")<<",";
			fout<<AllBins[i][j].getFCt("F")<<","<<AllBins[i][j].getFCt("CD2")<<",";
			fout<<AllBins[i][j].getNormNt("NH3")<<","<<AllBins[i][j].getNormNt("C")<<","<<AllBins[i][j].getNormNt("CH2")<<","<<AllBins[i][j].getNormNt("ET")<<",";
			fout<<AllBins[i][j].getNormNt("F")<<","<<AllBins[i][j].getNormNt("CD2")<<"\n";
		      }
		    }
	        }
	    }
	    fout.close();
	}
	void WriteRawAsymsTXT(string outFilePath){ // Creates an output file for all the raw asymmetries in a tab delimited format
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		// Write titles
		fout <<"Bin_Q2 Bin_X, AllrawNH3, ErrAllrawNH3, AllrawND3, ErrAllrawND3\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0){	
		      }
		    }
		}
	    }
	}
	void WriteToTXT(string outFilePath){ // Creates tab delimited file
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		// Write the titles
		fout<<"Bin_Q2,Bin_X,DF_NH3,Err_DF_NH3,DF_ND3,Err_DF_ND3,PF_bath_NH3,ErrPF_bath_NH3,FC_total,N_NH3,N_C,N_CH2,N_ET,N_F,N_CD2,FC_NH3,FC_C,FC_CH2,FC_ET,FC_F,FC_CD2\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0){
			fout<<AllBins[i][j].getBinQ2()<<"	"<<AllBins[i][j].getBinX()<<"	"<<AllBins[i][j].getDF_NH3()<<"	"<<AllBins[i][j].getErrDF_NH3()<<"	";
			fout<<AllBins[i][j].getDF_ND3()<<"	"<<AllBins[i][j].getErrDF_ND3()<<"	"<<AllBins[i][j].getPF_bath_NH3()<<"	"<<AllBins[i][j].getErrPF_bath_NH3();
			fout<<AllBins[i][j].getAllFCt()<<"	";
			fout<<AllBins[i][j].getNt("NH3")<<"	"<<AllBins[i][j].getNt("C")<<"	"<<AllBins[i][j].getNt("CH2")<<"	"<<AllBins[i][j].getNt("ET")<<"	";
			fout<<AllBins[i][j].getNt("F")<<"	"<<AllBins[i][j].getNt("CD2")<<"	";
			fout<<AllBins[i][j].getFCt("NH3")<<"	"<<AllBins[i][j].getFCt("C")<<"	"<<AllBins[i][j].getFCt("CH2")<<"	"<<AllBins[i][j].getFCt("ET")<<"	";
			fout<<AllBins[i][j].getFCt("F")<<"	"<<AllBins[i][j].getFCt("CD2")<<"\n";
		      }
		    }
	        }
	    }
	    fout.close();
	}
	void WriteRawDFInputs(string outFilePath){
	    ofstream fout(outFilePath);
	    if( !fout.fail() ){
		fout <<" Q2_Min   Q2_Max   X_Min   X_Max   Nm_Counts   Np_Counts   FCm_Charge   FCp_Charge\n";
		for(size_t i=0; i<AllBins.size(); i++){
		    for(size_t j=0; j<AllBins[i].size(); j++){
			fout << AllBins[i][j].ReturnBinContents();
		    }
		}
	    }
	    fout.close();
	}
	// Destructor
	~DataSet() = default;	
};

// This class is used for calculating the beam-target polarization on a run-by-run basis
class PbPt{

    private:
	vector<vector<Bin>> AllBins;
	int RunNumber = 0;
	double NMR_Tpol = 0;
	string Target = "N/A";
	double BT_Pol = 0.0;
	double BT_Pol_Err = 0.0;
	int EpochNum = 0;
	//double AllPhys = 0.0;
	//double AllPhysErr = 0.0;
    public:
	// Constructors
	PbPt() : AllBins( SetBinWidths() ){}

	// Mutators
	void SetRunInfo(int run, string target ){
	    RunNumber = run; Target = target;
	}
	// This function reads in from a particular run so PbPt can be calculated for it
	void ReadInThisRun(int thisRun, string TARGET_TYPE, RunPeriod& Period){
	    RunNumber = thisRun;
	    Target = TARGET_TYPE;
	    EpochNum = Period.getEpochNum( thisRun );
	    string infile = "Text_Files/"+Target+"_"+to_string(thisRun)+"_DF_Data.txt";
	    ifstream fin(infile.c_str());
	    if( !fin.fail() ){
		NMR_Tpol = Period.getTargetPolarization( thisRun );
		string line; getline(fin, line); // Throw away header row
		while( getline(fin, line) ){
		    stringstream sin(line);
		    double qmin, qmax, xmin, xmax, Nm, Np, FCm, FCp;
		    sin >> qmin >> qmax >> xmin >> xmax >> Nm >> Np >> FCm >> FCp;
		    double qmid = (qmin + qmax) / 2.0;
		    double xmid = (xmin + xmax) / 2.0;
		    for(size_t i=0; i<AllBins.size(); i++){
			for(size_t j=0; j<AllBins[i].size(); j++){
			    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
				AllBins[i][j].AddCounts( Nm, Np, Target );
				AllBins[i][j].AddFCCharge( FCm, FCp, Target );
				//cout << qmid <<" "<< xmid <<" "<< Nm <<" "<< Np <<" "<< FCm <<" "<< FCp << endl;
			    }
			}
		    }
		}
	    }
	    else{
		cout << "ERROR: Couldn't open file '"<<infile<<"'. Check path and try again.\n";
	    }
	    fin.close();
	}
	// Calculates the raw double-spin asymmetry for this bin
	void CalculateAllRaw(){
	    if( Target == "NH3" || Target == "ND3" ){
	      for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].CalculateAllRaw();
		    /*
	    	    double normNm = AllBins[i][j].getNormNm(Target);
		    double normNp = AllBins[i][j].getNormNp(Target);
	    	    double nm = AllBins[i][j].getNm(); double np = AllBins[i][j].getNp();
	    	    if( normNm > 0.0 && normNp > 0.0 ){ 
			AllRaw = (normNm - normNp) / (normNm + normNp);
			AllRawErr = 0.5*sqrt( (nm + np) / nm*np );
			AllBins[i][j].SetAllRaw( AllRaw, Target );
			AllBins[i][j].SetAllRawErr( AllRawErr, Target );
	    	    }
	    	    else{
			//cout << "ERROR: FC-normalized counts might be zero; check inputs. Set AllRawNH3 = 0\n";
	    	    }
		    */
		}
	      }
	    }
	    else cout <<"ERROR: Target type not properly set: "<< Target << endl;
	}
	void SetTarget( string target ){
	    Target = target;
	}

	// This copies a set of dilution factors from a given epoch; allows for selection between using a fixed PF for the
	// entire epoch or using a unique PF for each bin.
	void SetDFs( DataSet& Epoch ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    double qmid = AllBins[i][j].getBinQ2(); double xmid = AllBins[i][j].getBinX();
		    //cout << Epoch.getDF_FixedPF_NH3(qmid,xmid) <<" "<< Epoch.getErrDF_FixedPF_NH3(qmid,xmid) << endl;
		    AllBins[i][j].SetDF_FixedPF_NH3( Epoch.getDF_FixedPF_NH3(qmid,xmid), Epoch.getErrDF_FixedPF_NH3(qmid,xmid) ); 
		}
	    }
	}
	void SetDFs( double qval, double xval, double df, double dferr ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    double Qmid = AllBins[i][j].getBinQ2(); double Xmid = AllBins[i][j].getBinX();
		    //cout << Epoch.getDF_FixdPF_NH3(qmid,xmid) <<" "<< Epoch.getErrDF_FixedPF_NH3(qmid,xmid) << endl;
		    if( AllBins[i][j].IsInBin( qval, xval ) ){
		        AllBins[i][j].SetDF_FixedPF_NH3( df, dferr  ); 
		    }
		}
	    }
	}

	// After calculating the AllRaw and copying the DFs for each x, Q2 bin, find the value of PbPt for each Q2 bin
	void CalculatePbPt(){
	    //double total = 0.0; double counts = 0.0;
	    vector<double> pbptVals, pbptErrs;
	    // Account for the sign of the target polarization
	    int corrFactor = 0;
	    if( NMR_Tpol < 0 ) corrFactor = -1;
	    else if( NMR_Tpol > 0 ) corrFactor = 1; // Allow for corrFactor to be zero for error tracing purposes

	    double num = 0.0; double denom = 0.0;
	    double errTerm = 0.0; // Not the full error; used in the calculation of the error
	    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		//double total = 0.0; double counts = 0.0;
		for(size_t j=0; j<AllBins[i].size(); j++){ // Loop over x bins for a given Q2 bin
		    double thisDF = AllBins[i][j].getDF_FixedPF_NH3();
		    double thisDFErr = AllBins[i][j].getErrDF_FixedPF_NH3();
		    double AllRaw = AllBins[i][j].getAllRaw("NH3");
		    double AllRawErr = AllBins[i][j].getAllRawErr("NH3");
		    double AllTheory = AllBins[i][j].getAllTheory("NH3");
		    //cout <<"DF = "<< thisDF <<" +- "<< thisDFErr <<" "<< AllRaw <<" +- "<< AllRawErr<<" "<< AllTheory << endl;
		    if( thisDF > 0.0 && AllRaw != 0.0 && AllTheory != 0.0 && thisDFErr > 0.0 && AllRawErr > 0.0 ){
			/*
			num += AllRaw * corrFactor;
			denom += AllTheory*thisDF;
			errTerm += pow((AllRawErr/AllRaw),2) + pow((thisDFErr/thisDF),2);
			*/
			double polErr = sqrt( pow(AllRawErr,2) + pow( (AllRaw*thisDFErr/thisDF),2 )  ) / (thisDF*AllTheory);
			pbptVals.push_back( corrFactor * AllRaw / (thisDF * AllTheory) );
			pbptErrs.push_back( polErr );
			//total += AllRaw / (thisDF * AllTheory);
			//counts++;
		    }
		}
		//double totalPbPt = total / counts;
		//cout << "Beam-target polaization for Q^2 = "<< AllBins[i][0].getBinQ2() <<": PbPt = "<< totalPbPt << endl;
	    }
	    //double totalPbPt = total / counts;
	    //cout << "Beam-target polaization for this run: PbPt = "<< totalPbPt << endl;
	    // Now calculate the weighted average of the beam-target polarizations
	    //double num = 0.0; double denom = 0.0;
	    for(size_t i=0; i<pbptVals.size(); i++){
		num += pbptVals[i] / pow( pbptErrs[i], 2);
		denom += 1.0 / pow( pbptErrs[i], 2);
	    }
	    
	    if( denom != 0.0 ){
	        cout << "Beam-target polaization for run: PbPt = "<< num / denom <<" +- "<< sqrt(1.0 / denom) << endl;
		BT_Pol = num / denom;
		BT_Pol_Err = sqrt( 1.0 / denom );
		//BT_Pol_Err = BT_Pol*0.05;//sqrt( errTerm );
	    }
	    else
		cout << "ERROR: Invalid value of PbPt.\n";
	}
	// This calculates the PbPt using the elastic electron-proton channel that Noemie uses. I'm trying to use this for
	// DIS as well, but let's see how that works.
	void CalculateElasticPbPt(){
	    // Account for the sign of the target polarization
	    int corrFactor = 0;
	    if( NMR_Tpol < 0 ) corrFactor = -1;
	    else if( NMR_Tpol > 0 ) corrFactor = 1; // Allow for corrFactor to be zero for error tracing purposes

	    double numPbPt = 0.0; double denomPbPt = 0.0; // Numerator and denominator terms used for PbPt calculation
	    double numPbPtErr = 0.0; double denomPbPtErr = 0.0; // Numerator and denominator terms used for error in PbPt calculation
							        // denomPbPtErr term needs to be squared at the end of the calculation!!!
	    for(size_t i=0; i<AllBins.size(); i++){ // Q2 bin loop
		//double total = 0.0; double counts = 0.0;
		for(size_t j=0; j<AllBins[i].size(); j++){ // Loop over x bins for a given Q2 bin
		    double thisNm = AllBins[i][j].getNm("NH3"); // Raw counts NOT normalized to FC charge
		    double thisNp = AllBins[i][j].getNp("NH3"); // Same as above
		    double thisFCm= AllBins[i][j].getFCm("NH3");// The FCm charge from the NH3 targets
		    double thisFCp= AllBins[i][j].getFCp("NH3");// Same as above for FCp
		    double thisNormNm = AllBins[i][j].getNormNm("NH3"); // Nm counts normalized to FCm
		    double thisNormNp = AllBins[i][j].getNormNp("NH3"); // Np counts normalized to FCp
		    double thisDF = AllBins[i][j].getDF_FixedPF_NH3();
		    double thisDFErr = AllBins[i][j].getErrDF_FixedPF_NH3();
		    double AllTheory = AllBins[i][j].getAllTheory("NH3");
		    if( thisDF > 0.0 && AllTheory != 0.0 && thisDFErr > 0.0 && thisFCm > 0.0 && thisFCp > 0.0 && thisNm > 0.0 && thisNp > 0.0 ){
		    //if( thisDF > 0.0 && AllTheory != 0.0 && thisDFErr > 0.0 && thisNormNm > 0.0 && thisNormNp > 0.0 ){
			
			numPbPt += corrFactor * thisDF * AllTheory * (thisNm - thisNp);
			denomPbPt += thisDF*thisDF * AllTheory*AllTheory * (thisNm + thisNp);
			numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * ((thisNm/thisFCm) + (thisNp/thisFCp));
			denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(thisNm + thisNp);
			
			/*
			numPbPt += corrFactor * thisDF * AllTheory * (thisNormNm - thisNormNp);
			denomPbPt += thisDF*thisDF * AllTheory*AllTheory * (thisNormNm + thisNormNp);
			numPbPtErr += thisDF*thisDF*AllTheory*AllTheory * ((thisNormNm/thisFCm) + (thisNormNp/thisFCp));
			denomPbPtErr += thisDF*thisDF*AllTheory*AllTheory*(thisNormNm + thisNormNp);
			*/
		    }
		}
	    } 
	    if( denomPbPt != 0.0 && denomPbPtErr != 0.0 ){
		BT_Pol = numPbPt / denomPbPt;
		BT_Pol_Err = sqrt( numPbPtErr )/denomPbPtErr; // Pulled out of square root

	        cout << "Beam-target polaization for run: PbPt = "<< BT_Pol <<" +- "<< BT_Pol_Err << endl;
		//BT_Pol_Err = BT_Pol*0.05;//sqrt( errTerm );
	    }
	    else
		cout << "ERROR: Invalid value of PbPt.\n";
	}

	void SetAllPhys( double allphys, double allphyserr, double qmid, double xmid ){
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ) AllBins[i][j].SetAllPhys( allphys, allphyserr, "NH3" );
		}
	    }
	}
	void CalculateAllPhys( int Run, string Target, DataSet& Epoch, RunPeriod& Period ){
	    // Make sure that everything has been properly calculated.
	    // Read in the dilution factors from the input epoch
	    //cout << "Set DFs\n";
	    this->SetDFs( Epoch );
	    //cout << "Read in run "<< Run <<" "<< Target << endl;
	    this->ReadInThisRun( Run, Target, Period ); // Read in info for this run
	    //cout << "Calculate AllRaw\n";
	    this->CalculateAllRaw(); // Set the raw asymmetry
	    //cout << "Calculate PbPt\n";
	    //this->CalculatePbPt(); // Calculate the beam-target polarization for this run
	    this->CalculateElasticPbPt(); // Calculates PbPt using Noemie's formula
	    // Once we have the raw asymmetries and target polarizations, calculate the A1 values for each bin
	    for( size_t i=0; i<AllBins.size(); i++ ){ // Q2 bin loop
		for( size_t j=0; j<AllBins[i].size(); j++){ // x bin loop
		    Bin thisBin = AllBins[i][j];
		    if( thisBin.getDF_FixedPF_NH3() != 0.0 && thisBin.getDepolFactorNH3() != 0.0 ){
			AllBins[i][j].CalculateAllPhys( BT_Pol, BT_Pol_Err );
			//cout << "A_ll,raw = "<< AllBins[i][j].getAllRaw("NH3") <<" +- "<< AllBins[i][j].getAllRawErr("NH3") << endl;
			//cout << "A_ll,phys = "<< AllBins[i][j].getAllPhysNH3() <<" +- "<< AllBins[i][j].getAllPhysNH3Err() << endl;
		        //double a1 = (thisBin.getAllRaw("NH3")/(BT_Pol*thisBin.getDF_FixedPF_NH3()*thisBin.getDepolFactorNH3())) - thisBin.getEtaNH3()*thisBin.getA2NH3();
		        //double a1_err = a1*0.05; // temporary placeholder
		        //AllBins[i][j].SetA1( a1, a1_err, "NH3" );
		    }
		}
	    }
	}
	void CalculateA1(){
	    // Make sure that everything has been properly calculated.
	    // Read in the dilution factors from the input epoch
	    // Once we have the raw asymmetries and target polarizations, calculate the A1 values for each bin
	    for( size_t i=0; i<AllBins.size(); i++ ){ // Q2 bin loop
		for( size_t j=0; j<AllBins[i].size(); j++){ // x bin loop
		    Bin thisBin = AllBins[i][j];
		    if( thisBin.getDepolFactorNH3() != 0.0 && thisBin.getAllPhysNH3() > 0.0 ){
		        double a1 = (thisBin.getAllPhysNH3()/(thisBin.getDepolFactorNH3())) - thisBin.getEtaNH3()*thisBin.getA2NH3();
		        double a1_err = thisBin.getAllPhysNH3Err()/thisBin.getDepolFactorNH3(); // temporary placeholder
			//cout << "A1 = "<< a1 <<" +- "<< a1_err << endl;
		        AllBins[i][j].SetA1( a1, a1_err, "NH3" );
		    }
		}
	    }
	}


	// Accessors
	double getBT_Pol() const{ return BT_Pol; }
	double getBT_Pol_Err() const{ return BT_Pol_Err; }
	int getEpochNum() const{ return EpochNum; }
	double getAllPhys(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getAllPhysNH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getErrAllPhys(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getAllPhysNH3Err();
		    }
		}
	    }
	    return 0.0;
	}

	double getA1(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getA1NH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getErrA1(double qmid, double xmid) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getErrA1NH3();
		    }
		}
	    }
	    return 0.0;
	}
	double getA1Theory(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){
			return AllBins[i][j].getA1Theory( target );
		    }
		}
	    }
	    return 0.0;
	}
	double getAllRaw(double qmid, double xmid, string target) const{
	    for(size_t i=0; i<AllBins.size(); i++){
	 	for(size_t j=0; j<AllBins[i].size(); j++){
		    if( AllBins[i][j].IsInBin( qmid, xmid ) ){ 
			return AllBins[i][j].getAllRaw( target );
		    }
		}
	    }
	    cout <<"ERROR: Couldn't find bin in range specified: Q2 = "<<qmid<<", X = "<<xmid<<endl;
	    return 0.0;
	}
	void Print() const{
	    for(size_t i=0; i<AllBins.size(); i++){
		for(size_t j=0; j<AllBins[i].size(); j++){
		    AllBins[i][j].Print();
		}
	    }
	}
	// Destructors

};

// Reads in A1 data from previous experiments; used in the Make_A1_NH3_Plots function
// Data in the following form: Compass
vector<TGraphErrors*> Old_A1_Data(){
    
    vector<vector<double>> xVals(7);
    vector<vector<double>> A1data(7);
    vector<vector<double>> A1err(7);
    vector<vector<double>> zeros(7);
    vector<string> Titles = {"Compass","E143","E155","EMC","Hermes","SMC","EG1b"};
    vector<int> Styles = {3,4,8,21,26,30,34};
    string infile = THIS_DIR + "A1p_DIS_Pushpa.txt";
    ifstream fin(infile);
    string line; getline(fin,line); // Throw out header row
    while( getline( fin, line ) ){
	double Expt, x, q2, w2, a1data, a1stat, a1sys, a1err;
	stringstream sin(line);
	sin >> Expt >> x >> q2 >> w2 >> a1data >> a1stat >> a1sys >> a1err;
	if( Expt != 9 && Expt < 11 ){
	    int index; 
	    if( Expt <= 3 ) index = 0;
	    else if( Expt == 10 ) index = 6;
	    else index = Expt - 3;
	    xVals[index].push_back(x);
	    A1data[index].push_back(a1data);
	    A1err[index].push_back(a1err);
	    zeros[index].push_back(0);
	}
    }
    fin.close();
    // Now make the plots for each data set
    vector<TGraphErrors*> OldData;
    TCanvas* DataPlot = new TCanvas("DataPlot","DataPlot",800,600);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.1,0.7,0.5,0.9); leg->SetNColumns(2);
    for(size_t i=0; i<A1data.size(); i++){
	TGraphErrors* gr = new TGraphErrors(xVals[i].size(),xVals[i].data(),A1data[i].data(),zeros[i].data(),A1err[i].data());
	gr->SetMarkerColor(kGray); gr->SetLineColor(kGray);
	gr->SetMarkerStyle(Styles[i]);
	leg->AddEntry(gr,Titles[i].c_str(),"p");
	mg->Add(gr,"p");
	OldData.push_back(gr);
    }
    mg->SetTitle("Old A_{1,p} Data; X; A_{1,p}");
    mg->GetYaxis()->SetRangeUser(0,1.2); 
    DataPlot->cd(); mg->Draw("ap"); leg->Draw("same");
    return OldData;
}

// Makes a plot of both All,phys and A1 for the proton
vector<TCanvas*> Make_A1_NH3_Plots( vector<PbPt>& AllPhysVals, DataSet& Epoch, string sector, string color ){
//vector<TGraphErrors*> Make_A1_NH3_Plots( vector<PbPt>& AllPhysVals, DataSet& Epoch, string sector, string color ){

    PbPt AllA1Vals; AllA1Vals.SetTarget("NH3");

    for( size_t i=0; i<Q2_Bin_Bounds.size()-1; i++ ){ // q2 bin loop
        for( size_t j=0; j<X_Bin_Bounds.size()-1; j++ ){ // x bin loop
	    double thisBinNum = 0.0;
	    double thisBinDenom = 0.0;
	    double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    //cout << "Looking at bin "<< qmid <<" "<< xmid <<":\n";
            for( size_t k=0; k<AllPhysVals.size(); k++){
	        // Loop through each bin and get the average AllPhys
		double thisBinAllPhys = AllPhysVals[k].getAllPhys(qmid,xmid);
		double thisBinAllPhysErr = AllPhysVals[k].getErrAllPhys(qmid,xmid);
	        //cout << "Checking AllPhys for bin "<< qmid <<" "<< xmid <<" "<< thisBinAllPhys <<" +- "<< thisBinAllPhysErr << endl;
		if( thisBinAllPhys > 0 && thisBinAllPhysErr > 0 ){
		    thisBinNum += thisBinAllPhys / pow(thisBinAllPhysErr,2);
		    thisBinDenom += 1.0 / pow(thisBinAllPhysErr,2);
		}
	    }
	    if( thisBinDenom > 0.0 ){
	        double AvgAllPhys = thisBinNum / thisBinDenom;
	        double AvgAllPhysErr = 1.0 / sqrt( thisBinDenom );
	        //cout << "Entering AllPhys for bin "<< qmid <<" "<< xmid <<" "<< AvgAllPhys <<" +- "<< AvgAllPhysErr;
	        AllA1Vals.SetAllPhys( AvgAllPhys, AvgAllPhysErr, qmid, xmid ); 
	        AllA1Vals.CalculateA1();
		
	    }
	}
    }
    // Now make the plot
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgA1_NH3 = new TMultiGraph();
    TMultiGraph* mgAllPhys_NH3 = new TMultiGraph();
    TMultiGraph* mgA1_NH3_Q2Avg = new TMultiGraph();

    TLegend* leg = new TLegend(0.1,0.7,0.5,0.9);
    leg->SetNColumns(3);
    leg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Add this historical data to the plot
    vector<TGraphErrors*> OldA1Plots = Old_A1_Data(); // Historical A1 data
    for(size_t i=0; i<OldA1Plots.size(); i++){ 
	mgA1_NH3->Add(OldA1Plots[i],"p");
	mgA1_NH3_Q2Avg->Add(OldA1Plots[i],"p");
    }

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, a1nh3vals, zeros, a1nh3errs, allphysnh3vals, allphysnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisA1NH3 = AllA1Vals.getA1( qmid, xmid );
	    double thisErrA1NH3 = AllA1Vals.getErrA1( qmid, xmid);
	    double thisAllPhysNH3 = AllA1Vals.getAllPhys( qmid, xmid );
	    double thisAllPhysErrNH3 = AllA1Vals.getErrAllPhys( qmid, xmid );
	    if( thisA1NH3 != 0.0 && thisAllPhysNH3 != 0.0 /*&& thisAllPhysNH3 < 0.5*/ ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		a1nh3vals.push_back( thisA1NH3 );
		a1nh3errs.push_back( thisErrA1NH3 );
		allphysnh3vals.push_back( thisAllPhysNH3 );
		allphysnh3errs.push_back( thisAllPhysErrNH3 );
		//cout << thisA1NH3 <<" +- "<< thisErrA1NH3 << endl;
		//cout << thisAllPhysNH3 <<" +- "<< thisAllPhysErrNH3 << endl;
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), a1nh3vals.data(), zeros.data(), a1nh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    TGraphErrors* gr2 = new TGraphErrors(xbins.size(), xbins.data(), allphysnh3vals.data(), zeros.data(), allphysnh3errs.data() );
	    gr2->SetMarkerColor( palette[pint] );

	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); gr2->SetMarkerStyle(kFullCircle);  }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); gr2->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); gr2->SetMarkerStyle(kFullSquare); }
	    mgA1_NH3->Add( gr, "p" );
	    mgAllPhys_NH3->Add( gr2, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    leg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }
    // Now make the plot for the fixed Q2. For a fixed value in X, take the average over all Q2 bins for that value of X
    vector<double> avgA1vals, Xvals, Zeros, avgA1valsErr, avgA1theory;
    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	double thisXBinA1num = 0.0; double thisXBinA1denom = 0.0;
	double sumA1Theory = 0.0; double counts = 0.0;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisA1NH3 = AllA1Vals.getA1(qmid, xmid);
	    double thisErrA1NH3 = AllA1Vals.getErrA1( qmid, xmid );
	    double thisA1TheoryNH3 = AllA1Vals.getA1Theory( qmid, xmid, "NH3" );
	    if( thisA1NH3 > 0 && thisErrA1NH3 > 0 && thisErrA1NH3 < 0.5 ){
		thisXBinA1num += thisA1NH3 / pow(thisErrA1NH3,2);
		thisXBinA1denom += 1.0 / pow(thisErrA1NH3,2);
		sumA1Theory += thisA1TheoryNH3; counts++;
	    }
	}
	if( thisXBinA1denom != 0 && counts != 0 ){
	    avgA1vals.push_back( thisXBinA1num / thisXBinA1denom );
	    avgA1valsErr.push_back( 1.0/sqrt(thisXBinA1denom)  );
	    Xvals.push_back( xmid ); Zeros.push_back(0);
	    avgA1theory.push_back( sumA1Theory / counts );
	}
    }
    vector<TGraphErrors*> A1Plots;
    TGraphErrors* gr3 = new TGraphErrors(Xvals.size(), Xvals.data(), avgA1vals.data(), Zeros.data(), avgA1valsErr.data() );
    if( color == "Blue" ){
        gr3->SetMarkerColor( kBlue ); gr3->SetMarkerStyle(kFullCircle);
        TGraphErrors* gr4 = new TGraphErrors(Xvals.size(), Xvals.data(), avgA1theory.data(), Zeros.data(), Zeros.data() );
        gr4->SetLineColor( kRed ); gr4->SetLineWidth(4);
        mgA1_NH3_Q2Avg->Add( gr3, "p" );
        mgA1_NH3_Q2Avg->Add( gr4, "l" );
	mgA1_NH3->Add( gr4, "l" );
	A1Plots.push_back( gr3 ); A1Plots.push_back( gr4 );
    }
    else{
        gr3->SetMarkerColor( kBlack ); gr3->SetMarkerStyle(kFullSquare);
        mgA1_NH3_Q2Avg->Add( gr3, "p" );
	A1Plots.push_back( gr3 );
    }
    //mgA1_NH3_Q2Avg->Add( gr3, "p" );
    //mgA1_NH3_Q2Avg->Add( gr4, "l" );


    string title = "A_{1,p}(X,Q^{2}) for "+sector+"; X; A_{1,p}";
    mgA1_NH3->SetTitle(title.c_str());
    mgA1_NH3->GetXaxis()->SetLimits(0,0.8); //mgA1_NH3->GetYaxis()->SetRangeUser(0.0,0.4);
    mgA1_NH3->GetYaxis()->SetRangeUser(0,1.2);
    string title2 = "A_{||,phys} for "+sector+"; X; A_{||,phys}";
    mgAllPhys_NH3->SetTitle(title2.c_str());
    mgAllPhys_NH3->GetXaxis()->SetLimits(0,0.8); //mgAllPhys_NH3->GetYaxis()->SetRangeUser(0.0,0.4);
    mgAllPhys_NH3->GetYaxis()->SetRangeUser(0,1.2);
    string title3 = "A_{1,p}(X) for "+sector+"; X; A_{1,p}";
    mgA1_NH3_Q2Avg->SetTitle(title3.c_str());
    mgA1_NH3_Q2Avg->GetXaxis()->SetLimits(0,0.8);
    mgA1_NH3_Q2Avg->GetYaxis()->SetRangeUser(0,1.2);

    // Make a line showing the SU(6) prediction
    TF1* SU6 = new TF1("SU6","5.0/9.0",-1,1);
    SU6->SetLineWidth(3); SU6->SetLineColor(kBlack);
    SU6->SetLineStyle(2);

    //return mgA1_NH3;
    string c1title = "c1_"+sector;
    string c2title = "c2_"+sector;
    string c3title = "c3_"+sector;

    TCanvas* c1 = new TCanvas(c1title.c_str(),c1title.c_str(),800,600);
    mgA1_NH3->Draw("AP"); leg->Draw("same"); SU6->Draw("same");
    TCanvas* c2 = new TCanvas(c2title.c_str(),c2title.c_str(),800,600);
    mgAllPhys_NH3->Draw("AP"); leg->Draw("same");
    TCanvas* c3 = new TCanvas(c3title.c_str(),c3title.c_str(),800,600);
    mgA1_NH3_Q2Avg->Draw("ap"); SU6->Draw("same");
    vector<TCanvas*> Plots = {c1, c2, c3};
    
    return Plots;

    // return A1Plots;

}

// This makes a TGraphErrors object of all the target polarizations for the runs in a given epoch
// Make sure that the DFs have been calculated for the "DF_Data" set!!!
TGraphErrors* Run_Range_PbPt( vector<int>& Epoch, DataSet& DF_Data, string Target, RunPeriod& Period ){

    vector<double> runs, pbpt, pbptErr, zeros; // Used for plotting
    for(size_t i=0; i<Epoch.size(); i++){
	PbPt thisTpol;
	thisTpol.ReadInThisRun( Epoch[i], Target, Period );
	thisTpol.SetDFs( DF_Data );
	thisTpol.CalculateAllRaw();
	thisTpol.CalculatePbPt();
	double bt = thisTpol.getBT_Pol(); double btErr = thisTpol.getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i] ); pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	}
    }
    if( runs.size() > 0 ){
	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	return PbPt_Vals;
    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TGraphErrors...\n";
	TGraphErrors* empty = new TGraphErrors();
	return empty;
    }

}

// Plots PbPt over several epochs using the elastic method
TMultiGraph* All_Epochs_ElasticPbPt( vector<vector<int>>& Epoch, vector<DataSet>& DF_Data, string Target, RunPeriod& Period ){

    vector<double> runs, pbpt, pbptErr, zeros, NMRpbpt, NMRpbptErr; // Used for plotting
    for(size_t i=0; i<Epoch.size(); i++){ // Loop over epochs
     for(size_t j=0; j<Epoch[i].size(); j++){ // Loop over the runs within an epoch
	PbPt thisTpol;
	cout << "Starting PbPt for run "<< Epoch[i][j] << endl;
	thisTpol.ReadInThisRun( Epoch[i][j], Target, Period );
	thisTpol.SetDFs( DF_Data[i] ); // This assumes, of course, that the epochs in Epoch and DF_Data line up, so be careful of that...
	thisTpol.CalculateElasticPbPt();
	double bt = thisTpol.getBT_Pol(); double btErr = thisTpol.getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i][j] ); pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	    //NMRpbpt.push_back( 0.8 * Period.getTargetPolarization( Epoch[i][j] ) );
	    NMRpbpt.push_back( 0.8 * Period.getOfflineTPol( Epoch[i][j] ));
	    NMRpbptErr.push_back( 0.8 * Period.getOfflineTPolErr( Epoch[i][j] ));
	}
     }
    }
    if( runs.size() > 0 ){
	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	PbPt_Vals->SetMarkerColor(kBlue); PbPt_Vals->SetMarkerStyle(kFullCircle);
	//TGraphErrors* NMR_PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), NMRpbpt.data(), zeros.data(), NMRpbptErr.data() );
	//NMR_PbPt_Vals->SetMarkerColor(kViolet); NMR_PbPt_Vals->SetMarkerStyle(kFullSquare);
	TMultiGraph* mg = new TMultiGraph();
	mg->Add( PbPt_Vals, "p");
	//mg->Add( NMR_PbPt_Vals, "p");
	
	return mg;
    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TMultiGraph...\n";
	TMultiGraph* empty = new TMultiGraph();
	return empty;
    }

}

// Plots PbPt over several epochs
TMultiGraph* All_Epochs_PbPt( vector<vector<int>>& Epoch, vector<DataSet>& DF_Data, string Target, RunPeriod& Period ){

    vector<double> runs, pbpt, pbptErr, zeros, NMRpbpt, NMRpbptErr; // Used for plotting
    for(size_t i=0; i<Epoch.size(); i++){ // Loop over epochs
     for(size_t j=0; j<Epoch[i].size(); j++){ // Loop over the runs within an epoch
	PbPt thisTpol;
	cout << "Starting PbPt for run "<< Epoch[i][j] << endl;
	thisTpol.ReadInThisRun( Epoch[i][j], Target, Period );
	thisTpol.SetDFs( DF_Data[i] ); // This assumes, of course, that the epochs in Epoch and DF_Data line up, so be careful of that...
	thisTpol.CalculateAllRaw();
	thisTpol.CalculatePbPt();
	double bt = thisTpol.getBT_Pol(); double btErr = thisTpol.getBT_Pol_Err();
	if( bt != 0.0 && btErr > 0.0 ){
	    runs.push_back( Epoch[i][j] ); pbpt.push_back( bt );
	    pbptErr.push_back( btErr ); zeros.push_back(0.0);
	    //NMRpbpt.push_back( 0.8 * Period.getTargetPolarization( Epoch[i][j] ) );
	    NMRpbpt.push_back( 0.8 * Period.getOfflineTPol( Epoch[i][j] ));
	    NMRpbptErr.push_back( 0.8 * Period.getOfflineTPolErr( Epoch[i][j] ));
	}
     }
    }
    if( runs.size() > 0 ){
	//string name = to_string(runs[i]) + "-" + to_string(runs[ runs.size()-1 ]);
	TGraphErrors* PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), pbpt.data(), zeros.data(), pbptErr.data() );
	PbPt_Vals->SetMarkerColor(kBlue); PbPt_Vals->SetMarkerStyle(kFullCircle);
	//TGraphErrors* NMR_PbPt_Vals = new TGraphErrors( runs.size(), runs.data(), NMRpbpt.data(), zeros.data(), NMRpbptErr.data() );
	//NMR_PbPt_Vals->SetMarkerColor(kViolet); NMR_PbPt_Vals->SetMarkerStyle(kFullSquare);
	TMultiGraph* mg = new TMultiGraph();
	mg->Add( PbPt_Vals, "p");
	//mg->Add( NMR_PbPt_Vals, "p");
	
	return mg;
    }
    else{
	cout <<"ERROR: Couldn't calculate any PbPt for this run range. Returning empty TMultiGraph...\n";
	TMultiGraph* empty = new TMultiGraph();
	return empty;
    }

}


// This function takes in a DataSet object and makes plots for Q2 bins
TMultiGraph* Make_NH3_Bin_Plots( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNH3 = new TMultiGraph();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFNH3 = AllData.getDF_NH3(qmid, xmid);
	    double thisErrDFNH3 = AllData.getErrDF_NH3(qmid, xmid);
	    if( thisDFNH3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisDFNH3 );
		dfnh3errs.push_back( thisErrDFNH3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgNH3->Add( gr, "p" );
	    pint++;
	}
    }

    string title = "DF_{NH3} for "+sector+"; X; DF_{NH3}";
    mgNH3->SetTitle(title.c_str());
    mgNH3->GetXaxis()->SetLimits(0,0.8); mgNH3->GetYaxis()->SetRangeUser(0.0,0.4);

    return mgNH3;
}

TCanvas* NH3_Legend_Plot( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNH3 = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFNH3 = AllData.getDF_NH3(qmid, xmid);
	    double thisErrDFNH3 = AllData.getErrDF_NH3(qmid, xmid);
	    if( thisDFNH3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisDFNH3 );
		dfnh3errs.push_back( thisErrDFNH3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgNH3->Add( gr, "p" ); 
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "DF_{NH3} for "+sector+"; X; DF_{NH3}";
    mgNH3->SetTitle(title.c_str());
    mgNH3->GetXaxis()->SetLimits(0,0.8); mgNH3->GetYaxis()->SetRangeUser(0.0,0.4);

    string pltTitle = "NH3_Plot_"+sector;
    TCanvas* NH3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    NH3_Plot->cd();
    mgNH3->Draw("AP");
    mgLeg->Draw("same");
    return NH3_Plot;
}

TCanvas* A1_Legend_Plot( PbPt& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNH3 = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, a1nh3vals, zeros, a1nh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisA1NH3 = AllData.getA1(qmid, xmid);
	    double thisErrA1NH3 = AllData.getErrA1(qmid, xmid);
	    if( thisA1NH3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		a1nh3vals.push_back( thisA1NH3 );
		a1nh3errs.push_back( thisErrA1NH3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), a1nh3vals.data(), zeros.data(), a1nh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgNH3->Add( gr, "p" ); 
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "A1_{P} for "+sector+"; X; A1_{P}";
    mgNH3->SetTitle(title.c_str());
    mgNH3->GetXaxis()->SetLimits(0,0.8); mgNH3->GetYaxis()->SetRangeUser(0.0,0.4);

    string pltTitle = "A1p_Plot_"+sector;
    TCanvas* NH3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    NH3_Plot->cd();
    mgNH3->Draw("AP");
    mgLeg->Draw("same");
    return NH3_Plot;
}

// Makes a plot of the numerator term from the DF through PF calculation
TCanvas* DF_Numerator_Legend_Plot( DataSet& AllData, string sector, string targetType ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    //if( thruPF == "YES" ){
	AllData.CalculateDFThruPF( targetType );
	cout << "Calculating DF thru PF for " << targetType << endl;
    //}
    //else
	//AllData.CalculateDF();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfvals, zeros, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFnumerator = 0; double thisErrDFnumerator = 0;
	    
	    thisDFnumerator = AllData.getNumeratorDFThruPF( qmid, xmid, targetType );
	    thisErrDFnumerator = AllData.getNumeratorDFThruPFError( qmid, xmid, targetType );
	    
            if( thisDFnumerator != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    dfvals.push_back( thisDFnumerator );
		    dferrs.push_back( thisErrDFnumerator );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), zeros.data(), dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" DF_{"+targetType+"} 'Numerator' Term for "+sector+"; X; Numerator Term";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); //mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+"_DF_Num_Plots_"+sector;
    TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    DF_Plot->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    return DF_Plot;
}

// Makes a plot of the numerator term from the DF through PF calculation
TCanvas* DF_NA_Counts_Legend_Plot( DataSet& AllData, string sector, string targetType ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    //if( thruPF == "YES" ){
	AllData.CalculateDFThruPF( targetType );
	cout << "Calculating DF thru PF for " << targetType << endl;
    //}
    //else
	//AllData.CalculateDF();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfvals, zeros, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDFnacounts = 0; double thisErrDFnacounts = 0;
	    
	    thisDFnacounts = AllData.getNACountsDFThruPF( qmid, xmid, targetType );
	    thisErrDFnacounts = AllData.getNACountsDFThruPFError( qmid, xmid, targetType );
	    
            if( thisDFnacounts != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    dfvals.push_back( thisDFnacounts );
		    dferrs.push_back( thisErrDFnacounts );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), zeros.data(), dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" DF_{"+targetType+"} 'n_{A} Counts' Term for "+sector+"; X; n_{A} Counts Term";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); //mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+"_DF_NA_Counts_Plots_"+sector;
    TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    DF_Plot->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    return DF_Plot;
}

//TCanvas* Counts_Q2_Legend_Plot( DataSet& AllData, string sector, string numTarg, string denomTarg ){
TMultiGraph* Counts_Q2_Legend_Plot( DataSet& AllData, string sector, string numTarg, string denomTarg ){

    if( !( numTarg == "NH3" || numTarg == "ND3" || numTarg == "C" || numTarg == "CH2" || numTarg == "CD2" || numTarg == "F" || numTarg == "ET" || 
	   denomTarg == "NH3" || denomTarg == "ND3" || denomTarg == "C" || denomTarg == "CH2" || denomTarg == "CD2" || denomTarg == "F" || denomTarg == "ET" ) ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    //TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    TMultiGraph* emptyC;
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgRatio = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, ratiovals, zeros, ratioerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisRatio = 0; double thisErrRatio = 0;
	    thisRatio = AllData.getCountRatio( qmid, xmid, numTarg, denomTarg );
	    thisErrRatio = AllData.getCountRatioErr( qmid, xmid, numTarg, denomTarg );
            if( thisRatio != 0.0 && thisErrRatio){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    ratiovals.push_back( thisRatio );
		    ratioerrs.push_back( thisErrRatio );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), ratiovals.data(), zeros.data(), ratioerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgRatio->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "Ratio of n_{"+numTarg+"}/n_{"+denomTarg+"}; X; n_{"+numTarg+"}/n_{"+denomTarg+"}";
    mgRatio->SetTitle(title.c_str());
    //mgRatio->GetXaxis()->SetLimits(0,0.8); mgRatio->GetYaxis()->SetRangeUser(0.0,0.5);
   /* 
    string pltTitle = numTarg+"_to_"+denomTarg+"_Ratio_Plots_"+sector;
    TCanvas* Ratio_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    Ratio_Plot->cd();
    mgRatio->Draw("AP");
    mgLeg->Draw("same");
   */
    //return Ratio_Plot;
    return mgRatio;
}

// This makes a plot of the count ratios for all target types.
TCanvas* Counts_Q2_Legend_All_Plots( DataSet& AllData, string Runperiod, string denomTarg ){

    if( !( denomTarg == "NH3" || denomTarg == "ND3" || denomTarg == "C" || denomTarg == "CH2" || denomTarg == "CD2" || denomTarg == "F" || denomTarg == "ET" ) ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    //TMultiGraph* emptyC;
	    return emptyC;
    }

    string title = Runperiod+"_Data";
    TCanvas* All_Plots = new TCanvas(title.c_str(),title.c_str(),1400,700);
    All_Plots->Divide(3,2);
/*
    TCanvas* NH3_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "NH3", denomTarg );
    All_Plots->cd(1); NH3_Ratio->Draw();
    TCanvas* ND3_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "ND3", denomTarg );
    All_Plots->cd(2); ND3_Ratio->Draw();
    TCanvas* C_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "C", denomTarg );
    All_Plots->cd(3); C_Ratio->Draw();
    TCanvas* CH2_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "CH2", denomTarg );
    All_Plots->cd(4); CH2_Ratio->Draw();
    TCanvas* CD2_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "CD2", denomTarg );
    All_Plots->cd(5); CD2_Ratio->Draw();
    TCanvas* ET_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "ET", denomTarg );
    All_Plots->cd(6); ET_Ratio->Draw();
*/
    vector<string> AllTargs = {"NH3","ND3","C","CH2","CD2","ET","F"};
    int it = 1; // iterator used for plotting
    for( size_t i=0; i<AllTargs.size(); i++ ){
	string numTarg = AllTargs[i];
	if( numTarg != denomTarg ){
	    TMultiGraph* thisRatio = Counts_Q2_Legend_Plot( AllData, Runperiod, numTarg, denomTarg );
	    All_Plots->cd(it); thisRatio->Draw("AP");
	    it++;
	}
    }
/*
    TMultiGraph* NH3_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "NH3", denomTarg );
    All_Plots->cd(1); NH3_Ratio->Draw("AP");
    TMultiGraph* ND3_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "ND3", denomTarg );
    All_Plots->cd(2); ND3_Ratio->Draw("AP");
    TMultiGraph* C_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "C", denomTarg );
    All_Plots->cd(3); C_Ratio->Draw("AP");
    TMultiGraph* CH2_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "CH2", denomTarg );
    All_Plots->cd(4); CH2_Ratio->Draw("AP");
    TMultiGraph* CD2_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "CD2", denomTarg );
    All_Plots->cd(5); CD2_Ratio->Draw("AP");
    TMultiGraph* ET_Ratio = Counts_Q2_Legend_Plot( AllData, sector, "ET", denomTarg );
    All_Plots->cd(6); ET_Ratio->Draw("AP");
*/
    return All_Plots;
}

TCanvas* DF_Legend_Plot( DataSet& AllData, string sector, string targetType, string thruPF ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgDF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    if( thruPF == "YES" ){
	AllData.CalculateDFThruPF( targetType );
	cout << "Calculating DF thru PF for " << targetType << endl;
    }
    else
	AllData.CalculateDF();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfvals, zeros, dferrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisDF = 0; double thisErrDF = 0;
	    if( targetType == "NH3" && thruPF != "YES" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_NH3(qmid, xmid);
	    }
	    else if( targetType == "NH3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF != "YES" ){
	        thisDF = AllData.getDF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_ND3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_ND3(qmid, xmid);
	    }
            if( thisDF != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    dfvals.push_back( thisDF );
		    dferrs.push_back( thisErrDF );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfvals.data(), zeros.data(), dferrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgDF->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" DF_{"+targetType+"} for "+sector+"; X; DF_{"+targetType+"}";
    mgDF->SetTitle(title.c_str());
    mgDF->GetXaxis()->SetLimits(0,0.8); 
    if( targetType == "NH3" ) mgDF->GetYaxis()->SetRangeUser(0.1,0.3);
    else mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+"_DF_Plots_"+sector;
    TCanvas* DF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    DF_Plot->cd();
    mgDF->Draw("AP");
    mgLeg->Draw("same");

    return DF_Plot;
}

TCanvas* PF_Legend_Plot( DataSet& AllData, string sector, string targetType ){

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas("emptyC","emptyC",800,600);
	    return emptyC;
    }
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgPFbath = new TMultiGraph();
    TMultiGraph* mgPFcell = new TMultiGraph();
    TLegend* mgLegbath = new TLegend(0.1,0.7,0.4,0.9);
    TLegend* mgLegcell = new TLegend(0.1,0.7,0.4,0.9);
    mgLegbath->SetNColumns(3);
    mgLegbath->SetHeader("Q^{2} Bins (GeV^{2})","C");
    mgLegcell->SetNColumns(3);
    mgLegcell->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculatePF();

    // These are used to calculate the weighted average of all the PF values in the data set
    double totalPFbath = 0; double denominatorbath = 0;
    double totalPFcell = 0; double denominatorcell = 0;

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbinsbath, xbinscell, pfcellvals, cellzeros, pfcellerrs, bathzeros, pfbathvals, pfbatherrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisPFbath = 0; double thisErrPFbath = 0;
	    double thisPFcell = 0; double thisErrPFcell = 0;
	    if( targetType == "NH3" ){
	        thisPFbath = AllData.getPF_bath_NH3(qmid, xmid);
	        thisErrPFbath = AllData.getErrPF_bath_NH3(qmid, xmid);
	        thisPFcell = AllData.getPF_cell_NH3(qmid, xmid);
	        thisErrPFcell = AllData.getErrPF_cell_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" ){
	        thisPFbath = AllData.getPF_bath_ND3(qmid, xmid);
	        thisErrPFbath = AllData.getErrPF_bath_ND3(qmid, xmid);
	        thisPFcell = AllData.getPF_cell_ND3(qmid, xmid);
	        thisErrPFcell = AllData.getErrPF_cell_ND3(qmid, xmid);
	    }
            if( thisPFbath != 0.0 ){
		    xbinsbath.push_back(xmid); bathzeros.push_back(0.0);
		    pfbathvals.push_back( thisPFbath );
		    pfbatherrs.push_back( thisErrPFbath );
		    //if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPFbath += thisPFbath / pow(thisErrPFbath,2) ;
		    denominatorbath += 1.0 / pow(thisErrPFbath,2);
	    	    //}
	    }
	    if( thisPFcell != 0.0 ){
		    xbinscell.push_back(xmid); cellzeros.push_back(0.0);
		    pfcellvals.push_back( thisPFcell );
		    pfcellerrs.push_back( thisErrPFcell );
		    //if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPFcell += thisPFcell / pow(thisErrPFcell,2) ;
		    denominatorcell += 1.0 / pow(thisErrPFcell,2);
	    	    //}
	    }	    
	}
	if( xbinsbath.size() > 0 && xbinscell.size() > 0 ){
	    TGraphErrors* gr_cell = new TGraphErrors(xbinscell.size(), xbinscell.data(), pfcellvals.data(), cellzeros.data(), pfcellerrs.data() );
	    TGraphErrors* gr_bath = new TGraphErrors(xbinsbath.size(), xbinsbath.data(), pfbathvals.data(), bathzeros.data(), pfbatherrs.data() );
	    gr_cell->SetMarkerColor( palette[pint] );
	    gr_bath->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr_cell->SetMarkerStyle(kFullCircle); gr_bath->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr_cell->SetMarkerStyle(kFullTriangleUp); gr_bath->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr_cell->SetMarkerStyle(kFullSquare); gr_bath->SetMarkerStyle(kFullSquare); }
	    mgPFcell->Add( gr_cell, "p" );
	    mgPFbath->Add( gr_bath, "p" );
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLegcell->AddEntry( gr_cell, legTitle.str().c_str(), "p");
	    mgLegbath->AddEntry( gr_bath, legTitle.str().c_str(), "p");

	    pint++;
	}
    }
    if( denominatorcell != 0.0 && denominatorbath != 0.0 ){
        double avgPFcell = totalPFcell / denominatorcell;
        cout << "Average PF_cell for data set is: "<< avgPFcell <<" +/- "<< 1.0/sqrt(denominatorbath) << endl;
        double avgPFbath = totalPFbath / denominatorbath;
        cout << "Average PF_bath for data set is: "<< avgPFbath <<" +/- "<< 1.0/sqrt(denominatorbath) << endl;
    }

    string title = targetType+" PF_{bath} for "+sector+"; X; "+targetType+" PF_{bath}";
    mgPFbath->SetTitle(title.c_str());
    mgPFbath->GetXaxis()->SetLimits(0,0.8); mgPFbath->GetYaxis()->SetRangeUser(0.4,0.7);
    
    title = targetType+" PF_{cell} for "+sector+"; X; "+targetType+" PF_{cell}";
    mgPFcell->SetTitle(title.c_str());
    mgPFcell->GetXaxis()->SetLimits(0,0.8); mgPFcell->GetYaxis()->SetRangeUser(0.4,0.7);

    string pfbathTitle = targetType+"_PF_bath_avg_"+sector;
    TF1* pfbathavg = new TF1(pfbathTitle.c_str(),"[0]",-1,2);
    pfbathavg->SetLineColor(kBlack);
    if( denominatorbath != 0 ) pfbathavg->SetParameter( 0, totalPFbath/denominatorbath );
    else pfbathavg->SetParameter( 0, totalPFbath/denominatorbath );

    string pltTitle = targetType+"_PF_Plots_"+sector;
    //TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),1400,600);
    TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    //PF_Plot->Divide(2,1);
    //PF_Plot->cd(1);
    PF_Plot->cd();
    mgPFbath->Draw("AP");
    mgLegbath->Draw("same");
    pfbathavg->Draw("same");
/*    PF_Plot->cd(2);
    mgPFcell->Draw("AP");
    mgLegcell->Draw("same");
*/
    return PF_Plot;
}

// Plots the raw double-spin asymmetries for NH3 targets
TCanvas* AllRaw_Legend_Plot( DataSet& AllData, string sector, string targetType, string FCcorr ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    if( !( targetType == "NH3" || targetType == "ND3" ) ){
	cout << "Valid target types are 'NH3' or 'ND3'\n";
	TCanvas* empty = new TCanvas();
	return empty;
    }

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = 0.0; double thisErrAllRaw = 0.0;
	    if( FCcorr == "Yes" ){
	        thisAllRaw = AllData.getAllRawNoFC(qmid, xmid, targetType ); 
	        thisErrAllRaw = AllData.getAllRawNoFCErr(qmid, xmid, targetType);//0.05*thisAllRaw; // Set this to 5% error as a placeholder
	    }
	    else{
	        thisAllRaw = AllData.getAllRaw(qmid, xmid, targetType ); 
	        thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, targetType);//0.05*thisAllRaw; // Set this to 5% error as a placeholder
	    }
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = targetType+" A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRaw"+targetType+" Plots "+sector;
    TCanvas* AllRaw_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRaw_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRaw_Plot;
}

// Plots the raw double-spin asymmetries for NH3 targets
TCanvas* AllRawNH3_Legend_Plot( DataSet& AllData, string sector ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(qmid, xmid, "NH3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, "NH3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "NH3 A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRawNH3 Plots "+sector;
    TCanvas* AllRawNH3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRawNH3_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRawNH3_Plot;
}
// Plots the raw double-spin asymmetries for ND3 targets
TCanvas* AllRawND3_Legend_Plot( DataSet& AllData, string sector ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgAllRaw = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for ND3

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(qmid, xmid, "ND3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(qmid, xmid, "ND3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    xbins.push_back(xmid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), arawvals.data(), zeros.data(), arawerrs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ){ gr->SetMarkerStyle(kFullCircle); }
	    else if(pint >= 4 && pint < 8){ gr->SetMarkerStyle(kFullTriangleUp); }
	    else if(pint >= 8){ gr->SetMarkerStyle(kFullSquare); }
	    mgAllRaw->Add( gr, "p" );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }

    string title = "ND3 A_{||,raw} for "+sector+"; X; A_{||}";
    mgAllRaw->SetTitle(title.c_str());
    mgAllRaw->GetXaxis()->SetLimits(0,0.8); mgAllRaw->GetYaxis()->SetRangeUser(-0.07,0.07);
    
    string pltTitle = "AllRawND3 Plots "+sector;
    TCanvas* AllRawND3_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    AllRawND3_Plot->cd();
    //gPad->SetLogx();
    mgAllRaw->Draw("AP");
    mgLeg->Draw("same");

    return AllRawND3_Plot;
}

// Plots the raw double-spin asymmetries for NH3 targets
vector<TCanvas*> AllRawNH3_Q2_Bin_Plot( DataSet& AllData, string sector ){
    
    string title1 = sector + " 1";
    string title2 = sector + " 2";
    TCanvas* q2BinPlots1 = new TCanvas(title1.c_str(),title1.c_str(),1900,1000); // Plots the first 15 x bins
    TCanvas* q2BinPlots2 = new TCanvas(title2.c_str(),title2.c_str(),1900,1000); // Plots the last 15 x bins
    q2BinPlots1->Divide(5,3,0,0); q2BinPlots2->Divide(5,3,0,0);
    // Iterators used for drawing to the correct canvas
    int plotIt = 1;

    AllData.CalculateAllRaw(); // Calculate the raw double-spin asymmetries for NH3

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, arawvals, zeros, arawerrs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double q2mid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisAllRaw = AllData.getAllRaw(q2mid, xmid, "NH3" ); 
	    double thisErrAllRaw = AllData.getAllRawErr(q2mid, xmid, "NH3");//0.05*thisAllRaw; // Set this to 5% error as a placeholder
            if( thisAllRaw != 0.0 ){
		    q2bins.push_back(q2mid); zeros.push_back(0.0);
		    arawvals.push_back( thisAllRaw );
		    arawerrs.push_back( thisErrAllRaw );
	    }
	}
	if( q2bins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(q2bins.size(), q2bins.data(), arawvals.data(), zeros.data(), arawerrs.data() ); gr->SetTitle("");
	    gr->SetMarkerStyle(kFullCircle);
	    gr->GetXaxis()->SetLimits(0,11); // 0 to 11 in GeV^2
	    gr->GetYaxis()->SetRangeUser(0,0.2);// Raw asymmetry from 0 to 5%
	    string title = "A_{||} at X = "+to_string(xmid);
	    TLegend* leg = new TLegend(0.1,0.5,0.7,0.9);
	    leg->AddEntry(gr, title.c_str(), ""); leg->SetFillStyle(0); leg->SetBorderSize(0);
	    if( plotIt > 15 ){
		q2BinPlots2->cd(plotIt-15); //gr->SetTitle(title.c_str());
		//gPad->SetLogx();
		gr->Draw("AP"); leg->Draw("same"); plotIt++;
	    }		
	    else{
		q2BinPlots1->cd(plotIt); //gr->SetTitle(title.c_str());
		//gPad->SetLogx();
		gr->Draw("AP"); leg->Draw("same"); plotIt++;
	    }
	}
    }

    vector<TCanvas*> Plots = { q2BinPlots1, q2BinPlots2 };
    return Plots;
}

// Makes a plot of the PF for the target cell only, excluding the LHe bath
TCanvas* PFcell_Legend_Plot( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgPF = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Q^{2} Bins (GeV^{2})","C");
    AllData.CalculatePF();

    // These are used to calculate the weighted average of all the PF values in the data set
    double totalPF = 0; double denominator = 0;

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    // Multiplying by LHe / 5 cm gives the PF for the cell alone
	    double thisPF = (LHe/5.0)*AllData.getPF_bath_NH3(qmid, xmid);
	    double thisErrPF = (LHe/5.0)*AllData.getErrPF_bath_NH3(qmid, xmid);
	    if( thisPF != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisPF );
		dfnh3errs.push_back( thisErrPF );
		if( thisPF > 0.0 && thisErrPF != 0.0 && abs(thisErrPF) < (0.1*thisPF) ){ 
		    totalPF += thisPF / pow(thisErrPF,2) ;
		    denominator += 1.0 / pow(thisErrPF,2);
	    	}
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerColor( palette[pint] );
	    if(pint < 4 ) gr->SetMarkerStyle(kFullCircle); 
	    else if(pint >= 4 && pint < 8) gr->SetMarkerStyle(kFullTriangleUp);
	    else if(pint >= 8) gr->SetMarkerStyle(kFullSquare);
	    mgPF->Add( gr, "p" ); 
	    //string legTitle = "Q^{2}=" + to_string( trunc(qmid*1000)/1000 );
	    stringstream legTitle; legTitle << "Q^{2}=" << fixed << setprecision(3) << qmid;
	    mgLeg->AddEntry( gr, legTitle.str().c_str(), "p");
	    pint++;
	}
    }
    if( denominator != 0.0 ){
        double avgPF = totalPF / denominator;
        cout << "Average PF_cell for data set is: "<< avgPF <<" +/- "<< 1.0/sqrt(denominator) << endl;
    }
    string title = "NH3 PF_{cell} for "+sector+"; X; NH3 PF_{cell}";
    mgPF->SetTitle(title.c_str());
    mgPF->GetXaxis()->SetLimits(0,0.8); mgPF->GetYaxis()->SetRangeUser(0.0,1.0);

    string pltTitle = "PF_{cell} Plot "+sector;
    TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),800,600);
    PF_Plot->cd();
    mgPF->Draw("AP");
    mgLeg->Draw("same");
    return PF_Plot;
}

// This function creates a single plot for a given epoch, creating DF values that are averaged over
// all Q^2 bins for a fixed bin in Bjorken x.
TCanvas* DF_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, string thruPF ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas();
	    return emptyC;
    }

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    /*
    if( thruPF == "YES" )
	Epochs.CalculateDFThruPF( targetType );
    else
	Epochs.CalculateDF();
    */
    //vector<double> xbins, zeros, avg_dfs, avg_df_errs;
    /*
    vector<double> modXbin, modZeros, modAvgdfs; // Model values for comparison
    string toFile = THIS_DIR + "Darren_Model_DF.txt";
    // Loop over the model data first
    DataSet ModelDFs;
    ifstream fin(toFile.c_str()); string line;
    while(getline(fin,line)){
	stringstream sin(line);
	double q2, x, qbin, xbin, df, dferr;
	sin >> q2 >> x >> qbin >> xbin >> df >> dferr;
        ModelDFs.SetDFs( q2+0.1, x+0.01, df, 0 );
    }
    fin.close();
    */
    // Now read in the data
 TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
 leg->SetHeader("Epochs","C"); leg->SetNColumns(2);
 TMultiGraph* mg = new TMultiGraph();
 mg->SetTitle(epoch.c_str());
 vector<TGraphErrors*> DF_Epoch_Plots;
 for(size_t k=0; k<Epochs.size(); k++){
    if( thruPF == "YES" )
	Epochs[k].CalculateDFThruPF( targetType );
    else
	Epochs[k].CalculateDF();

    DataSet AllData = Epochs[k]; // Very inefficient!!! Needs improvement!!!!11!
    vector<double> xbins, zeros, avg_dfs, avg_df_errs;

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, dfvals, dferrs, modelq2bins, modeldfvals, modeldferrs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisDF = 0; double thisErrDF = 0;
	    // double modelDF = ModelDFs.getDF_FixedPF_NH3( qmid, xmid );
	    // double modelDFerr = 0;
	    if( targetType == "NH3" && thruPF != "YES" ){
	        thisDF = AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_NH3(qmid, xmid);
	    }
	    else if( targetType == "NH3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF != "YES" ){
	        thisDF = AllData.getDF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_ND3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_ND3(qmid, xmid);
	    }
            if( thisDF != 0.0 ){
		    q2bins.push_back(qmid);
		    dfvals.push_back( thisDF );
		    dferrs.push_back( thisErrDF );
	    }
/*
	    if( modelDF != 0.0 && k==0 ){
		    modelq2bins.push_back(qmid);
		    modeldfvals.push_back( modelDF );
		    modeldferrs.push_back( modelDFerr );	
	    }
*/
	}
	if( q2bins.size() > 0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<q2bins.size(); k++){
		num += dfvals[k] / pow(dferrs[k],2);
		den += 1.0 / pow(dferrs[k],2);
	    }
	    if( den != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		avg_dfs.push_back( num / den ); avg_df_errs.push_back( 1.0 / sqrt(den) );
	    }
	}
	/*
	if( modelq2bins.size() > 0 && k==0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<modelq2bins.size(); k++){
		num += modeldfvals[k] ;/// pow(modeldferrs[k],2);
		den += 1.0 ;/// pow(modeldferrs[k],2);
	    }
	    if( den != 0.0 ){
		modXbin.push_back(xmid); modAvgdfs.push_back( num/den ); modZeros.push_back(0.0);
	    }
	}
	*/
    }
    //cout << "Epoch "<< epoch <<" first bin value: "<< avg_dfs[0] << endl;
    TGraphErrors* thisGR = new TGraphErrors(xbins.size(), xbins.data(), avg_dfs.data(), zeros.data(), avg_df_errs.data() );
    DF_Epoch_Plots.push_back(thisGR);
/*
    string thisTitle = "Epoch "+to_string(k+2); thisGR->SetTitle(thisTitle.c_str());
    thisGR->SetMarkerColor( palette[pint+1] ); thisGR->SetMarkerStyle(kFullCircle); pint++;
    mg->Add( thisGR , "p");
    leg->AddEntry( thisGR , thisTitle.c_str(), "p");
*/
 }
 for( size_t i=0; i<DF_Epoch_Plots.size(); i++){
/*
    DF_Epoch_Plots[0]->SetMarkerColor(kRed);   DF_Epoch_Plots[0]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[1]->SetMarkerColor(kOrange);DF_Epoch_Plots[1]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[2]->SetMarkerColor(kGreen); DF_Epoch_Plots[2]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[3]->SetMarkerColor(kBlue);  DF_Epoch_Plots[3]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[4]->SetMarkerColor(kViolet);DF_Epoch_Plots[4]->SetMarkerStyle(kFullTriangleDown);
    DF_Epoch_Plots[5]->SetMarkerColor(kTeal);  DF_Epoch_Plots[5]->SetMarkerStyle(kFullTriangleDown);
*/
    DF_Epoch_Plots[i]->SetMarkerColor(palette[i]); 
    if( i<2 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullCircle);
    else if( i<4 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleUp);
    else if( i<6 ) DF_Epoch_Plots[i]->SetMarkerStyle(kFullTriangleDown);
    else DF_Epoch_Plots[i]->SetMarkerStyle(kFullSquare);
 }

    for(size_t i=0; i<DF_Epoch_Plots.size(); i++){
	string thisTitle;
	if( epoch == "Summer" ) thisTitle = "Epoch "+to_string(i+2);
	else if( epoch == "Fall" ) thisTitle = "Epoch "+to_string(i+8);
	else if( epoch == "Spring" ) thisTitle = "Epoch "+to_string(i+15);
        //DF_Epoch_Plots[i]->SetMarkerColor( palette[i+1] );
        //DF_Epoch_Plots[i]->SetMarkerStyle(i);
        mg->Add( DF_Epoch_Plots[i] , "p");
        leg->AddEntry( DF_Epoch_Plots[i] , thisTitle.c_str(), "p");
    }
/*
    TGraphErrors* modDF = new TGraphErrors(modXbin.size(), modXbin.data(), modAvgdfs.data(), modZeros.data(), modZeros.data() );
    modDF->SetMarkerColor(kBlack); modDF->SetMarkerStyle(kFullSquare);
    leg->AddEntry( modDF, "Model", "p");
    mg->Add( modDF, "p");
*/
    mg->GetYaxis()->SetRangeUser(0,0.4);
    mg->GetXaxis()->SetRangeUser(0,1);
    TCanvas* mgplt = new TCanvas(epoch.c_str(),epoch.c_str(),800,600);
    mgplt->cd();
    string mgTitle;
    if( epoch == "Summer" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; D_{F}";
    else if( epoch == "Fall" ) mgTitle = epoch+" 2022 "+ targetType+" Data; X; D_{F}";
    else if( epoch == "Spring" ) mgTitle = epoch+" 2023 "+ targetType+" Data; X; D_{F}";
    mg->SetTitle(mgTitle.c_str()); 
    mg->Draw("ap"); leg->Draw("same");

    return mgplt;

}


TCanvas* AllRaw_Q2_Averaged_Plot( vector<DataSet>& Epochs, string epoch, string targetType, string thruPF ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    if( !(targetType == "NH3" || targetType == "ND3") ){
	    cout <<"ERROR: enter a valid target type to plot: 'NH3' or 'ND3'\n";
	    TCanvas* emptyC = new TCanvas();
	    return emptyC;
    }

    // Check to see if we're using the average PF across all bins (thruPF = "YES") or if the
    // PF is being calculated for each bin separately (thruPF = else)
    /*
    if( thruPF == "YES" )
	Epochs.CalculateDFThruPF( targetType );
    else
	Epochs.CalculateDF();
    */
    //vector<double> xbins, zeros, avg_dfs, avg_df_errs;
    vector<double> modXbin, modZeros, modAvgdfs; // Model values for comparison
    string toFile = THIS_DIR + "Darren_Model_DF.txt";
    // Loop over the model data first
    DataSet ModelDFs;
    ifstream fin(toFile.c_str()); string line;
    while(getline(fin,line)){
	stringstream sin(line);
	double q2, x, qbin, xbin, df, dferr;
	sin >> q2 >> x >> qbin >> xbin >> df >> dferr;
        ModelDFs.SetDFs( q2+0.1, x+0.01, df, 0 );
    }
    fin.close();
    // Now read in the data
 TLegend* leg = new TLegend(0.1,0.7,0.3,0.9);
 leg->SetHeader("Epochs","C"); leg->SetNColumns(2);
 TMultiGraph* mg = new TMultiGraph();
 mg->SetTitle(epoch.c_str());
 vector<TGraphErrors*> DF_Epoch_Plots;
 for(size_t k=0; k<Epochs.size(); k++){
    Epochs[k].CalculateAllRaw();
    if( thruPF == "YES" ){
	Epochs[k].CalculateDFThruPF( targetType );
    }
    else{
	Epochs[k].CalculateDF();
    }

    DataSet AllData = Epochs[k]; // Very inefficient!!! Needs improvement!!!!11!
    vector<double> xbins, zeros, avg_dfs, avg_df_errs;

    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+X_Bin_Bounds[i+1])/2.0;
	vector<double> q2bins, dfvals, dferrs, modelq2bins, modeldfvals, modeldferrs;
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    double thisDF = 0; double thisErrDF = 0;
	    double modelDF = ModelDFs.getDF_FixedPF_NH3( qmid, xmid );
	    double modelDFerr = 0;
	    //if( targetType == "NH3" && thruPF != "YES" ){
	        thisDF = AllData.getAllRaw(qmid, xmid, targetType);//AllData.getDF_NH3(qmid, xmid);
	        thisErrDF = AllData.getAllRawErr(qmid, xmid, targetType);//getErrDF_NH3(qmid, xmid);
	    //}
/*
	    else if( targetType == "NH3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_NH3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_NH3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF != "YES" ){
	        thisDF = AllData.getDF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_ND3(qmid, xmid);
	    }
	    else if( targetType == "ND3" && thruPF == "YES" ){
	        thisDF = AllData.getDF_FixedPF_ND3(qmid, xmid);
	        thisErrDF = AllData.getErrDF_FixedPF_ND3(qmid, xmid);
	    }
*/
            if( thisDF != 0.0 ){
		    q2bins.push_back(qmid);
		    dfvals.push_back( thisDF );
		    dferrs.push_back( thisErrDF );
	    }
	    if( modelDF != 0.0 && k==0 ){
		    modelq2bins.push_back(qmid);
		    modeldfvals.push_back( modelDF );
		    modeldferrs.push_back( modelDFerr );	
	    }
	}
	if( q2bins.size() > 0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<q2bins.size(); k++){
		num += dfvals[k] / pow(dferrs[k],2);
		den += 1.0 / pow(dferrs[k],2);
	    }
	    if( den != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		avg_dfs.push_back( num / den ); avg_df_errs.push_back( 1.0 / sqrt(den) );
	    }
	}
	if( modelq2bins.size() > 0 && k==0 ){
	    // Calculate the weighted average of all the DFs in the given xbin
	    double num = 0.0; double den = 0.0;
	    for(size_t k=0; k<modelq2bins.size(); k++){
		num += modeldfvals[k] ;/// pow(modeldferrs[k],2);
		den += 1.0 ;/// pow(modeldferrs[k],2);
	    }
	    if( den != 0.0 ){
		modXbin.push_back(xmid); modAvgdfs.push_back( num/den ); modZeros.push_back(0.0);
	    }
	}

    }
    //cout << "Epoch "<< epoch <<" first bin value: "<< avg_dfs[0] << endl;
    TGraphErrors* thisGR = new TGraphErrors(xbins.size(), xbins.data(), avg_dfs.data(), zeros.data(), avg_df_errs.data() );
    DF_Epoch_Plots.push_back(thisGR);
/*
    string thisTitle = "Epoch "+to_string(k+2); thisGR->SetTitle(thisTitle.c_str());
    thisGR->SetMarkerColor( palette[pint+1] ); thisGR->SetMarkerStyle(kFullCircle); pint++;
    mg->Add( thisGR , "p");
    leg->AddEntry( thisGR , thisTitle.c_str(), "p");
*/
 }
    DF_Epoch_Plots[0]->SetMarkerColor(kRed);   DF_Epoch_Plots[0]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[1]->SetMarkerColor(kOrange);DF_Epoch_Plots[1]->SetMarkerStyle(kFullCircle);
    DF_Epoch_Plots[2]->SetMarkerColor(kGreen); DF_Epoch_Plots[2]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[3]->SetMarkerColor(kBlue);  DF_Epoch_Plots[3]->SetMarkerStyle(kFullTriangleUp);
    DF_Epoch_Plots[4]->SetMarkerColor(kViolet);DF_Epoch_Plots[4]->SetMarkerStyle(kFullTriangleDown);
    DF_Epoch_Plots[5]->SetMarkerColor(kTeal);  DF_Epoch_Plots[5]->SetMarkerStyle(kFullTriangleDown);
    DF_Epoch_Plots[6]->SetMarkerColor(kBlack);  DF_Epoch_Plots[6]->SetMarkerStyle(kFullSquare);


    for(size_t i=0; i<DF_Epoch_Plots.size()-1; i++){
        string thisTitle = "Epoch "+to_string(i+2); //DF_Epoch_Plots[i]->SetTitle(thisTitle.c_str());
        //DF_Epoch_Plots[i]->SetMarkerColor( palette[i+1] );
        //DF_Epoch_Plots[i]->SetMarkerStyle(i);
        mg->Add( DF_Epoch_Plots[i] , "p");
        leg->AddEntry( DF_Epoch_Plots[i] , thisTitle.c_str(), "p");
    }
/*
    TGraphErrors* modDF = new TGraphErrors(modXbin.size(), modXbin.data(), modAvgdfs.data(), modZeros.data(), modZeros.data() );
    modDF->SetMarkerColor(kBlack); modDF->SetMarkerStyle(kFullSquare);
    leg->AddEntry( modDF, "Model", "p");
    mg->Add( modDF, "p");
*/
    leg->AddEntry( DF_Epoch_Plots[6], "All Su22 Data","p");
    mg->Add( DF_Epoch_Plots[6], "p");
    mg->GetYaxis()->SetRangeUser(0,0.4);
    mg->GetXaxis()->SetLimits(0,1);
    TCanvas* mgplt = new TCanvas("mgplt2","mgplt2",800,600);
    mgplt->cd();
    mg->Draw("ap"); leg->Draw("same");

    return mgplt;

}

// Using the outputs from the proceeding "DF_Q2_Averaged_Plot" function, this function takes
// a vector of pointers to the TGraphErrors objects produced by the previous function and
// plots them together.
TCanvas* Averaged_DF_By_Sector( vector<TGraphErrors*>& Graphs, string epoch ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mg = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Run Range Epochs","C");

    for(size_t i=0; i<Graphs.size(); i++){
        Graphs[i]->SetMarkerColor( palette[pint] );
	//Graphs[i]->GetXaxis()->SetRangeUser(0,0.8);
	if(pint < 4 ){ Graphs[i]->SetMarkerStyle(kFullCircle); }
	else if(pint >= 4 && pint < 8){ Graphs[i]->SetMarkerStyle(kFullTriangleUp); }
	else if(pint >= 8){ Graphs[i]->SetMarkerStyle(kFullSquare); }
	Graphs[i]->GetYaxis()->SetRangeUser(0,0.4);
	mg->Add( Graphs[i], "p" );
	string legTitle = Graphs[i]->GetTitle();
        mgLeg->AddEntry( Graphs[i], legTitle.c_str(), "p");
	pint++;
    }

    string title = "Average D_{F} for "+ epoch +"; X; D_{F}";
    //mg->SetTitle("Average D_{F} for 2022 Epochs; X; D_{F}");
    mg->SetTitle(title.c_str());
    
    TCanvas* AvgDF = new TCanvas("AvgDF","AvgDF",800,600);
    AvgDF->cd();
    mg->Draw("ap");
    mgLeg->Draw("same");

    return AvgDF;
}

// This function creates plots of the ratio of other target counts normalized to the NH3 target counts.
// Each target type is first normalized to the respective total gated FC charge.
TCanvas* PlotNH3NormalizedCounts( DataSet& AllData ){

    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    int pint = 0;

    TMultiGraph* mgNormCounts = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.1,0.7,0.4,0.9);
    //mgLeg->SetNColumns(3);
    mgLeg->SetHeader("Counts Normalized to NH3","C");

    // These vectors hold the total number normalized counts for each bin in Bjorken X.	
    vector<double> NH3_Xbins, ND3_Xbins, CH2_Xbins, C_Xbins, F_Xbins, ET_Xbins, CD2_Xbins, Xbins, Zeros;
    for(size_t i=0; i<X_Bin_Bounds.size()-1; i++){
	double xmid = (X_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
        // Get the normalized counts for each bin
	double n_NH3 = 0; double n_ND3 = 0; double n_CH2 = 0; double n_C = 0; double n_F = 0; double n_ET = 0; double n_CD2 = 0;
	
	for(size_t j=0; j<Q2_Bin_Bounds.size()-1; j++){
	    double qmid = (Q2_Bin_Bounds[j]+Q2_Bin_Bounds[j+1])/2.0;
	    Bin thisBin = AllData.getThisBin( qmid, xmid ); // Get the current bin and bin contents
	    n_NH3 += thisBin.getNormNt("NH3");
	    n_ND3 += thisBin.getNormNt("ND3");
	    n_CH2 += thisBin.getNormNt("CH2");
	    n_C   += thisBin.getNormNt("C");
	    n_F   += thisBin.getNormNt("F");
	    n_ET  += thisBin.getNormNt("ET");
	    n_CD2 += thisBin.getNormNt("CD2");
	}
	if( n_NH3 > 0 ){
	    NH3_Xbins.push_back( n_NH3/n_NH3 ); // Should just be one, but this is a sanity check.
	    ND3_Xbins.push_back( n_ND3/n_NH3 ); CH2_Xbins.push_back( n_CH2/n_NH3 ); C_Xbins.push_back( n_C/n_NH3 ); 
	    F_Xbins.push_back( n_F/n_NH3 ); ET_Xbins.push_back( n_ET/n_NH3 ); CD2_Xbins.push_back( n_CD2/n_NH3 );
	    Xbins.push_back( xmid ); Zeros.push_back(0.0);
	}
    }
    if( Xbins.size() > 0 ){
	TGraphErrors* gr_NH3 = new TGraphErrors(Xbins.size(), Xbins.data(), NH3_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_NH3->SetMarkerColor( palette[0] ); gr_NH3->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_NH3, "p" ); mgLeg->AddEntry( gr_NH3, "NH3", "p" ); 
	//gr_NH3->GetYaxis()->SetLogy();
 	TGraphErrors* gr_ND3 = new TGraphErrors(Xbins.size(), Xbins.data(), ND3_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_ND3->SetMarkerColor( palette[1] ); gr_ND3->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_ND3, "p" ); mgLeg->AddEntry( gr_ND3, "ND3", "p" ); 
	//gr_ND3->GetYaxis()->SetLogy();
 	TGraphErrors* gr_CH2 = new TGraphErrors(Xbins.size(), Xbins.data(), CH2_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_CH2->SetMarkerColor( palette[2] ); gr_CH2->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_CH2, "p" ); mgLeg->AddEntry( gr_CH2, "CH2", "p" ); 
	//gr_CH2->GetYaxis()->SetLogy();
 	TGraphErrors* gr_C = new TGraphErrors(Xbins.size(), Xbins.data(), C_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_C->SetMarkerColor( palette[3] ); gr_C->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_C, "p" ); mgLeg->AddEntry( gr_C, "C", "p" ); 
	//gr_C->GetYaxis()->SetLogy();
 	TGraphErrors* gr_F = new TGraphErrors(Xbins.size(), Xbins.data(), F_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_F->SetMarkerColor( palette[4] ); gr_F->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_F, "p" ); mgLeg->AddEntry( gr_F, "F", "p" ); 
	//gr_F->GetYaxis()->SetLogy();
 	TGraphErrors* gr_ET = new TGraphErrors(Xbins.size(), Xbins.data(), ET_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_ET->SetMarkerColor( palette[5] ); gr_ET->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_ET, "p" ); mgLeg->AddEntry( gr_ET, "ET", "p" ); 
	//gr_ET->GetYaxis()->SetLogy();
 	TGraphErrors* gr_CD2 = new TGraphErrors(Xbins.size(), Xbins.data(), CD2_Xbins.data(), Zeros.data(), Zeros.data() );
	gr_CD2->SetMarkerColor( palette[6] ); gr_CD2->SetMarkerStyle(kFullCircle);
	mgNormCounts->Add( gr_CD2, "p" ); mgLeg->AddEntry( gr_CD2, "CD2", "p" ); 
	//gr_CD2->GetYaxis()->SetLogy();
	TCanvas* NH3_Norm_Plots = new TCanvas("NH3_Norm_Plots","NH3_Norm_Plots",800,600);
	NH3_Norm_Plots->cd();
	gPad->SetLogy();
	mgNormCounts->Draw("ap");
	mgLeg->Draw("same");
	return NH3_Norm_Plots;
    }
    else{
	cout << "Couldn't read in any counts; check inputs and try again.\n";
	TCanvas* nullPlt = new TCanvas();
	return nullPlt;
    }
}

