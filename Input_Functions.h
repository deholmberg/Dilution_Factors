// This header file contains the information required to make the bins for
// the dilution factors in x, Q2.

using namespace std;

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

const double lC = 1.678; // cm (length of carbon target)
const double LHe = 5.86; // cm (length of LHe bath)
const double lCH = 3.18; // cm (length of CH2 target)
const double lCD = 2.686;// cm (length of CD2 target)
const double pC = 1.7926/12.0; // mol/cm^3 (carbon target density)
const double pCH = 0.9425/14.027; // mol/cm^3 (CH2 target density)
const double pCD = 1.0979/16.039; // mol/cm^3 (CD2 target density)
const double pA = 0.867/17.031; // mol/cm^3 (NH3 target density)
const double pD = 1.007/20.048; // mol/cm^3 (ND3 target density)

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
double Coeffnd3 = abs( -9.0*pA*aa/cc );

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
	double numerator = Coeffnd3*(nA-nET)*(uCHnd3*nCD - uCBnd3*nC - uFOnd3*nF + nET);
	double denominator = nA*(dCHnd3*nCD - dCBnd3*nC - dFOnd3*nF - nET);
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

// These constants are used to make the calulation of the DF_NH3 using the PF a little easier
const double qq = 3*pA / (2*pC*pCH); const double rr = pC - pCH + LHe*(pCH/lC - pC/lCH);
const double ss = -LHe*pCH/lC; const double tt = LHe*pC/lCH; const double uu = pCH - pC;

// Calculates the dilution factor using a previously calculated value of the packing fraction (PF)
double ND3DF_Thru_PF(double pf, double nA, double nCD, double nC, double nET, double nF){
  if( nA != 0.0 ) return (qq*pf/nA)*(rr*nET + ss*nC + tt*nCD + uu*nF );
  else return 0.0; // error case
}
// Calculates the error in the dilution factor calculated with the PF as input; does NOT include error
// propagation from the PF right now
double ND3DF_Thru_PF_Error(double pf, double nA, double nCD, double nC, double nET, double nF, double dpf, double dnA, double dnCD, double dnC, double dnET, double dnF){
  // These terms collect items and makes things more collected
  double DFA = (-qq*pf/(nA*nA))*( rr*nET + ss*nC + tt*nCD + uu*nF );
  double DFET = qq*pf*rr/nA;
  double DFNC = qq*pf*ss/nA;
  double DFCH = qq*pf*tt/nA;
  double DFNF = qq*pf*uu/nA;
  if( nA != 0.0 ){
    return sqrt( dnA*dnA*DFA*DFA + dnET*dnET*DFET*DFET + dnC*dnC*DFNC*DFNC + dnCD*dnCD*DFCH*DFCH + dnF*dnF*DFNF*DFNF );
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
	double SolenoidScale = 0; // Current direction in solenoid
	double TorusScale = 0; // Current direction in torus
	string Epoch; // Run range 'epoch' that this run is a part of; used for DF/PF calculation
    public:
	// Constructor
	Run() = default;
	// Mutators
	void SetRunInfo(int rn, string tt, string bc, double be, int ne, int hwp, double tp, double ss, double ts, string epoch){
	    RunNumber = rn; TargetType = tt; BeamCurrent = bc; BeamEnergy = be; NumberEvents = ne;
	    HWPStatus = hwp; TargetPolarization = tp; SolenoidScale = ss; TorusScale = ts; Epoch = epoch;
	}
	// Accessors
	int getRunNumber() const{ return RunNumber;}
	string getTargetType() const{ return TargetType;}
	string getBeamCurrent() const{ return BeamCurrent;}
	double getBeamEnergy() const{ return BeamEnergy;}
	int getNumberEvents() const{ return NumberEvents;}
	int getHWPStatus() const{ return HWPStatus;}
	double getTargetPolarization() const{ return TargetPolarization;}
	double getSolenoidScale() const{ return SolenoidScale;}
	double getTorusScale() const{ return TorusScale;}
	string getEpoch() const{ return Epoch;}
	// Destructor
	~Run() =  default;
};

// This class holds the run for a given run period
class RunPeriod{
    private:
	vector<Run> Runs;
    public:
	// Constructor
	RunPeriod() = default;
	// Mutator
	void SetRunPeriod(){
	    ifstream fin( "Input_Text_Files/RGC_Run_Info.txt" );
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
	vector<int> getRunEpoch(string epoch) const{
	    vector<int> RunsInEpoch;
	    for(size_t i=0; i<Runs.size(); i++){
		if( Runs[i].getEpoch() == epoch ) RunsInEpoch.push_back( Runs[i].getRunNumber() );
	    }
	    if( RunsInEpoch.size() == 0 ) cout << "Couldn't find runs for epoch "<<epoch<<endl;
	    return RunsInEpoch;
	}
	// Destructor
	~RunPeriod() = default;
};

// This vector contains the bin boundaries for the Q2 bins
const vector<double> Q2_Bin_Bounds = { /*0.9188, 1.0969,*/ 1.3094, 1.5632, 1.8661, 2.2277,
	2.6594, 3.1747, 3.7899, 4.5243, 5.4009, 6.4475, 7.6969, 9.1884, 10.9689};
// This vector contains the bin boundaries for the x bins
const vector<double> X_Bin_Bounds = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
	0.5, 0.55, 0.6, 0.65, 0.7, 0.75};
//const vector<double> X_Bin_Bounds = {   0.0875, 0.1125, 0.1375, 0.1625, 0.1875, 0.2125, 0.2375, 0.2625,
//	0.2875, 0.3125, 0.3375, 0.3625, 0.3875, 0.4125, 0.4375, 0.4625, 0.4875, 0.5125, 0.5375, 0.5625, 0.5875, 
//	0.6125, 0.6375, 0.6625, 0.6875, 0.7125, 0.7375, 0.7625, 0.7875 };
/*
const vector<double> X_Bin_Bounds = {0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
	0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,
	0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675,
	0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875 };
*/
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
	double getNormNm() const{return NormNm;}
	double getNormNp() const{return NormNp;}
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
		else{
			cout <<"ERROR: Some FC charges are zero for NH3 calculation; check inputs.\n";
		}
		// Calculate ND3 dilution factor for this bin
		// Using CH2 in place of CD2 for the time being
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 ){
			double dfnd3 = ND3DF( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errdfnd3 = ND3DFError( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_ND3 = dfnd3; ErrDF_ND3 = errdfnd3;
		}
		else{
			cout <<"ERROR: Some FC charges are zero for ND3 calculation; check inputs.\n";
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
		else{
			cout <<"ERROR: Some FC charges are zero for NH3 calculation; check inputs.\n";
		}
		// Calculate ND3 packing fraction for this bin
		if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 ){
			double pfnd3 = NH3PF( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errpfnd3 = NH3PFError( nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF, sqrt(nD)/fD, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			PF_bath_ND3 = pfnd3; ErrPF_bath_ND3 = errpfnd3;
			PF_cell_ND3 = (LHe / 5.0)*pfnd3; ErrPF_cell_ND3 = (LHe / 5.0)*errpfnd3;
		}
		else{
			cout <<"ERROR: Some FC charges are zero for ND3 calculation; check inputs.\n";
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
		    if( fA>0 && fCH>0 && fC>0 && fET>0 && fF>0 ){
			double dfnh3 = NH3DF_Thru_PF( thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF );
			double errdfnh3 = NH3DF_Thru_PF_Error(thisPF, nA/fA, nCH/fCH, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCH)/fCH, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_FixedPF_NH3 = dfnh3; ErrDF_FixedPF_NH3 = errdfnh3;
		    }
		    else{
			cout <<"ERROR: Some FC charges are zero for NH3 calculation; check inputs.\n";
		    }
		}
		// Calculate ND3 packing fraction for this bin
		else if( targetType == "ND3" ){
		    if( fD>0 && fCD>0 && fC>0 && fET>0 && fF>0 ){
			double dfnd3 = ND3DF_Thru_PF( thisPF, nD/fD, nCD/fCD, nC/fC, nET/fET, nF/fF );
			double errdfnd3 = ND3DF_Thru_PF_Error(thisPF, nA/fA, nCD/fCD, nC/fC, nET/fET, nF/fF, 0, sqrt(nA)/fA, sqrt(nCD)/fCD, sqrt(nC)/fC, sqrt(nET)/fET, sqrt(nF)/fF );
			DF_FixedPF_ND3 = dfnd3; ErrDF_FixedPF_ND3 = errdfnd3;
			cout << dfnd3 <<" +- "<< errdfnd3 << endl;
		    }
		    else{
			cout <<"ERROR: Some FC charges are zero for ND3 calculation; check inputs.\n";
		    }
		}
		else cout <<"ERROR: Invalid target type. Can only calculate DF(PF) for NH3 and ND3 targets.\n";
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
	    // Now get the average PF for the dataset
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
	    if( denominator != 0.0 ){
	        double avgPF = totalPF / denominator;
	        cout << "Average PF for data set is: "<< avgPF << endl;
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
		fout<<"Bin_Q2,Bin_X,DF_NH3,Err_DF_NH3,DF_Fixed_PF_NH3,ErrDF_FixedPF_NH3,DF_ND3,Err_DF_ND3,PF_bath_NH3,ErrPF_bath_NH3,PF_cell_NH3,ErrPF_cell_NH3,";
		fout<<"FC_total,N_NH3,N_C,N_CH2,N_ET,N_F,N_CD2,FC_NH3,FC_C,FC_CH2,FC_ET,FC_F,FC_CD2,";
		fout<<"normN_normNH3,normN_C,normN_CH2,normN_ET,normN_F,normN_CD2\n";
		for(size_t i=0; i<AllBins.size(); i++){
	 	    for(size_t j=0; j<AllBins[i].size(); j++){
		      if( AllBins[i][j].getDF_NH3() != 0.0 || AllBins[i][j].getDF_ND3() != 0.0){
			fout<<AllBins[i][j].getBinQ2()<<","<<AllBins[i][j].getBinX()<<","<<AllBins[i][j].getDF_NH3()<<","<<AllBins[i][j].getErrDF_NH3()<<",";
			fout<<AllBins[i][j].getDF_FixedPF_NH3()<<","<<AllBins[i][j].getErrDF_FixedPF_NH3()<<",";
			fout<<AllBins[i][j].getDF_ND3()<<","<<AllBins[i][j].getErrDF_ND3()<<","<<AllBins[i][j].getPF_bath_NH3()<<","<<AllBins[i][j].getErrPF_bath_NH3()<<",";
			fout<<AllBins[i][j].getPF_cell_NH3()<<","<<AllBins[i][j].getErrPF_cell_NH3()<<",";
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
    if( thruPF == "YES" )
	AllData.CalculateDFThruPF( targetType );
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
    mgDF->GetXaxis()->SetLimits(0,0.8); mgDF->GetYaxis()->SetRangeUser(0.0,0.5);
    
    string pltTitle = targetType+" DF Plots "+sector;
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
    mgPFbath->GetXaxis()->SetLimits(0,0.8); mgPFbath->GetYaxis()->SetRangeUser(0.0,1.0);
    
    title = targetType+" PF_{cell} for "+sector+"; X; "+targetType+" PF_{cell}";
    mgPFcell->SetTitle(title.c_str());
    mgPFcell->GetXaxis()->SetLimits(0,0.8); mgPFcell->GetYaxis()->SetRangeUser(0.0,1.0);

    string pltTitle = targetType+" PF Plots "+sector;
    TCanvas* PF_Plot = new TCanvas(pltTitle.c_str(),pltTitle.c_str(),1400,600);
    PF_Plot->Divide(2,1);
    PF_Plot->cd(1);
    mgPFbath->Draw("AP");
    mgLegbath->Draw("same");
    PF_Plot->cd(2);
    mgPFcell->Draw("AP");
    mgLegcell->Draw("same");

    return PF_Plot;
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

TGraphErrors* Make_NH3_Sector_Plots( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600};//, 840, 880, 900, 616, 820, 432, 920, 860 };
    //vector<TGraphErrors*> All_NH3_Plots;

    TMultiGraph* mgNH3 = new TMultiGraph();
    int s = stoi(sector);
    vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
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
    }
    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
    gr->SetMarkerStyle(kFullCircle); gr->SetMarkerColor( palette[s-1] );
    //All_NH3_Plots.push_back( gr );

    string title = "DF_{NH3} for "+sector+"; X; DF_{NH3}";
    gr->SetTitle(title.c_str());
    gr->GetXaxis()->SetRangeUser(0,0.8); mgNH3->GetYaxis()->SetRangeUser(-0.25,0.8);

    //TCanvas* Plot = new TCanvas("Plot","Plot",800,600);
    //Plot->cd();
    //mgNH3->Draw("AP");
    return gr;
}

TMultiGraph* Make_ND3_Bin_Plots( DataSet& AllData, string sector ){
 
    vector<int> palette = { 1, 632, 800, 400, 416, 600, 840, 880, 900, 616, 820, 432, 920, 860 };
    //vector<TGraphErrors*> All_ND3_Plots;

    TMultiGraph* mgND3 = new TMultiGraph();

    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	//string title = "Q2_Bin = "+to_string(qmid)+" "+sector;
	vector<double> xbins, dfnh3vals, zeros, dfnh3errs;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){
	    double xmid = (X_Bin_Bounds[j]+X_Bin_Bounds[j+1])/2.0;
	    //Bin* thisBin = AllData.getThisBin( qmid, xmid );
	    //double thisDFND3 = thisBin->getDF_ND3();
	    //double thisErrDFND3 = thisBin->getErrDF_ND3();
	    double thisDFND3 = AllData.getDF_ND3(qmid, xmid);
	    double thisErrDFND3 = AllData.getErrDF_ND3(qmid, xmid);
	    if( thisDFND3 != 0.0 ){
		xbins.push_back(xmid); zeros.push_back(0.0);
		dfnh3vals.push_back( thisDFND3 );
		dfnh3errs.push_back( thisErrDFND3 );
	    }
	}
	if( xbins.size() > 0 ){
	    TGraphErrors* gr = new TGraphErrors(xbins.size(), xbins.data(), dfnh3vals.data(), zeros.data(), dfnh3errs.data() );
	    gr->SetMarkerStyle(kFullCircle); gr->SetMarkerColor( palette[i] );
	    //All_ND3_Plots.push_back( gr );
	    mgND3->Add( gr, "p" );
	}
    }

    string title = "DF_{ND3} for "+sector+"; X; DF_{ND3}";
    mgND3->SetTitle(title.c_str());
    mgND3->GetXaxis()->SetRangeUser(0,0.8); mgND3->GetYaxis()->SetRangeUser(-0.25,0.8);

    //TCanvas* Plot = new TCanvas("Plot","Plot",800,600);
    //Plot->cd();
    //mgND3->Draw("AP");
    return mgND3;

}
