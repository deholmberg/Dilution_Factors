//***********************************************************************
//
// Date Created: 8/30/2024
// Author: Derek Holmberg 
// Last Modified: 6/20/2025
//
// The purpose of this program is to do my own analysis of the dilution
// factors (DF) following the procedure outlined in Sebastian's paper.
// This program specifically checks the DF for each sector in the CLAS12
// detector.
// 
//***********************************************************************

#include "Input_Functions.h"
#include "QADB.h"
#include "Fiducial_Cuts.h"

using namespace std;
using namespace QA;

// Returns angle in degrees
double getPhi(double px, double py){
    if( px > 0 && py > 0 ) return 180.0 * atan(py / px) / M_PI;       // first quadrant
    else if( px<0 && py>0 ) return 180.0+ 180.0*atan(py / px) / M_PI;// second quadrant
    else if( px<0 && py<0 ) return 180.0+ 180.0*atan(py / px) / M_PI; // third quadrant
    else if( px>0 && py<0 ) return 360.0+ 180.0*atan(py / px) / M_PI;// fourth quadrant
    else{ // either px or py is zero 
	std::cout << "Crap dammit there was a zero! I put the angle at -20000 degrees for error purposes lol.\n";
	return -20000;
    }
}

// Checks if the particle passed the HTCC cut of nphe > 2
double getHTCCnphe(int pIndex, bank& HTCC ){
    for(int i=0; i<HTCC.getRows(); i++){
	if( HTCC.getInt("pindex",i) == pIndex && HTCC.getInt("detector",i) == 15 ){
	    return HTCC.getDouble("nphe",i);
	}
    }
    return 0; // If there's no particle there
}

// The following function is used to write each sector's data set
void WriteSectorData( DataSet& AllBins, vector<DataSet>& SectorData, string TARGET_TYPE, int runToAnalyze, string thisDir ){
  string outname;
  cout << outname << endl;
  //outname = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Text_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_DF_Data.txt";
  outname = thisDir+"/Text_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_DF_Data.txt";
  AllBins.WriteRawDFInputs(outname);
  
  for(int i=1; i<7; i++){
    //outname = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Text_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_DF_Data_S"+to_string(i)+".txt";
    outname = thisDir+"/Text_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_DF_Data_S"+to_string(i)+".txt";
    cout << outname << endl;
    SectorData[i-1].WriteRawDFInputs(outname);
  }
}
/*
// This function is used for slotting counts into the correct histogram based on Q2 binning
void SlotIntoQ2Bin(double vz, double theta, double momentum, double& X_POS_DC[3], double& Y_POS_DC[3], double& Z_POS_DC[3],
vector<TH1D*>& VZ, vector<TH1D*>& THETA, vector<TH1D*> MOMEN, vector<TH2D*>& DCXYR1 double qmid){
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	if( qmid >= Q2_Bin_Bounds[i] && qmid < Q2_Bin_Bounds[i+1] ){
	    VZ[i]->Fill(vz);
	    THETA[i]->Fill(theta);
	    MOMEN[i]->Fill(momentum);
	    DCXYR1[i]->Fill( X_POS_DC[0], Y_POS_DC[0] );
	    break;
	}
    }
}
// This function is used to write the above histograms to a TFile object
void WriteQ2Bin( vector<TH1D*>& VZ, vector<TH1D*>& THETA, vector<TH1D*> MOMEN, vector<TH2D*>& DCXYR1 ){
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	VZ[i]->Write();
	THETA[i]->Write();
	MOMEN[i]->Write();
	DCXYR1[i]->Write();
    }
}
*/
// This function is used for slotting counts into the correct histogram based on Q2 binning
void SlotIntoQ2Bin(double vz, double theta, double momentum, double X_POS_DC[3], double Y_POS_DC[3], double Z_POS_DC[3], double Edge[3], double TRK_CHI2,
vector<TH1D*>& VZ, vector<TH1D*>& THETA, vector<TH1D*> MOMEN, vector<TH2D*>& DCXYR1, vector<TH2D*>& TRKCHI2EDGE, double qmid){
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	if( qmid >= Q2_Bin_Bounds[i] && qmid < Q2_Bin_Bounds[i+1] ){
	    VZ[i]->Fill(vz);
	    THETA[i]->Fill(theta);
	    MOMEN[i]->Fill(momentum);
	    DCXYR1[i]->Fill( X_POS_DC[0], Y_POS_DC[0] );
	    if( TRK_CHI2 != 0 && Edge[0] != 0 ) TRKCHI2EDGE[i]->Fill( Edge[0], TRK_CHI2 );
	    break;
	}
    }
}
// This function is used to write the above histograms to a TFile object
void WriteQ2Bin( vector<TH1D*>& VZ, vector<TH1D*>& THETA, vector<TH1D*> MOMEN, vector<TH2D*>& DCXYR1, vector<TH2D*>& TRKCHI2EDGER1 ){
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	VZ[i]->Write();
	THETA[i]->Write();
	MOMEN[i]->Write();
	DCXYR1[i]->Write();
	TRKCHI2EDGER1[i]->Write();
    }
}
// This function takes one of the kinematic distributions and plots all of them on the same canvas
template <class T>
void PlotTogether( vector<T>& Plots, string CanvasName ){

    // First, adjust the number of entries in the tcanvas
    int R = 0; int C = 0; // rows and columns, respectively
    while( R*C < Plots.size() ){
	C++;
	if( R*C < Plots.size() ) R++;
    }

    TCanvas* All_Plots = new TCanvas(CanvasName.c_str(), CanvasName.c_str(), 1400, 1000 );
    All_Plots->Divide(C,R);

    for(size_t i=0; i<Plots.size(); i++){
	All_Plots->cd(i+1);
	Plots[i]->Draw("colz");
    }
    All_Plots->Write();

}

// This function grabs the DC hit position for a given particle
void GetDCHitPos( double X_POS_DC[3], double Y_POS_DC[3], double Z_POS_DC[3], double edgeRegion[3], double& TRK_CHI2,  bank& TRAJ, bank& TRACK, int pindex ){
                //cout <<" Loop through TRACK and get the reduced chi2 for the DC track"<<endl;
    double chi2 = 0; double ndf = 0;
                for(int j=0; j<TRACK.getRows(); j++){
                    if( TRACK.getInt("pindex",j) == pindex ){
                        chi2 = TRACK.getFloat("chi2",j); ndf = 1.0*TRACK.getInt("NDF",j);
                        //sectorDC = TRACK.getInt("sector",j);
                    }
                }
    
    if( ndf != 0 ) TRK_CHI2 = chi2/ndf;

    //double X_POS_DC[3] = {0,0,0}; double Y_POS_DC[3] = {0,0,0}; double Z_POS_DC[3] = {0,0,0};
            for(int j=0; j<TRAJ.getRows(); j++){
                if( TRAJ.getInt("pindex",j) == pindex ){
                    int detector = TRAJ.getInt("detector",j); // 6 corresponds to DC
                    int layer = TRAJ.getInt("layer",j); // 6 is R1, 18 is R2, 36 is R3
                    if( detector == 6 && layer == 6 ){
                        edgeRegion[0] = TRAJ.getFloat("edge",j);
			X_POS_DC[0] = TRAJ.getFloat("x",j);
			Y_POS_DC[0] = TRAJ.getFloat("y",j);
			Z_POS_DC[0] = TRAJ.getFloat("z",j);
                    }
                    else if( detector == 6 && layer == 18 ){
                        edgeRegion[1] = TRAJ.getFloat("edge",j);
			X_POS_DC[1] = TRAJ.getFloat("x",j);
			Y_POS_DC[1] = TRAJ.getFloat("y",j);
			Z_POS_DC[1] = TRAJ.getFloat("z",j);
                    }
                    else if( detector == 6 && layer == 36 ){
                        edgeRegion[2] = TRAJ.getFloat("edge",j);
			X_POS_DC[2] = TRAJ.getFloat("x",j);
			Y_POS_DC[2] = TRAJ.getFloat("y",j);
			Z_POS_DC[2] = TRAJ.getFloat("z",j);
                    }
                }
            }
}

// These two TLorentzVector objects hold the information for the beam electron and proton target
// (px, py, pz, E)
// Everything in this program ignores the electron mass, which is negligible in this GeV range
const TLorentzVector BeamElectron = TLorentzVector(0, 0, beam_energy, beam_energy);
const TLorentzVector ProtonTarget = TLorentzVector(0, 0, 0, p_mass);

// This is the main loop that goes through the data
void Get_DF_By_Sector(){

   string path_to_sidis, TARGET_TYPE, thisDir;
   int runToAnalyze;
   cin >> path_to_sidis >> TARGET_TYPE >> runToAnalyze >> thisDir;
   // The following cuts off the "/slurm" part of the string
   cout << thisDir << endl;
   thisDir.erase( thisDir.size() - 6, thisDir.size() );
   cout << thisDir << endl;
   int event_counter = 0; // Counts total number of events

   // These diagnostic plots are used to keep track of the electron vertex positions,
   // energy distributions, and raster positions
   // Vertex Plots:
   string VXname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_VX"; string VXtitle = "Plot of e^{-} V_{X}; V_{X} (cm);;";
   string VYname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_VY"; string VYtitle = "Plot of e^{-} V_{Y}; V_{Y} (cm);;";
   string VZname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_VZ"; string VZtitle = "Plot of e^{-} V_{Z}; V_{Z} (cm);;";
   string Wname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_W"; string Wtitle = "Plot of e^{-} Missing Mass; W (GeV);;";
   string Q2name = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_Q2"; string Q2title = "Plot of e^{-} Q^{2}; Q^{2} (GeV^{2});;";
   string Thetaname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_Theta"; string Thetatitle = "Plot of e^{-} Polar Angle; #theta (degrees);;";
   string Phiname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_Phi"; string Phititle = "Plot of e^{-} Azimuthal Angle; #phi (degrees);;";
   string Xname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_X"; string Xtitle = "Plot of e^{-} Bjorken X; X;;";
   string VSname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_Q2_vs_X"; string VStitle = "e^{-} Q^{2} vs. Bjorken X;X; Q^{2} (GeV^{2})";
   string RASTERname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_RASTER"; string RASTERtitle = "Plot of Beam Raster Position; V_{X} (cm); V_{Y} (cm);";
   string PvsNphename = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_PvsNphe";
   string ECALname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_ECAL"; string ECALtitle = "E_{PCAL}/P vs. (E_{ECinner}+E_{ECouter})/P; (E_{ECinner}+E_{ECouter})/P; E_{PCAL}/P";
   string ECAL_DiagCutname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_ECAL_DiagCut"; string ECAL_DiagCuttitle = "E_{PCAL}/P vs. (E_{ECinner}+E_{ECouter})/P (Diag. Cut); (E_{ECinner}+E_{ECouter})/P; E_{PCAL}/P";
   string NPHEname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_NPHE"; string NPHEtitle = "Photoelectrons in HTCC; Nphe;;";
   string Chi2PIDname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_Chi2PID"; string Chi2PIDtitle = "#chi^{2} PID of Trigger e^{-}; #chi^{2} PID;;";
   string VxVyName = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_VxVy"; string VxVyTitle = "Electron V_{Y} vz. V_{X}; V_{X} (cm); V_{Y} (cm)";

   TH1D* ElecVX = new TH1D(VXname.c_str(), VXtitle.c_str(),60,-3,3);
   TH1D* ElecVY = new TH1D(VYname.c_str(), VYtitle.c_str(),60,-3,3);
   TH1D* ElecVZ = new TH1D(VZname.c_str(), VZtitle.c_str(),400,-20,20);
   TH1D* ElecW  = new TH1D(Wname.c_str(),  Wtitle.c_str(), 500, 0,5);
   TH1D* ElecQ2 = new TH1D(Q2name.c_str(), Q2title.c_str(),1200,0,12);
   TH1D* ElecTheta = new TH1D(Thetaname.c_str(), Thetatitle.c_str(),180,0,180);
   TH1D* ElecPhi = new TH1D(Phiname.c_str(), Phititle.c_str(),360,0,360);
   TH1D* ElecX = new TH1D(Xname.c_str(), Xtitle.c_str(), 100,0,1);
   TH2D* ElecQ2vsX = new TH2D(VSname.c_str(), VStitle.c_str(), 100,0,1,1200,0,12);
   TH2D* ElecRASTER = new TH2D(RASTERname.c_str(), RASTERtitle.c_str(),90,-1.5,1.5,90,-1.5,1.5);
   TH2D* ElecECAL_DiagCut = new TH2D(ECAL_DiagCutname.c_str(), ECAL_DiagCuttitle.c_str(),80,0,0.4,80,0,0.4);
   TH2D* ElecECAL = new TH2D(ECALname.c_str(), ECALtitle.c_str(),80,0,0.4,80,0,0.4);

   TH1D* ElecNPHE = new TH1D(NPHEname.c_str(), NPHEtitle.c_str(),80,0,40);
   TH1D* ElecChi2PID = new TH1D(Chi2PIDname.c_str(), Chi2PIDtitle.c_str(),200,-10,10);
   TH2D* ElecVxVy = new TH2D(VxVyName.c_str(), VxVyTitle.c_str(), 60,-3,3, 60,-3,3);

   // These create loops over the bins in Q2 that I use, and track what several kinematics look
   // like for each bin in Q2
   vector<TH1D*> Vz_Per_Q2, Theta_Per_Q2, Momentum_Per_Q2;
   vector<TH2D*> DC_Hit_Per_Q2, TrackChi2Edge_Per_Q2;
   for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){
	double qmid = (Q2_Bin_Bounds[i]+Q2_Bin_Bounds[i+1])/2.0;
	string vzName = "Vz_for_Qbin"+to_string(i+1);
	string vzTitle = "V_{Z} Dist. for "+to_string(Q2_Bin_Bounds[i])+"#leq Q^{2} <"+to_string(Q2_Bin_Bounds[i+1])+"; V_{Z} (cm); Counts";
	TH1D* vzBin = new TH1D(vzName.c_str(), vzTitle.c_str(), 400, -20,20); Vz_Per_Q2.push_back( vzBin );
	string thetaName = "Theta_for_Qbin"+to_string(i+1);
	string thetaTitle = "#theta Dist. for "+to_string(Q2_Bin_Bounds[i])+"#leq Q^{2} <"+to_string(Q2_Bin_Bounds[i+1])+"; #theta (degree); Counts";
	TH1D* thetaBin = new TH1D(thetaName.c_str(), thetaTitle.c_str(), 40,0 ,40); Theta_Per_Q2.push_back( thetaBin );
	string MomentumName = "P_for_Qbin"+to_string(i+1);
	string MomentumTitle = "P Dist. for "+to_string(Q2_Bin_Bounds[i])+"#leq Q^{2} <"+to_string(Q2_Bin_Bounds[i+1])+"; P (GeV/c); Counts";
	TH1D* MomentumBin = new TH1D(MomentumName.c_str(), MomentumTitle.c_str(), 100, 0,10); Momentum_Per_Q2.push_back( MomentumBin );
	string DC_Hit_Name = "DC_Hit_for_Qbin"+to_string(i+1);
	string DC_Hit_Title = "DC Hit Dist. for "+to_string(Q2_Bin_Bounds[i])+"#leq Q^{2} <"+to_string(Q2_Bin_Bounds[i+1])+"; X (cm); Y (cm)";
	TH2D* DC_Hit_Bin = new TH2D(DC_Hit_Name.c_str(), DC_Hit_Title.c_str(), 400, -200,200, 400, -200,200); DC_Hit_Per_Q2.push_back( DC_Hit_Bin );
	string TrackChi2Edge_Name = "TrackChi2Edge_for_Qbin"+to_string(i+1);
	string TrackChi2Edge_Title = "#chi^{2}_{Track} Dist. for "+to_string(Q2_Bin_Bounds[i])+"#leq Q^{2} <"+to_string(Q2_Bin_Bounds[i+1])+"; Edge (cm); #chi^{2}_{Track}";
	TH2D* TrackChi2Edge_Bin = new TH2D(TrackChi2Edge_Name.c_str(), TrackChi2Edge_Title.c_str(), 100, 0,100, 100,0,100); TrackChi2Edge_Per_Q2.push_back( TrackChi2Edge_Bin );

   }

   int NBins[3] = {1000, 100, 8};
   double Bounds[6] = {0, 0, 0, 200, 10, 40};
   //HistoHorde PvsNphe_NegTrack; PvsNphe_NegTrack.Initialize("NegTrackNphe","NPhe","P","#theta", NBins, Bounds);
   //HistoHorde PvsNphe_Electron; PvsNphe_Electron.Initialize(PvsNphename.c_str(),"Nphe","P","#theta", NBins, Bounds);
   //HistoHorde PvsNphe_NegTrackQ2; PvsNphe_NegTrackQ2.Initialize("NegTrackNpheQ2","NPhe","P","Q^{2}", NBins, XYBounds, Q2_Bin_Bounds );
   //HistoHorde PvsNphe_ElectronQ2; PvsNphe_ElectronQ2.Initialize("ElectronNpheQ2","Nphe","P","Q^{2}", NBins, XYBounds, Q2_Bin_Bounds );

   // Create the ROOT file for storing the diagnostic plots
   string tfile_path = thisDir+"/ROOT_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_Diagnostic_Plots.root";
   TFile* file = new TFile(tfile_path.c_str(), "RECREATE");

   // Read in the information for the runs
   //RunPeriod Su22; Su22.SetRunPeriod("Su22");
   RunPeriod Period; //Period.SetRunPeriod();
   double thisTPol = Period.getTargetPolarization( runToAnalyze );
   // This is a correction factor for the target polarization
   double corr_factor = 1;
   if( thisTPol < 0 ) corr_factor = -1; //Test if this is causing issues

   DataSet AllData, AllDataS1, AllDataS2, AllDataS3, AllDataS4, AllDataS5, AllDataS6;
   vector<DataSet> SectorData = {AllDataS1, AllDataS2, AllDataS3, AllDataS4, AllDataS5, AllDataS6};

   // Gated FCup counts from scalers
   double FCup_hel_n = 0;
   double FCup_hel_p = 0;      
   double FCup_hel_0 = 0;
   double FCup_run = 0;

   // instantiate QADB
   //QADB * qa = new QADB(15000,18000,true);
   QADB* qa = new QADB("latest");

   // custom QA cut definition
   // - decide which defects you want to check for; an event will not pass
   //   the QA cut if the associated QA bin has any of the specified defects
   // - set to true to check the bit
   // - set to false to ignore the bit (by default, all bits are ignored)
   qa->SetMaskBit("TotalOutlier",true);
   qa->SetMaskBit("TerminalOutlier",true);
   qa->SetMaskBit("MarginalOutlier",true);
   qa->SetMaskBit("SectorLoss",true); 
   qa->SetMaskBit("LowLiveTime",true);
   qa->SetMaskBit("Misc", Period.applyMiscBit(runToAnalyze) );
   qa->SetMaskBit("TotalOutlierFT",true);
   qa->SetMaskBit("TerminalOutlierFT",true);
   qa->SetMaskBit("MarginalOutlierFT",true);
   qa->SetMaskBit("LossFT",true);
   qa->SetMaskBit("BSAWrong",true);
   qa->SetMaskBit("BSAUnknown",true);
   qa->SetMaskBit("ChargeHigh",true);
   qa->SetMaskBit("ChargeNegative",true);
   qa->SetMaskBit("ChargeUnknown",true);
   qa->SetMaskBit("PossiblyNoBeam",true);

   hipo::reader  reader;

   string fullFilePath = path_to_sidis;// + to_string(runToAnalyze) + ".hipo";
   ifstream FIN(fullFilePath.c_str());
   if( !FIN.fail() ){
    reader.open(fullFilePath.c_str());
    cout << "Opened file " << fullFilePath << endl;

    // This is the main analysis loop
    hipo::dictionary  factory;
    reader.readDictionary(factory);
    hipo::structure  particles;
    hipo::structure  detectors;
    hipo::event      event;
    hipo::bank  dataPART;

    hipo::bank PART(factory.getSchema("REC::Particle"));
    hipo::bank HEL(factory.getSchema("REC::Event"));
    hipo::bank TRACK(factory.getSchema("REC::Track"));
    hipo::bank TRAJ(factory.getSchema("REC::Traj"));
    hipo::bank CAL(factory.getSchema("REC::Calorimeter"));
    hipo::bank CHKV(factory.getSchema("REC::Cherenkov"));
    hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
    hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
    hipo::bank CONF(factory.getSchema("RUN::config"));
    hipo::bank RASTER(factory.getSchema("RASTER::position"));

    while(reader.next()==true){ // #1
	reader.read(event);
	event.getStructure(PART);
	event.getStructure(HEL);
	event.getStructure(TRACK);
	event.getStructure(TRAJ);
	event.getStructure(CAL);
	event.getStructure(CHKV);
	event.getStructure(CONF);
        event.getStructure(HEL_SCALER);
        event.getStructure(RUN_SCALER);
	event.getStructure(RASTER);
	int evnum = CONF.getInt("event",0);
	int runnum= CONF.getInt("run",0); // Already read in, but this is a double-check

      if(qa->Pass(runnum,evnum)) {
	//PART.show();
	 // Temporarily commented out so I can test the QADB FC charge cuts
	 
	if( HEL_SCALER.getRows() > 0 ){
          //for(int row=0; row<HEL_SCALER.getRows(); row++){
	      int hel = corr_factor * HEL_SCALER.getByte("helicity",0);
	      if(hel == -1){
		     FCup_hel_n += HEL_SCALER.getFloat("fcupgated",0);
	      }
	      else if(hel == 1){
		     FCup_hel_p += HEL_SCALER.getFloat("fcupgated",0);
	      }
	      else if(hel == 0){
		     FCup_hel_0 += HEL_SCALER.getFloat("fcupgated",0);
	      }
	  //}
        }
	
	int trigger_status = PART.getInt("status",0); // Trigger particle status variable
	int trigger_pid = PART.getInt("pid",0); // Trigger particle pid; looking for electron, so pid = 11
	double trigger_chi2pid = PART.getDouble("chi2pid",0);
	int hel = corr_factor * HEL.getInt("helicity",0); // Get the helicity for this event
	int trk_sector = TRACK.getInt("sector",0); // Get the sector of the DC that the electron went thru
	int trk_pindex = TRACK.getInt("pindex",0); // Gets pindex to make sure we're actually grabbing the electron in the track bank
	double calEnergy[3] = {0,0,0}; // Holds the energy deposited by trigger particle in: {PCAL, ECinner, ECouter}
	// These check if the trigger electron meets the fiducial cuts
	bool pcalCutPass = true;//pcal_fiducial_cut( 0, 3, PART, CAL );
	bool dcCutPass = true; //dc_fiducial_cut( 0, PART, TRAJ, 1 );
	if( trigger_pid == 11 && trigger_status <= -2000 && trigger_status > -4000 && abs(hel) == 1 && trk_pindex == 0 && pcalCutPass && dcCutPass && abs(trigger_chi2pid) < 3 ){ // (#2) Trigger electron in Forward detector w/ defined helicity state
	   double px = PART.getDouble("px",0); double py = PART.getDouble("py",0); double pz = PART.getDouble("pz",0);
	   double vx = PART.getDouble("vx",0); double vy = PART.getDouble("vy",0);  // electron x,y vertices
	   double vz = PART.getDouble("vz",0);
	   double pt2 = px*px + py*py + pz*pz; // Total momentum squared
	   double E = sqrt(pt2);
	   TLorentzVector ScatteredElectron = TLorentzVector(px, py, pz, E);
	   // Now calculate the kinematics to find the appropriate bin
	   TLorentzVector VirtualPhoton = BeamElectron - ScatteredElectron;
	   double Q2 = abs( VirtualPhoton * VirtualPhoton ); // Magnitude of the four-momentum transfer squared
	   double x  = Q2 / (2.0* (ProtonTarget*VirtualPhoton));
	   double w  = sqrt( (ProtonTarget+VirtualPhoton)*(ProtonTarget+VirtualPhoton) );
	   double theta = acos(pz/sqrt(pt2)) * 180.0/M_PI ; // Theta of scattered electron in degrees!!!
	   double y = ( beam_energy - E ) / beam_energy;
	   double phi = getPhi( px, py ); // Phi of scattered electron in degrees!!!
	   double nphe = getHTCCnphe( 0, CHKV );
	   bool ecalPass = ecalEnergyCheck( 0, E, calEnergy, CAL ); // Checks if there's energy in the ECAL and sets the calorimeter energy.
	   // The following applies kinematic cuts that Gregory used in his analysis:
	   // E >= 2.6 GeV, 5 deg < theta < 35 deg for scattered e-, abs(Vz+4.5) > 4, W >= 2.0
	   //bool VxVy_cut = vx < 0.2218+0.6493 && vx > 0.2218-0.6493 && vy < 0.1149+0.6647 && vy > 0.1149-0.6647;
	   if( theta >= 5.0 && theta < 40.0 && w > 2.0 && Q2 > 1.0 && E > 2.6 /*&& vz > -10 && vz < 2*/ ){
		qa->AccumulateCharge(); // Accumulate FC charge from QADB
		double rx = RASTER.getDouble("x",0);double ry = RASTER.getDouble("y",0); // beam raster positions
		// Figure out which bin to put the count in
		AllData.SlotCountsIntoBin( Q2, x, hel );
		if( trk_sector == 1 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		else if( trk_sector == 2 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		else if( trk_sector == 3 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		else if( trk_sector == 4 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		else if( trk_sector == 5 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		else if( trk_sector == 6 ) SectorData[trk_sector-1].SlotCountsIntoBin( Q2, x, hel );
		// Now slot into the the histograms
		ElecVX->Fill(vx); ElecVY->Fill(vy); ElecVZ->Fill(vz);
		ElecW->Fill(w); ElecQ2->Fill(Q2); ElecRASTER->Fill(rx,ry);
		ElecTheta->Fill(theta); ElecPhi->Fill(phi);
		ElecX->Fill(x); ElecQ2vsX->Fill( x, Q2 );
		ElecVxVy->Fill( vx, vy );

		double X_POS_DC[3] = {0,0,0}; double Y_POS_DC[3] = {0,0,0}; double Z_POS_DC[3] = {0,0,0};
		double Edge[3] = {0,0,0}; double TRK_CHI2 = 0;
		GetDCHitPos( X_POS_DC, Y_POS_DC, Z_POS_DC, Edge, TRK_CHI2, TRAJ, TRACK, 0 );
		SlotIntoQ2Bin( vz, theta, E, X_POS_DC, Y_POS_DC, Z_POS_DC, Edge, TRK_CHI2, Vz_Per_Q2, Theta_Per_Q2, Momentum_Per_Q2, DC_Hit_Per_Q2, TrackChi2Edge_Per_Q2, Q2 );

		if( E > 4.5 ){ 
		    if( ((calEnergy[1]+calEnergy[2])/E) < (0.2 - (calEnergy[0]/E)) ) ElecECAL_DiagCut->Fill( (calEnergy[1]+calEnergy[2])/E, calEnergy[0]/E );
		}
		else ElecECAL_DiagCut->Fill( (calEnergy[1]+calEnergy[2])/E, calEnergy[0]/E );
		
		ElecECAL->Fill( (calEnergy[1]+calEnergy[2])/E, calEnergy[0]/E );
		ElecNPHE->Fill( nphe ); 
		ElecChi2PID->Fill( trigger_chi2pid );
	   } // end of "if" checking the kinematic cuts
	} // end #2
	event_counter++;
      } // End of QA pass "if" statement
    } // End of "reader.next()" while loop #1
   } // end of "if" statement checking if the file opened
   else cout << "Failed to open file " << fullFilePath <<". Please check path.\n";
   FIN.close();
   cout << "Now starting to write to the output files.\n";
/*
   // Write FC charge to all bins using the FC charge readout at the end of each helicity window.
   AllData.SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
   for(size_t i=0; i<SectorData.size(); i++) SectorData[i].SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
*/
   // Using the QADB doesn't distinguish between different helicity states. To keep agreement with the
   // data structures I wrote in my code, I'll just divide the total FC charge in two, passing each 
   // half into my data structures; this basically treats both the + and - helicity states as having
   // the same measured FC charge.
   
   //double fullFC = (qa->GetAccumulatedCharge()) / 2.0;
   //AllData.SlotChargeIntoBin( fullFC, fullFC );
   //for(size_t i=0; i<SectorData.size(); i++) SectorData[i].SlotChargeIntoBin( fullFC, fullFC );
 
   AllData.SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
   for(size_t i=0; i<SectorData.size(); i++) SectorData[i].SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
  
   // Write the data
   WriteSectorData( AllData, SectorData, TARGET_TYPE, runToAnalyze, thisDir );

   // Write histograms to Tfile
   ElecVX->Write(); ElecVY->Write(); ElecVZ->Write();
   ElecW->Write(); ElecQ2->Write(); ElecRASTER->Write();
   ElecTheta->Write(); ElecPhi->Write();
   ElecX->Write(); ElecQ2vsX->Write();
   //PvsNphe_Electron.Write( PvsNphename, file );
   ElecECAL->Write(); ElecNPHE->Write();
   ElecECAL_DiagCut->Write();
   ElecChi2PID->Write();
   ElecVxVy->Write();
   WriteQ2Bin( Vz_Per_Q2, Theta_Per_Q2, Momentum_Per_Q2, DC_Hit_Per_Q2, TrackChi2Edge_Per_Q2 );
   PlotTogether( DC_Hit_Per_Q2, "DC_Per_Q2" );

   file->Close();

}
