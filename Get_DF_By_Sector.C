//***********************************************************************
//
// Date Created: 8/30/2024
// Author: Derek Holmberg 
// Last Modified: 1/10/2025
//
// The purpose of this program is to do my own analysis of the dilution
// factors (DF) following the procedure outlined in Sebastian's paper.
// This program specifically checks the DF for each sector in the CLAS12
// detector.
// 
//***********************************************************************

#include "Input_Functions.h"
#include "QADB.h"

using namespace std;
using namespace QA;

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
   string RASTERname = TARGET_TYPE +"_"+to_string(runToAnalyze)+"_RASTER"; string RASTERtitle = "Plot of Beam Raster Position; V_{X} (cm); V_{Y} (cm);";
   TH1D* ElecVX = new TH1D(VXname.c_str(), VXtitle.c_str(),60,-3,3);
   TH1D* ElecVY = new TH1D(VYname.c_str(), VYtitle.c_str(),60,-3,3);
   TH1D* ElecVZ = new TH1D(VZname.c_str(), VZtitle.c_str(),400,-20,20);
   TH1D* ElecW  = new TH1D(Wname.c_str(),  Wtitle.c_str(), 500, 0,5);
   TH1D* ElecQ2 = new TH1D(Q2name.c_str(), Q2title.c_str(),1000,0,10);
   TH2D* ElecRASTER = new TH2D(RASTERname.c_str(), RASTERtitle.c_str(),90,-1.5,1.5,90,-1.5,1.5);
   // Create the ROOT file for storing the diagnostic plots
   string tfile_path = thisDir+"/ROOT_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_Diagnostic_Plots.root";
   TFile* file = new TFile(tfile_path.c_str(), "RECREATE");

   // Read in the information for the runs
   //RunPeriod Su22; Su22.SetRunPeriod("Su22");
   RunPeriod Period; Period.SetRunPeriod();
   double thisTPol = Period.getTargetPolarization( runToAnalyze );
   // This is a correction factor for the target polarization
   double corr_factor = 1;
   if( thisTPol < 0 ) corr_factor = -1; 

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
   qa->SetMaskBit("SectorLoss",true); // this is the only bit we check here
   qa->SetMaskBit("LowLiveTime",true);
   qa->SetMaskBit("Misc",false);

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
    hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
    hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
    hipo::bank CONF(factory.getSchema("RUN::config"));
    hipo::bank RASTER(factory.getSchema("RASTER::position"));

    while(reader.next()==true){ // #1
	reader.read(event);
	event.getStructure(PART);
	event.getStructure(HEL);
	event.getStructure(TRACK);
	event.getStructure(CONF);
        event.getStructure(HEL_SCALER);
        event.getStructure(RUN_SCALER);
	event.getStructure(RASTER);
	int evnum = CONF.getInt("event",0);
	int runnum= CONF.getInt("run",0); // Already read in, but this is a double-check

      if(qa->Pass(runnum,evnum)) {
	//PART.show();
	if( HEL_SCALER.getRows() > 0 ){
          for(int row=0; row<HEL_SCALER.getRows(); row++){
	      int hel = corr_factor * HEL_SCALER.getByte("helicity",row);
	      if(hel == -1){
		     FCup_hel_n += HEL_SCALER.getFloat("fcupgated",row);
	      }
	      else if(hel == 1){
		     FCup_hel_p += HEL_SCALER.getFloat("fcupgated",row);
	      }
	      else if(hel == 0){
		     FCup_hel_0 += HEL_SCALER.getFloat("fcupgated",row);
	      }
	  }
        }
	
	int trigger_status = PART.getInt("status",0); // Trigger particle status variable
	int trigger_pid = PART.getInt("pid",0); // Trigger particle pid; looking for electron, so pid = 11
	double trigger_chi2pid = PART.getDouble("chi2pid",0);
	int hel = corr_factor * HEL.getInt("helicity",0); // Get the helicity for this event
	int trk_sector = TRACK.getInt("sector",0); // Get the sector of the DC that the electron went thru
	int trk_pindex = TRACK.getInt("pindex",0); // Gets pindex to make sure we're actually grabbing the electron in the track bank

	if( trigger_pid == 11 && trigger_status < -2000 && abs(hel) == 1 && trk_pindex == 0 ){ // (#2) Trigger electron in Forward detector w/ defined helicity state
	   double px = PART.getDouble("px",0); double py = PART.getDouble("py",0); double pz = PART.getDouble("pz",0);
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
	   // The following applies kinematic cuts that Gregory used in his analysis:
	   // E >= 2.6 GeV, 5 deg < theta < 35 deg for scattered e-, abs(Vz+4.5) > 4, W >= 2.0
	   if( E >= 2.6 && theta > 5.0 && theta < 35.0 && w > 2.0 && Q2 > 1.0 ){
		qa->AccumulateCharge(); // Accumulate FC charge from QADB
		double vx = PART.getDouble("vx",0); double vy = PART.getDouble("vy",0);  // electron x,y vertices
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
	   } // end of "if" checking the kinematic cuts
	} // end #2
	event_counter++;
      } // End of QA pass "if" statement
    } // End of "reader.next()" while loop #1
   } // end of "if" statement checking if the file opened
   else cout << "Failed to open file " << fullFilePath <<". Please check path.\n";
   FIN.close();
   cout << "Now starting to write to the output files.\n";

   // Write FC charge to all bins
   //double fullFC = (qa->GetAccumulatedCharge()) / 2.0;
   AllData.SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
   for(size_t i=0; i<SectorData.size(); i++) SectorData[i].SlotChargeIntoBin( FCup_hel_n, FCup_hel_p );
   //AllData.SlotChargeIntoBin( fullFC, fullFC );
   //for(size_t i=0; i<SectorData.size(); i++) SectorData[i].SlotChargeIntoBin( fullFC, fullFC );
   // Write the data
   WriteSectorData( AllData, SectorData, TARGET_TYPE, runToAnalyze, thisDir );

   // Write histograms to Tfile
   ElecVX->Write(); ElecVY->Write(); ElecVZ->Write();
   ElecW->Write(); ElecQ2->Write(); ElecRASTER->Write();
   file->Close();

}
