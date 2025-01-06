//***********************************************************************
//***********************************************************************
//  HIPO 4.0
//=======================================================================
// Date Created: 8/30/2024
// // Author: Derek Holmberg 
// Last Modified: 11/21/2024
//
// The purpose of this program is to do my own analysis of the dilution
// factors (DF) following the procedure outlined in Sebastian's paper.
// 
//***********************************************************************

#include "/w/hallb-scshelf2102/clas12/holmberg/Beam_Raster_Calibration/Headers/functions.h"
#include "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Test_Header.h"
#include "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Input_Functions.h"
#include "QADB.h"

using namespace QA;
using namespace std;

// These two TLorentzVector objects hold the information for the beam electron and proton target
// (px, py, pz, E)
// Everything in this program ignores the electron mass, which is negligible in this GeV range
const TLorentzVector BeamElectron = TLorentzVector(0, 0, beam_energy, beam_energy);
const TLorentzVector ProtonTarget = TLorentzVector(0, 0, 0, p_mass);

// This is the main loop that goes through the data
void Get_DF(){

   string path_to_sidis, TARGET_TYPE;//, outfile_path;
   int runToAnalyze;
   cin >> path_to_sidis >> TARGET_TYPE >> runToAnalyze;// >> outfile_path;
   int event_counter = 0; // Counts total number of events

   // This reads in the gated Faraday Cup information
   cout << "Read in the FC information\n";
   string line;

   vector< vector<Bin>> All_Bins = SetKinematicBins(); // Sets up the bins to be used in the analysis
   
   vector<string> q2mid; // Used for storing the mid-values of Q2 in each of the bins based on the bin boundaries

   // Create the latest QADB for this particular run
   QADB* qa = new QADB("latest", runToAnalyze, runToAnalyze);
   // Definitions of QA cuts
   qa->CheckForDefect("TotalOutlier",true);
   qa->CheckForDefect("TerminalOutlier",true);
   qa->CheckForDefect("MarginalOutlier",true);
   qa->CheckForDefect("SectorLoss",true); // this is the only bit we check here
   qa->CheckForDefect("LowLiveTime",true);
   qa->CheckForDefect("Misc",true);
   cout << "\ndefect mask = 0b" << bitset<16>(qa->GetMask()) << endl;

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

    while(reader.next()==true){ // #1
	reader.read(event);
	event.getStructure(PART);
	event.getStructure(HEL);
	//PART.show();

	int trigger_status = PART.getInt("status",0); // Trigger particle status variable
	int trigger_pid = PART.getInt("pid",0); // Trigger particle pid; looking for electron, so pid = 11
	int hel = HEL.getInt("helicity",0); // Get the helicity for this event
	//double trigger_chi2pid = PART.getDouble("chi2pid",0);
	if( trigger_pid == 11 && trigger_status < -2000 && abs(hel) == 1 /*&& trigger_chi2pid > 3.0*/ ){ // (#2) Trigger electron in Forward detector w/ defined helicity state
	   double px = PART.getDouble("px",0); double py = PART.getDouble("py",0); double pz = PART.getDouble("pz",0);
	   double vz = PART.getDouble("vz",0);
	   double pt2 = px*px + py*py + pz*pz; // Total momentum squared
	   //double E = sqrt( pt2 + pow(e_mass,2) );
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
		// Figure out which bin to put the count in
		    for(size_t k=0; k<All_Bins.size(); k++){
			for(size_t l=0; l<All_Bins[k].size(); l++){
			  if(All_Bins[k][l].X_Min != 0.0){
			    if( x <= All_Bins[k][l].X_Max && x > All_Bins[k][l].X_Min && Q2 <= All_Bins[k][l].Q2_Max && Q2 > All_Bins[k][l].Q2_Min ){
			      if( hel == 1 ) All_Bins[k][l].AddCounts( 1.0, 0 ); // Add to the anti-aligned state
			      else if( hel == -1 ) All_Bins[k][l].AddCounts( 0, 1.0 ); // Add to the aligned state
			    }
			  }
			}
		    }
	   } // end of "if" checking the kinematic cuts
	} // end #2

	event_counter++;
    } // End of "reader.next()" while loop #1
   } // end of "if" statement checking if the file opened
   else cout << "Failed to open file " << fullFilePath <<". Please check path.\n";
   FIN.close();



// Now, go through the same file again, getting the FCup charges
//Gated FCup counts from scalers
 double FCup_hel_n = 0;
 double FCup_hel_p = 0;      
 double FCup_hel_0 = 0;
 double FCup_run = 0;

 ifstream FIN2(fullFilePath.c_str());
 if( !FIN2.fail() ){
    reader.open(fullFilePath.c_str());
    cout << "Opened file " << fullFilePath << endl;
    hipo::dictionary  factory;
    reader.readDictionary(factory);
    hipo::structure  particles;
    hipo::structure  detectors;
    hipo::event      event;
    hipo::bank  dataPART;

    hipo::bank HEL_SCALER(factory.getSchema("HEL::scaler"));
    hipo::bank RUN_SCALER(factory.getSchema("RUN::scaler"));
    hipo::bank CONF(factory.getSchema("RUN::config"));
    hipo::bank HEL(factory.getSchema("REC::Event"));
      
      //Start loop on events
      while (reader.next() == true){
      	  //Get hipo info
	  reader.read(event);  
	  event.getStructure(CONF);
	  event.getStructure(HEL_SCALER);
	  event.getStructure(RUN_SCALER);
	  event.getStructure(HEL);

	  //Event info
	  //RunNumber = CONF.getInt("run", 0);
	  //Helicity = HEL.getInt("helicity", 0);
	  
	  //Scalers information
	  //NOTE: "fcup" refers (at least I think) to the ungated FCup charge
          for(int row=0; row<HEL_SCALER.getRows(); row++){
	      int hel = HEL_SCALER.getByte("helicity",row);
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
	  if(RUN_SCALER.getRows()>0){
              FCup_run = RUN_SCALER.getFloat("fcupgated",0);
          }
	}//END LOOP ON EVENTS
 }
 FIN2.close();


  //string outname = "Text_Files/"+TARGET_TYPE +"_"+ to_string(runToAnalyze) +"_DF_Data.txt";
  string outname = "/w/hallb-scshelf2102/clas12/holmberg/Dilution_Factors/Text_Files/"+TARGET_TYPE+"_"+to_string(runToAnalyze)+"_DF_Data.txt";
  ofstream fout(outname.c_str());
  fout << "Q2_Min   Q2_Max   X_Min   X_Max   Nm_Counts   Np_Counts   FCm_Charge   FCp_Charge\n";
  for(size_t i=0; i<All_Bins.size(); i++){
	  for(size_t j=0; j<All_Bins[i].size(); j++){
		  if( All_Bins[i][j].X_Max != 0.0 ){
		      Bin bin = All_Bins[i][j];
		      //string out = All_Bins[i][j].WriteBins();
		      fout << bin.Q2_Min<<"	"<< bin.Q2_Max <<"	"<< bin.X_Min <<"	"<< bin.X_Max <<"	";
		      fout << bin.Nm <<"	"<<bin.Np <<"	"<< FCup_hel_n <<"	"<< FCup_hel_p << endl;
		      //fout << out << endl;
		  }
	  }
  }
  fout.close();



  //return 0;

}
