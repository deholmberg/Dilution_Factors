#include "Input_Functions.h"

void Test(){
//int main(){

	cout <<"Test the Count class\n";
	Count Blank; Count Ammonia("NH3");
	Blank.SetNormCounts(); // should fail
	Blank.AddCounts(100,90);
	Blank.AddFCCharge( 4000, 3500);
	Blank.SetNormCounts(); // should work
	Blank.Print();
	Ammonia.AddCounts(100,90);
	Ammonia.AddFCCharge( 4000, 3500);
	Ammonia.Print();
	// Change target type
	Blank.SetTarget("ND3");
	Blank.Print();

    cout << "Making the DataSet object\n";
    DataSet AllData;//, DataS1, DataS2, DataS3, DataS4, DataS5, DataS6;
    //vector<DataSet> Sectors = {DataS1, DataS2, DataS3, DataS4, DataS5, DataS6};
    cout << "Finished making the DataSet object\n";

    // This vector is used to loop over the runs used by Darren
    //vector<double> StableRuns = { 16318, 16323, 16337, 16341, 16308, 16309, 16290, 16296, 16292, 16300, 16301, 16302, 16303, 16194 };
    vector<double> StableRuns = {
    16194, 16290, 16291, 16292, 16293, 16296, 16297, 16298, 16299, 16300, 
    16301, 16302, 16303, 16307, 16308, 16309,  16697, 16698, 16699, 16700, 16701, 
    16702, 16704 };

    vector<double> NH3_2A2 = {
    16138, 16144, 16145, 16146, 16148, 16156, 16157, 16158, 16164, 16166, 16167, 16168, 16169, 16170, 
    16178 };

    vector<double> NH3_2B1 = {16211, 16213, 16214, 16221, 16222, 16223, 16224, 16225, 16226, 16228, 16231, 16232, 
    16233, 16234, 16235, 16236, 16238, 16243, 16244, 16245, 16246, 16248, 16249, 16250, 16251, 16252, 
    16253, 16256, 16257, 16259, 16260 };

    vector<double> NH3_2C1 = { 16318, 16320, 16321, 16322, 16323, 16325, 16326, 16237, 
    16328, 16329, 16330, 16331, 16332, 16333, 16335, 16336, 16337, 16338, 16339, 16341, 16343, 16345, 
    16346, 16348, 16350, 16352, 16353, 16354, 16355, 16356, 16357 };

    vector<double> NH3_7A1 = { 
    16664, 16665, 16666, 16671, 16672, 16673, 16674,//};
    
    //vector<double> NH3_7A2 = {
    16675, 16676, 16678, 16679, 16681, 16682, 16683, 16685, 16686, 16687, 
    16688, 16689, 16690, 16692, 16693, 16695};

    vector<double> NH3_7B1 = {
    16709, 16710, 16711, 16712, 16713, 16715, 16716, 16717, 16718, 16719, 16720, 16721, 16722, 16723, 
    16726, 16727, 16728, 16729, 16730, 16731, 16732, 16733, 16734, 16736, 16738,//};

    //vector<double> NH3_7B2 = {
    16742, 16743, 16744, 16746, 16747, 16748, 16749, 16750, 16751, 16752, //};

    //vector<double> NH3_7B3 = {
    16753, 16754, 16755, 16756, 16757, 16758, 16759, 16761, 16762, 16763, 16765, 16766,//};

    //vector<double> NH3_7B4 = {
    16767, 16768, 16769, 16770, 16771, 16772};


    //vector<double> StableRuns = { 16308, 16309, 16290, 16296, 16292, 16300, 16301, 16302, 16303, 16194 };
/*
    for(int i=1; i<7; i++ ){
        //for(int j=16100; j<16773; j++){
	for(size_t k=0; k<StableRuns.size(); k++){
	    double j = StableRuns[k];
	    stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data_S"<< i <<".txt";
	    stringstream s2; s2 << "Text_Files/C_"<< j <<"_DF_Data_S"<< i <<".txt";
	    stringstream s3; s3 << "Text_Files/CH2_"<< j <<"_DF_Data_S"<< i <<".txt";
	    stringstream s4; s4 << "Text_Files/ET_"<< j <<"_DF_Data_S"<< i <<".txt";
	    stringstream s5; s5 << "Text_Files/F_"<< j <<"_DF_Data_S"<< i <<".txt";
	    stringstream s6; s6 << "Text_Files/ND3_"<< j <<"_DF_Data_S"<< i <<".txt";
	    //cout << j << endl;
	    Sectors[i-1].AddToBins( s1.str(), "NH3" );
	    Sectors[i-1].AddToBins( s6.str(), "ND3" );
	    Sectors[i-1].AddToBins( s2.str(), "C" );
    	    Sectors[i-1].AddToBins( s3.str(), "CH2"   );
	    Sectors[i-1].AddToBins( s4.str(), "ET"  );
	    Sectors[i-1].AddToBins( s5.str(), "F"   );
        }

    }
*/
/*
    //for(int j=16100; j<16773; j++){
    for(size_t k=0; k<StableRuns.size(); k++){
        double j = StableRuns[k];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        stringstream s2; s2 << "Text_Files/C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << "Text_Files/CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << "Text_Files/ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << "Text_Files/F_"<< j <<"_DF_Data.txt";
        AllData.AddToBins( s1.str(), "NH3" );
        AllData.AddToBins( s2.str(), "C" );
        AllData.AddToBins( s3.str(), "CH2"   );
        AllData.AddToBins( s4.str(), "ET"  );
        AllData.AddToBins( s5.str(), "F"   );
    }
    
    for(size_t i=0; i<NH3_2A2.size(); i++){
	double j = NH3_2A2[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        AllData.AddToBins( s1.str(), "NH3" );
    }
 */

    DataSet FirstEp, SecondEp, ThirdEp, FourthEp, FifthEp;
    vector<DataSet> Epochs = {FirstEp, SecondEp, ThirdEp, FourthEp, FifthEp};
    //for(int j=16100; j<16773; j++){
    for(size_t k=0; k<StableRuns.size(); k++){
        double j = StableRuns[k];
        //stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        stringstream s2; s2 << "Text_Files/C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << "Text_Files/CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << "Text_Files/ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << "Text_Files/F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<Epochs.size(); i++){
          //AllData.AddToBins( s1.str(), "NH3" );
          Epochs[i].AddToBins( s2.str(), "C" );
          Epochs[i].AddToBins( s3.str(), "CH2"   );
	  Epochs[i].AddToBins( s3.str(), "CD2" );
          Epochs[i].AddToBins( s4.str(), "ET"  );
          Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }

    DataSet ND3_First, ND3_Second;
    vector<DataSet> ND3_Epochs = {ND3_First, ND3_Second};
    for(size_t k=0; k<StableRuns.size(); k++){
        double j = StableRuns[k];
        //stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        stringstream s2; s2 << "Text_Files/C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << "Text_Files/CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << "Text_Files/ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << "Text_Files/F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<ND3_Epochs.size(); i++){
          //AllData.AddToBins( s1.str(), "NH3" );
          ND3_Epochs[i].AddToBins( s2.str(), "C" );
          ND3_Epochs[i].AddToBins( s3.str(), "CH2"   );
	  ND3_Epochs[i].AddToBins( s3.str(), "CD2" );
          ND3_Epochs[i].AddToBins( s4.str(), "ET"  );
          ND3_Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }
    for(int i=16261; i<=16272; i++){
	stringstream s1; s1 << "Text_Files/ND3_"<< i <<"_DF_Data.txt";
	ND3_Epochs[0].AddToBins( s1.str(), "ND3" );
    }
    for(int i=16273; i<=16289; i++){
	stringstream s1; s1 << "Text_Files/ND3_"<< i <<"_DF_Data.txt";
	ND3_Epochs[1].AddToBins( s1.str(), "ND3" );
    }

    
    for(size_t i=0; i<NH3_2A2.size(); i++){
	double j = NH3_2A2[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[0].AddToBins( s1.str(), "NH3" );
    }
     for(size_t i=0; i<NH3_2B1.size(); i++){
	double j = NH3_2B1[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[1].AddToBins( s1.str(), "NH3" );
    }
     for(size_t i=0; i<NH3_2C1.size(); i++){
	double j = NH3_2C1[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[2].AddToBins( s1.str(), "NH3" );
    }
    for(size_t i=0; i<NH3_7A1.size(); i++){
	double j = NH3_7A1[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[3].AddToBins( s1.str(), "NH3" );
    }
    for(size_t i=0; i<NH3_7B1.size(); i++){
	double j = NH3_7B1[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[4].AddToBins( s1.str(), "NH3" );
    }

    //TFile* outfile = new TFile("testout.root","RECREATE");
    /*
    for(size_t i=0; i<Epochs.size(); i++){
	// These two commands will calculate the DF using a unique PF for each cell
	//Epochs[i].CalculateDF();
	//Epochs[i].CalculateDFThruPF();
	Epochs[i].CalculatePF();
	PF_Legend_Plot( Epochs[i], "NH3 All Sectors", "NH3")->Write();
	//PF_Legend_Plot( Epochs[i], "All Sectors", "ND3")->Write();

	//string outname = "TS/Epoch_Number_"+to_string(i+1)+".csv";
	//Epochs[i].WriteToCSV(outname);
    }
    */
    // Make plots for the ND3 data
    auto DFPF_ND3_E1 = DF_Legend_Plot( ND3_Epochs[1], "All Sectors", "ND3", "YES");
    DFPF_ND3_E1->Draw();
    //auto DF_ND3_E1 = DF_Legend_Plot( ND3_Epochs[1], "All Sectors", "ND3", "NO");
    //DF_ND3_E1->Draw();
/*
    auto PF_ND3_E1 = PF_Legend_Plot( ND3_Epochs[1], "All Sectors", "ND3");
    PF_ND3_E1->Draw();
    auto DF_NH3_E1 = DF_Legend_Plot( Epochs[1], "All Sectors", "NH3", "YES");
    DF_NH3_E1->Draw();
    auto PF_NH3_E1 = PF_Legend_Plot( Epochs[1], "All Sectors", "NH3");
    PF_NH3_E1->Draw();

    for(size_t i=0; i<ND3_Epochs.size(); i++){
	ND3_Epochs[i].CalculatePF();
	PF_Legend_Plot( ND3_Epochs[i], "ND3 All Sectors", "ND3")->Write();
    }
    outfile->Close();
*/
    // Now see if everything was read in
    //AllData.NormalizeAllCounts();
    //AllData.CalculateDF();
    //AllData.CalculatePF();
    //AllData.Print();

    //for(size_t i=0; i<Sectors.size(); i++) Sectors[i].CalculateDF();

    //AllData.WriteToTXT("All_DF_Data.txt");

    // Make the plots
/*    
    TMultiGraph* mgNH3 = Make_NH3_Bin_Plots( AllData, "All Sectors" );
    TCanvas* PlotNH3 = new TCanvas("PlotNH3","PlotNH3",800,600);
    PlotNH3->cd();
    mgNH3->Draw("AP");
    
    AllData.CalculateDFThruPF();
    //AllData.CalculateDF();
    //AllData.CalculatePF();
    AllData.PrintTotalCounts();
    AllData.WriteToCSV("TS/FirstEpoch.csv");

    TCanvas* NH3_Plot = NH3_Legend_Plot( AllData, "All Sectors");
    NH3_Plot->Draw();

    TCanvas* PF_Plot = PF_Legend_Plot( AllData, "All Sectors");
    PF_Plot->Draw();

    TCanvas* PFcell_Plot = PFcell_Legend_Plot( AllData, "All Sectors");
    PFcell_Plot->Draw();
*/
    //cout << AllData.getAveragePF_NH3() << endl;

    //ultiGraph* mgNH3 = Make_NH3_Bin_Plots( AllData, "All Sectors" );
    /*
    TMultiGraph* mgNH3 = Make_PF_NH3_Bin_Plots( AllData, "All Sectors" );
    TCanvas* PlotNH3 = new TCanvas("PlotNH3","PlotNH3",800,600);
    PlotNH3->cd();
    mgNH3->Draw("AP");
    cout << "PF_avg = " << AllData.getAvgPF() << endl;
    
    TMultiGraph* mgND3 = Make_ND3_Bin_Plots( AllData, "All Sectors" );
    TCanvas* PlotND3 = new TCanvas("PlotND3","PlotND3",800,600);
    PlotND3->cd();
    mgND3->Draw("AP");
     
    // This creates TGraphErrors objects for each of the sectors and then draws them together
    TMultiGraph* mgAllSec = new TMultiGraph();
    TLegend* mgLeg = new TLegend(0.6,0.6,0.9,0.9);
    mgLeg->SetNColumns(3);
    for(size_t i=0; i<Sectors.size(); i++){
	auto gr = Make_NH3_Sector_Plots(Sectors[i], to_string(i+1));
	string entry = "Sector "+to_string(i+1);
	mgAllSec->Add( gr, "p");
	mgLeg->AddEntry( gr, entry.c_str(), "p");
    }
    mgAllSec->SetTitle("DF Across All Sectors; X; DF_{NH3}");
    mgAllSec->GetXaxis()->SetRangeUser(0,0.8);
    mgAllSec->GetYaxis()->SetRangeUser(-0.25,0.8);
    TCanvas* ST = new TCanvas("ST","ST",800,600);
    ST->cd();
    mgAllSec->Draw("ap");
    mgLeg->Draw("same");
    
    TCanvas* More = new TCanvas("More","More",1400,800);
    More->Divide(3,2);
    for(size_t i=0; i<Sectors.size(); i++){
	string title = "Sector "+to_string(i+1);
	TMultiGraph* mgSec = Make_NH3_Bin_Plots( Sectors[i], title );
	More->cd(i+1);
	mgSec->Draw("AP");
    }
    More->Draw();
    */
    //return 0;
    cout << "All done!\n";

}
