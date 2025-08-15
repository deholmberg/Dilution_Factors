//***********************************************************************
//
// Date Created: 1/10/2024
// Author: Derek Holmberg 
// Last Modified: 1/10/2025
//
// This program calculates the dilution factors and packing fractions
// for various "epochs" (ranges of runs) in the RG-C dataset. Currently,
// my program has information for only the data taken in the summer 2022
// dataset (runs 16042-16772), so keep that in mind.
// Also, the code for the ND3 dilution factors and packing fractions
// doesn't work, so don't trust it lol.
//
//***********************************************************************

#include "Input_Functions.h"

void Dilution_Factor_Analysis(){

    // The "RunPeriod" class holds the information for each of the runs and to which epoch they belong.
    RunPeriod Period; Period.SetRunPeriod(); // Holds the information for all the runs in the dataset.
    // These vectors hold all the run numbers in each of the epochs.
    vector<int> Epoch2 = Period.getRunEpoch("P02");
    vector<int> Epoch3 = Period.getRunEpoch("P03");
    vector<int> Epoch4 = Period.getRunEpoch("P04");
    vector<int> Epoch5 = Period.getRunEpoch("P05");
    vector<int> Epoch6 = Period.getRunEpoch("P06");

    // These are a handful of runs that we've determined to have stable charge-normalized counts.
    // These runs should be good for a preliminary analysis; they consist of CH2, Carbon, Empty,
    // and Foil targets. 
    // This will likely be replaced with something else later on.
    vector<int> StableRuns = {
    16194, 16290, 16291, 16292, 16293, 16296, 16297, 16298, 16299, 16300, 
    16301, 16302, 16303, 16307, 16308, 16309,  16697, 16698, 16699, 16700, 16701, 
    16702, 16704 };

    // The "DataSet" class holds all the information for each of the bins in x, Q^2, like the
    // FC charge and counts. Member functions of this class are used to calculate DF and PF.
    DataSet FirstEp, SecondEp, ThirdEp, FourthEp, FifthEp;
    vector<DataSet> Epochs = {FirstEp, SecondEp, ThirdEp, FourthEp, FifthEp};
   
    // Don't worry about double-counting, since the DataSet class skips any invalid file entries.
    for(size_t k=0; k<StableRuns.size(); k++){
        int j = StableRuns[k];
        stringstream s2; s2 << "Text_Files/C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << "Text_Files/CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << "Text_Files/ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << "Text_Files/F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<Epochs.size(); i++){
          Epochs[i].AddToBins( s2.str(), "C" );
          Epochs[i].AddToBins( s3.str(), "CH2"   );
          Epochs[i].AddToBins( s4.str(), "ET"  );
          Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }

    // We'll just plot the third epoch right now.
    for(size_t i=0; i<Epoch3.size(); i++){
	int j = Epoch3[i];
        stringstream s1; s1 << "Text_Files/NH3_"<< j <<"_DF_Data.txt";
        Epochs[2].AddToBins( s1.str(), "NH3" );
    }

/*
    // Plot the DF for Epoch2
    auto DF_NH3_E2 = DF_Legend_Plot( Epochs[0], "All Sectors", "NH3", "YES");
    DF_NH3_E2->Draw();
    // Plot the PF for Epoch2
    auto PF_NH3_E2 = PF_Legend_Plot( Epochs[0], "All Sectors", "NH3");
    PF_NH3_E2->Draw();
*/
    //Plot the raw double-spin asymmetries for the NH3 runs in the epoch
    auto AllRaw_E2 = AllRawNH3_Legend_Plot( Epochs[2], "All Sectors");
    AllRaw_E2->Draw();

    // You can also save everything to a .txt or .csv file. 
    // They used to work but I broke something, so hold off on this for now...
    // Epochs[0].WriteToCSV("Text_Files/FirstEpoch.csv");
   
    cout << "All done!\n";

}
