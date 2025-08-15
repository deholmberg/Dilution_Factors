//***********************************************************************
//
// Date Created: 1/10/2024
// Author: Derek Holmberg 
// Last Modified: 1/30/2025
//
// This program calculates the dilution factors and packing fractions
// for various "epochs" (ranges of runs) in the RG-C dataset. Currently,
// my program has information for only the data taken in the summer 2022
// dataset (runs 16042-16772) and the fall 2022 dataset (runs 16843-17408).
// Also, the code for the ND3 dilution factors and packing fractions
// doesn't work, so don't trust it lol.
//
//***********************************************************************

#include "Input_Functions.h"

vector<double> WeightedAverage(vector<double>& PbPtVals, vector<double>& PbPtErrs){
    double num = 0; double denom = 0;
    double numNeg = 0; double denomNeg = 0; // For the negative polarizations
    vector<double> Results;
    if( PbPtVals.size() != PbPtErrs.size() ){
	cout << "Error: PbPtVals and PbPtErrs vectors not the same size.\n";
	return Results;
    }
    for(size_t i=0; i<PbPtVals.size(); i++){ // Assumes, of course, that PbPtVals have the same size
      if( PbPtVals[i] > 0 ){
	num += PbPtVals[i] / pow(PbPtErrs[i],2) ;
	denom += 1.0 / pow(PbPtErrs[i],2);
      }
      else if( PbPtVals[i] < 0 ){
	numNeg += PbPtVals[i] / pow(PbPtErrs[i],2) ;
	denomNeg += 1.0 / pow(PbPtErrs[i],2);
      }
    }
    if( denom != 0 && denomNeg != 0 ){
	double avgPbPt = num / denom; Results.push_back(avgPbPt);
	double avgPbPtErr = 1.0 / sqrt( denom ); Results.push_back(avgPbPtErr);
	double avgPbPtNeg = numNeg / denomNeg; Results.push_back(avgPbPtNeg);
	double avgPbPtErrNeg = 1.0 / sqrt( denomNeg ); Results.push_back(avgPbPtErrNeg);
	return Results;
    }
    return Results;
}

void PbPt_Extraction(){

    // The "RunPeriod" class holds the information for each of the runs and to which epoch they belong.
    cout << "Setting run period\n";
    RunPeriod Period; Period.SetRunPeriod(); // Holds the information for all the runs in the dataset.
    cout << "Set run period\n";
    // These vectors hold all the run numbers in each of the epochs.
    vector<int> Epoch2 = Period.getRunEpoch("P02");
    vector<int> Epoch3 = Period.getRunEpoch("P03");
    vector<int> Epoch4 = Period.getRunEpoch("P04");
    vector<int> Epoch5 = Period.getRunEpoch("P05");
    vector<int> Epoch6 = Period.getRunEpoch("P06");
    vector<int> Epoch7 = Period.getRunEpoch("P07");
    vector<int> Epoch8 = Period.getRunEpoch("P08");
    vector<int> Epoch9 = Period.getRunEpoch("P09");
    vector<int> Epoch10 = Period.getRunEpoch("P10");
    vector<int> Epoch11 = Period.getRunEpoch("P11");
    vector<int> Epoch12 = Period.getRunEpoch("P12");
    vector<int> Epoch13 = Period.getRunEpoch("P13");
    vector<int> Epoch14 = Period.getRunEpoch("P14");
    vector<int> Epoch15 = Period.getRunEpoch("P15");
    vector<int> Epoch16 = Period.getRunEpoch("P16");
    vector<int> Epoch17 = Period.getRunEpoch("P17");
    vector<int> Epoch18 = Period.getRunEpoch("P18");
    //vector<int> Epoch19 = Period.getRunEpoch("P19");


    vector<vector<int>> Summer22_NH3_Epochs = {Epoch2, Epoch3, Epoch4, Epoch5, Epoch6, Epoch7};
    vector<vector<int>> Fall22_NH3_Epochs = {Epoch8, Epoch9, Epoch10, Epoch11, Epoch12, Epoch13, Epoch14};
    vector<vector<int>> Spring23_NH3_Epochs = {Epoch15, Epoch16, Epoch17, Epoch18};
    vector<vector<int>> RGC_NH3_Epochs = {Epoch2, Epoch3, Epoch4, Epoch5, Epoch6, Epoch7, Epoch8, Epoch9, Epoch10, Epoch11, Epoch12, Epoch13, Epoch14, Epoch15, Epoch16, Epoch17, Epoch18};

    // These are a handful of runs that we've determined to have stable charge-normalized counts.
    // These runs should be good for a preliminary analysis; they consist of CH2, Carbon, Empty,
    // and Foil targets. 
    // This will likely be replaced with something else later on.
    vector<int> StableRuns = {
        16194, 16290, 16291, 16292, 16293, 16296, 16297, 16298, 16299, 16300, 
        16302, 16303, 16307, 16308, 16309, 16698, 16700, 16701, 
        16702, 16704 
    };
    
    /*
    vector<int> StableRuns = {
    	16194, 16290, 16296, 16292, 16300, 16301, 16302, 16303, 16308, 16309
    };
    */
    
    vector<int> StableRunsFall22 = {
	16975, 16872, 16979, 17180, 17181, 17182, 17183, 17134, 17135
    };

    vector<int> StableRunsSpring23 = {
	17763, 17764, 17765, // Foil runs
	17766, 17767, 17768, // Empty target runs
	//17584, 17585, 17586, 17587, 17588, 17589, 17593, 17594, 17595, 17596, 17597, 17604, 17605, 17606, // Darren's NH3 runs
	//17493, 17494, 17497, 17498, 17502, 17505, 17506, // Darren's ND3 runs
	17566, 17567, 17569, 17744, 17745, 17746, // Darren's C runs
	17625, 17626, 17627, 17635, 17637, 17638, // Darren's CH2 runs
	17516, 17517, 17520, 17669, 17670, 17672  // Darren's CD2 runs
    };

    vector<int> StableRunsRGC = {
        16194, 16290, 16291, 16292, 16293, 16296, 16297, 16298, 16299, 16300, 
        16302, 16303, 16307, 16308, 16309, 16698, 16700, 16701, 
        16702, 16704,
	16975, 16872, 16979, 17180, 17181, 17182, 17183, 17134, 17135,
	17763, 17764, 17765, // Foil runs
	17766, 17767, 17768, // Empty target runs
	17566, 17567, 17569, 17744, 17745, 17746, // Darren's C runs
	17625, 17626, 17627, 17635, 17637, 17638, // Darren's CH2 runs
	17516, 17517, 17520, 17669, 17670, 17672  // Darren's CD2 runs
    };

    cout <<"Making the dataset objects\n";
    // The "DataSet" class holds all the information for each of the bins in x, Q^2, like the
    // FC charge and counts. Member functions of this class are used to calculate DF and PF.
    DataSet SecondEp, ThirdEp, FourthEp, FifthEp, SixthEp, SeventhEp, EighthEp, NinthEp, TenthEp, EleventhEp, TwelfthEp, ThirteenthEp, FourteenthEp;
    DataSet FifteenthEp, SixteenthEp, SeventeenthEp, EighteenthEp;//, NineteenthEp;
    vector<DataSet> Epochs = {SecondEp, ThirdEp, FourthEp, FifthEp, SixthEp, SeventhEp};
    vector<DataSet> Fall22_Epochs = {EighthEp, NinthEp, TenthEp, EleventhEp, TwelfthEp, ThirteenthEp, FourteenthEp};
    vector<DataSet> Spring23_Epochs = {FifteenthEp, SixteenthEp, SeventeenthEp, EighteenthEp};
//    vector<DataSet> RGC_Epochs = {SecondEp, ThirdEp, FourthEp, FifthEp, SixthEp, SeventhEp, EighthEp, NinthEp, TenthEp, EleventhEp, TwelfthEp, 
//				  ThirteenthEp, FourteenthEp, FifteenthEp, SixteenthEp, SeventeenthEp, EighteenthEp};


    //string pathToInput = "Extra_Text_Files/Basic_Cuts/";
    string pathToInput = "Text_Files/";


    // Don't worry about double-counting since the DataSet class skips any invalid file entries.
    cout << "Starting to loop over stable background runs\n";
    for(size_t k=0; k<StableRuns.size(); k++){
        int j = StableRuns[k];
        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<Epochs.size(); i++){
          Epochs[i].AddToBins( s2.str(), "C" );
          Epochs[i].AddToBins( s3.str(), "CH2"   );
          Epochs[i].AddToBins( s4.str(), "ET"  );
          Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }
    
    for(size_t k=0; k<StableRunsFall22.size(); k++){
        int j = StableRunsFall22[k];
        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<Fall22_Epochs.size(); i++){
          Fall22_Epochs[i].AddToBins( s2.str(), "C" );
          Fall22_Epochs[i].AddToBins( s3.str(), "CH2"   );
          Fall22_Epochs[i].AddToBins( s4.str(), "ET"  );
          Fall22_Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }

    for(size_t k=0; k<StableRunsSpring23.size(); k++){
        int j = StableRunsSpring23[k];
        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<Spring23_Epochs.size(); i++){
          Spring23_Epochs[i].AddToBins( s2.str(), "C" );
          Spring23_Epochs[i].AddToBins( s3.str(), "CH2"   );
          Spring23_Epochs[i].AddToBins( s4.str(), "ET"  );
          Spring23_Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }

    // Loop through and add the appropriate NH3 runs to each of the epochs
    cout <<"Starting to loop over NH3 runs\n";
    for(size_t i=0; i<Summer22_NH3_Epochs.size(); i++){
	for(size_t j=0; j<Summer22_NH3_Epochs[i].size(); j++){
	    int run = Summer22_NH3_Epochs[i][j];
	    cout << run <<" ";
            stringstream s1; s1 << pathToInput<<"NH3_"<< run <<"_DF_Data.txt";
            Epochs[i].AddToBins( s1.str(), "NH3" );
	}
	cout << endl;
    }
    
    for(size_t i=0; i<Fall22_NH3_Epochs.size(); i++){
	for(size_t j=0; j<Fall22_NH3_Epochs[i].size(); j++){
	    int run = Fall22_NH3_Epochs[i][j];
	    cout << run <<" ";
            stringstream s1; s1 << pathToInput<<"NH3_"<< run <<"_DF_Data.txt";
            Fall22_Epochs[i].AddToBins( s1.str(), "NH3" );
	}
	cout << endl;
    }

    for(size_t i=0; i<Spring23_NH3_Epochs.size(); i++){
	for(size_t j=0; j<Spring23_NH3_Epochs[i].size(); j++){
	    int run = Spring23_NH3_Epochs[i][j];
	    cout << run <<" ";
            stringstream s1; s1 << pathToInput<<"NH3_"<< run <<"_DF_Data.txt";
            Spring23_Epochs[i].AddToBins( s1.str(), "NH3" );
	}
	cout << endl;
    }

    // Make sure the DF is set using the PF
    for(size_t i=0; i<Epochs.size(); i++) Epochs[i].CalculateDFThruPF("NH3");
    for(size_t i=0; i<Fall22_Epochs.size(); i++) Fall22_Epochs[i].CalculateDFThruPF("NH3");
    for(size_t i=0; i<Spring23_Epochs.size(); i++) Spring23_Epochs[i].CalculateDFThruPF("NH3");

/*
    // Do calculate the DFs for the entire RGC dataset (summer, fall, and spring)
    for(size_t i=0; i<RGC_NH3_Epochs.size(); i++){
	for(size_t j=0; j<RGC_NH3_Epochs[i].size(); j++){
	    int run = RGC_NH3_Epochs[i][j];
	    cout << run <<" ";
            stringstream s1; s1 << pathToInput<<"NH3_"<< run <<"_DF_Data.txt";
            RGC_Epochs[i].AddToBins( s1.str(), "NH3" );
	}
	cout << endl;
    }
*/
/*
    // Use this loop to use background run data from the entire RGC data set
    for(size_t k=0; k<StableRunsRGC.size(); k++){
        int j = StableRunsRGC[k];
        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
	for(size_t i=0; i<RGC_Epochs.size(); i++){
          RGC_Epochs[i].AddToBins( s2.str(), "C" );
          RGC_Epochs[i].AddToBins( s3.str(), "CH2"   );
          RGC_Epochs[i].AddToBins( s4.str(), "ET"  );
          RGC_Epochs[i].AddToBins( s5.str(), "F"   );
	}
    }
*/
/*
    // Use this loop to correlate background runs from summer, fall, spring to an epoch
    // from that same period
    for(size_t i=0; i<RGC_Epochs.size(); i++){
	if( i+2 <= 7 ){ // Summer 22 run period
	    for(size_t k=0; k<StableRuns.size(); k++){
        	int j = StableRuns[k];
	        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        	stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
	        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
	        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
    		RGC_Epochs[i].AddToBins( s2.str(), "C" );
	        RGC_Epochs[i].AddToBins( s3.str(), "CH2"   );
        	RGC_Epochs[i].AddToBins( s4.str(), "ET"  );
	        RGC_Epochs[i].AddToBins( s5.str(), "F"   );
    	    }
	}
	else if( i+2 <= 14 ){ // Fall 22 run period
	    for(size_t k=0; k<StableRunsFall22.size(); k++){
        	int j = StableRunsFall22[k];
	        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        	stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
	        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
	        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
    		RGC_Epochs[i].AddToBins( s2.str(), "C" );
	        RGC_Epochs[i].AddToBins( s3.str(), "CH2"   );
        	RGC_Epochs[i].AddToBins( s4.str(), "ET"  );
	        RGC_Epochs[i].AddToBins( s5.str(), "F"   );
    	    }
	}
	else if( i+2 <= 18 ){ // Spring 23 (excludes outbending)
	    for(size_t k=0; k<StableRunsSpring23.size(); k++){
        	int j = StableRunsSpring23[k];
	        stringstream s2; s2 << pathToInput<<"C_"<< j <<"_DF_Data.txt";
        	stringstream s3; s3 << pathToInput<<"CH2_"<< j <<"_DF_Data.txt";
	        stringstream s4; s4 << pathToInput<<"ET_"<< j <<"_DF_Data.txt";
	        stringstream s5; s5 << pathToInput<<"F_"<< j <<"_DF_Data.txt";
    		RGC_Epochs[i].AddToBins( s2.str(), "C" );
	        RGC_Epochs[i].AddToBins( s3.str(), "CH2"   );
        	RGC_Epochs[i].AddToBins( s4.str(), "ET"  );
	        RGC_Epochs[i].AddToBins( s5.str(), "F"   );
    	    }
	}
    }
    for(size_t i=0; i<RGC_Epochs.size(); i++) RGC_Epochs[i].CalculateDFThruPF("NH3");
*/
    

    //Epochs[1].CalculateDFThruPF("NH3");
   // Epochs[2].CalculateDFThruPF("NH3");

    
    // These lines are used to plot Noemie's elastic values of PbPt
    TF1* PbPt_Pos_Su22 = new TF1("PbPt_Pos","0.71",16100,16773); PbPt_Pos_Su22->SetLineColor(kRed);
    TF1* PbPt_Neg_Su22 = new TF1("PbPt_Neg","-0.66",16100,16773);PbPt_Neg_Su22->SetLineColor(kRed);
    PbPt_Pos_Su22->SetLineWidth(5);
    PbPt_Neg_Su22->SetLineWidth(5);
    TF1* PbPt_Pos_Fa22 = new TF1("PbPt_Pos","0.71",16800,17408); PbPt_Pos_Fa22->SetLineColor(kRed);
    TF1* PbPt_Neg_Fa22 = new TF1("PbPt_Neg","-0.70",16800,17408);PbPt_Neg_Fa22->SetLineColor(kRed);
    PbPt_Pos_Fa22->SetLineWidth(5);
    PbPt_Neg_Fa22->SetLineWidth(5);
    TF1* PbPt_Pos_Sp23 = new TF1("PbPt_Pos","0.67",17500,18100); PbPt_Pos_Sp23->SetLineColor(kRed);
    TF1* PbPt_Neg_Sp23 = new TF1("PbPt_Neg","-0.62",17500,18100);PbPt_Neg_Sp23->SetLineColor(kRed);
    PbPt_Pos_Sp23->SetLineWidth(5);
    PbPt_Neg_Sp23->SetLineWidth(5);

/*
    TMultiGraph* EveryEpochSu22_PbPt = All_Epochs_ElasticPbPt( Summer22_NH3_Epochs, Epochs, "NH3", Period );
    EveryEpochSu22_PbPt->SetTitle("P_{b}P_{t} for Summer 2022; Run Number; P_{b}P_{t}");
    EveryEpochSu22_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);


    TCanvas* EP_Su22 = new TCanvas("EP_Su22","EP_Su22",800,600);
    EP_Su22->cd();
    EveryEpochSu22_PbPt->Draw("ap");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");

    TMultiGraph* EveryEpochFa22_PbPt = All_Epochs_ElasticPbPt( Fall22_NH3_Epochs, Fall22_Epochs, "NH3", Period );
    EveryEpochFa22_PbPt->SetTitle("P_{b}P_{t} for Fall 2022; Run Number; P_{b}P_{t}");
    EveryEpochFa22_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);

    TCanvas* EP_Fa22 = new TCanvas("EP_Fa22","EP_Fa22",800,600);
    EP_Fa22->cd();
    EveryEpochFa22_PbPt->Draw("ap");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");

    TMultiGraph* EveryEpochSp23_PbPt = All_Epochs_ElasticPbPt( Spring23_NH3_Epochs, Spring23_Epochs, "NH3", Period );
    EveryEpochSp23_PbPt->SetTitle("P_{b}P_{t} for Spring 2023; Run Number; P_{b}P_{t}");
    EveryEpochSp23_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);

    TCanvas* EP_Sp23 = new TCanvas("EP_Sp23","EP_Sp23",800,600);
    EP_Sp23->cd();
    EveryEpochSp23_PbPt->Draw("ap");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");

    TMultiGraph* EveryEpochRGC_PbPt = All_Epochs_ElasticPbPt( RGC_NH3_Epochs, RGC_Epochs, "NH3", Period );
    EveryEpochRGC_PbPt->SetTitle("P_{b}P_{t} for All RGC Data; Run Number; P_{b}P_{t}");
    EveryEpochRGC_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);

    TCanvas* EP_RGC = new TCanvas("EP_RGC","EP_RGC",800,600);
    EP_RGC->cd();
    EveryEpochRGC_PbPt->Draw("ap");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");
    PbPt_Pos_Fa22->Draw("same");
    PbPt_Neg_Fa22->Draw("same");
    PbPt_Pos_Sp23->Draw("same");
    PbPt_Neg_Sp23->Draw("same");


    cout << "All done!\n";


    //vector<DataSet> ModelEpochs = Epochs;
    vector<DataSet> ModelEpochs = RGC_Epochs;
    for(size_t i=0; i<ModelEpochs.size(); i++){
        ModelEpochs[i].CalculateDFThruPF("NH3");
    }
 
    string toFile = THIS_DIR + "Darren_Model_DF.txt";
    // Loop over the model data first
    ifstream fin2(toFile.c_str()); string line;
    while(getline(fin2,line)){
	stringstream sin(line);
	double q2, x, qbin, xbin, df, dferr;
	sin >> q2 >> x >> qbin >> xbin >> df >> dferr;
	for(size_t i=0; i<ModelEpochs.size(); i++){
            ModelEpochs[i].SetDFs( q2+0.1, x+0.01, df, 0 );
	}
    }
    fin2.close();
*/
    // Calculate the PbPt values for the summer 2022 data and get the weighted average for that dataset
    vector<PbPt> AllA1ValsSu22, ModelA1ValsSu22;
    vector<double> PbPtValsSu22, PbPtErrsSu22;
    for( auto vec : Summer22_NH3_Epochs ){
      for( int run : vec ){
	PbPt thisRunPbPt, thisModelPbPt;
	int thisEpoch = Period.getEpochNum(run);
	cout <<"Calculate PbPt and AllPhys for run "<< run <<" "<< thisEpoch<< " using data DF\n";
	thisRunPbPt.CalculateAllPhys( run, "NH3", Epochs[thisEpoch-2], Period );
	AllA1ValsSu22.push_back( thisRunPbPt );
	PbPtValsSu22.push_back( thisRunPbPt.getBT_Pol() ); PbPtErrsSu22.push_back( thisRunPbPt.getBT_Pol_Err() );

	//cout <<"Calculate PbPt and AllPhys for run "<< run <<" "<< thisEpoch<< " using model DF\n";
 	//thisModelPbPt.CalculateAllPhys( run, "NH3", ModelEpochs[thisEpoch-2], Period );
	//ModelA1ValsSu22.push_back( thisModelPbPt );       
      }
    }
    //TCanvas* AllA1Su22 = new TCanvas("AllA1Su22","AllA1Su22",800,600);
    //vector<TCanvas*> PlotsSu22 = Make_A1_NH3_Plots( AllA1ValsSu22, Epochs[2], "Summer 2022 Data", "Blue" );
    vector<double> WeightedAvgSu22 = WeightedAverage( PbPtValsSu22, PbPtErrsSu22 );
    string AvgPbPtSu22Val = to_string( WeightedAvgSu22[0] );
    TF1* AvgPbPtSu22 = new TF1("AvgPbPtSu22",AvgPbPtSu22Val.c_str(),16100,16800);
    AvgPbPtSu22->SetLineColor(kBlack); AvgPbPtSu22->SetLineWidth(3);
    string AvgPbPtSu22ValNeg = to_string( WeightedAvgSu22[2] );
    TF1* AvgPbPtSu22Neg = new TF1("AvgPbPtSu22Neg",AvgPbPtSu22ValNeg.c_str(),16100,16800);
    AvgPbPtSu22Neg->SetLineColor(kBlack); AvgPbPtSu22Neg->SetLineWidth(3);
   
    // Get PbPt for the fall 2022 dataset
    vector<PbPt> AllA1ValsFa22, ModelA1ValsFa22;
    vector<double> PbPtValsFa22, PbPtErrsFa22;
    for( auto vec : Fall22_NH3_Epochs ){
      for( int run : vec ){
	PbPt thisRunPbPt, thisModelPbPt;
	int thisEpoch = Period.getEpochNum(run);
	cout <<"Calculate PbPt and AllPhys for run "<< run <<" "<< thisEpoch<< " using data DF\n";
	thisRunPbPt.CalculateAllPhys( run, "NH3", Fall22_Epochs[thisEpoch-8], Period );
	AllA1ValsFa22.push_back( thisRunPbPt );
	PbPtValsFa22.push_back( thisRunPbPt.getBT_Pol() ); PbPtErrsFa22.push_back( thisRunPbPt.getBT_Pol_Err() );
      }
    }
    //TCanvas* AllA1Fa22 = new TCanvas("AllA1Fa22","AllA1Fa22",800,600);
    //vector<TCanvas*> PlotsFa22 = Make_A1_NH3_Plots( AllA1ValsFa22, Epochs[2], "Fall 2022 Data", "Blue" );
    vector<double> WeightedAvgFa22 = WeightedAverage( PbPtValsFa22, PbPtErrsFa22 );
    string AvgPbPtFa22Val = to_string( WeightedAvgFa22[0] );
    TF1* AvgPbPtFa22 = new TF1("AvgPbPtFa22",AvgPbPtFa22Val.c_str(),16800,17500);
    AvgPbPtFa22->SetLineColor(kBlack); AvgPbPtFa22->SetLineWidth(3);
    string AvgPbPtFa22ValNeg = to_string( WeightedAvgFa22[2] );
    TF1* AvgPbPtFa22Neg = new TF1("AvgPbPtFa22Neg",AvgPbPtFa22ValNeg.c_str(),16800,17500);
    AvgPbPtFa22Neg->SetLineColor(kBlack); AvgPbPtFa22Neg->SetLineWidth(3);


    vector<PbPt> AllA1ValsSp23, ModelA1ValsSp23;
    vector<double> PbPtValsSp23, PbPtErrsSp23;
    for( auto vec : Spring23_NH3_Epochs ){
      for( int run : vec ){
	PbPt thisRunPbPt, thisModelPbPt;
	int thisEpoch = Period.getEpochNum(run);
	cout <<"Calculate PbPt and AllPhys for run "<< run <<" "<< thisEpoch<< " using data DF\n";
	thisRunPbPt.CalculateAllPhys( run, "NH3", Spring23_Epochs[thisEpoch-15], Period );
	AllA1ValsSp23.push_back( thisRunPbPt );
	PbPtValsSp23.push_back( thisRunPbPt.getBT_Pol() ); PbPtErrsSp23.push_back( thisRunPbPt.getBT_Pol_Err() );
      }
    }
    //TCanvas* AllA1Sp23 = new TCanvas("AllA1Sp23","AllA1Sp23",800,600);
    //vector<TCanvas*> PlotsSp23 = Make_A1_NH3_Plots( AllA1ValsSp23, Epochs[2], "Spring 2023 Data", "Blue" );
    vector<double> WeightedAvgSp23 = WeightedAverage( PbPtValsSp23, PbPtErrsSp23 );
    string AvgPbPtSp23Val = to_string( WeightedAvgSp23[0] );
    TF1* AvgPbPtSp23 = new TF1("AvgPbPtSp23",AvgPbPtSp23Val.c_str(),17500,18100);
    AvgPbPtSp23->SetLineColor(kBlack); AvgPbPtSp23->SetLineWidth(3);
    string AvgPbPtSp23ValNeg = to_string( WeightedAvgSp23[2] );
    TF1* AvgPbPtSp23Neg = new TF1("AvgPbPtSp23Neg",AvgPbPtSp23ValNeg.c_str(),17500,18100);
    AvgPbPtSp23Neg->SetLineColor(kBlack); AvgPbPtSp23Neg->SetLineWidth(3);

/*
    TMultiGraph* EveryEpochRGC_PbPt = All_Epochs_ElasticPbPt( RGC_NH3_Epochs, RGC_Epochs, "NH3", Period );
    EveryEpochRGC_PbPt->SetTitle("P_{b}P_{t} for All RGC Data; Run Number; P_{b}P_{t}");
    EveryEpochRGC_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
*/

    TCanvas* EP_Su22 = new TCanvas("EP_Su22","EP_Su22",800,600);
    EP_Su22->cd();
    TMultiGraph* EveryEpochSu22_PbPt = All_Epochs_ElasticPbPt( Summer22_NH3_Epochs, Epochs, "NH3", Period );
    EveryEpochSu22_PbPt->SetTitle("P_{b}P_{t} for Summer 2022; Run Number; P_{b}P_{t}");
    EveryEpochSu22_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
    EveryEpochSu22_PbPt->Draw("ap");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");
    AvgPbPtSu22->Draw("same");
    AvgPbPtSu22Neg->Draw("same");

    TCanvas* EP_Fa22 = new TCanvas("EP_Fa22","EP_Fa22",800,600);
    EP_Fa22->cd();
    TMultiGraph* EveryEpochFa22_PbPt = All_Epochs_ElasticPbPt( Fall22_NH3_Epochs, Fall22_Epochs, "NH3", Period );
    EveryEpochFa22_PbPt->SetTitle("P_{b}P_{t} for Fall 2022; Run Number; P_{b}P_{t}");
    EveryEpochFa22_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
    EveryEpochFa22_PbPt->Draw("ap");
    PbPt_Pos_Fa22->Draw("same");
    PbPt_Neg_Fa22->Draw("same");
    AvgPbPtFa22->Draw("same");
    AvgPbPtFa22Neg->Draw("same");

    TCanvas* EP_Sp23 = new TCanvas("EP_Sp23","EP_Sp23",800,600);
    EP_Sp23->cd();
    TMultiGraph* EveryEpochSp23_PbPt = All_Epochs_ElasticPbPt( Spring23_NH3_Epochs, Spring23_Epochs, "NH3", Period );
    EveryEpochSp23_PbPt->SetTitle("P_{b}P_{t} for Spring 2023; Run Number; P_{b}P_{t}");
    EveryEpochSp23_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
    EveryEpochSp23_PbPt->Draw("ap");
    PbPt_Pos_Sp23->Draw("same");
    PbPt_Neg_Sp23->Draw("same");
    AvgPbPtSp23->Draw("same");
    AvgPbPtSp23Neg->Draw("same");

/*
    TMultiGraph* EveryEpochFa22_PbPt = All_Epochs_ElasticPbPt( Fall22_NH3_Epochs, Fall22_Epochs, "NH3", Period );
    EveryEpochFa22_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
    TMultiGraph* EveryEpochSp23_PbPt = All_Epochs_ElasticPbPt( Spring23_NH3_Epochs, Spring23_Epochs, "NH3", Period );
    EveryEpochSp23_PbPt->GetYaxis()->SetRangeUser(-0.8,0.9);
    EveryEpochSu22_PbPt->Draw("ap");
    EveryEpochFa22_PbPt->Draw("same");
    EveryEpochSp23_PbPt->Draw("same");
    PbPt_Pos_Su22->Draw("same");
    PbPt_Neg_Su22->Draw("same");
    PbPt_Pos_Fa22->Draw("same");
    PbPt_Neg_Fa22->Draw("same");
    PbPt_Pos_Sp23->Draw("same");
    PbPt_Neg_Sp23->Draw("same");
    AvgPbPtSu22->Draw("same");
    AvgPbPtFa22->Draw("same");
    AvgPbPtSp23->Draw("same");
*/


    vector<PbPt> AllA1ValsSummer22;
    for( auto vec : Summer22_NH3_Epochs ){
      for( int run : vec ){
	PbPt thisRunPbPt;
	int thisEpoch = Period.getEpochNum(run);
	cout <<"Calculate PbPt and AllPhys for run "<< run <<" "<< thisEpoch<< " using data DF\n";
	thisRunPbPt.CalculateAllPhys( run, "NH3", Epochs[thisEpoch-2], Period );
	AllA1ValsSummer22.push_back( thisRunPbPt );
      }
    }
    TCanvas* AllA1Summer22 = new TCanvas("AllA1Summer22","AllA1Summer22",800,600);
    vector<TCanvas*> PlotsSummer22 = Make_A1_NH3_Plots( AllA1ValsSummer22, Epochs[2], "All Summer22 Data", "Blue" );


//vector<TGraphErrors*> Plots = Make_A1_NH3_Plots( AllA1Vals, Epochs[2], "Summer 2022 Data", "Blue" );

    //Plots[0]->Draw();
    //Plots[1]->Draw();
    //AllA1->cd(); //Plots->SetMarkerColor(kBlue); Plots->SetMarkerStyle(kFullCircle);
    //Plots[0]->Draw("ap"); 
    //Plots[1]->Draw("same");
    //Plots[1]->SetLineWidth(5);
    
    //vector<TCanvas*> ModelPlots = Make_A1_NH3_Plots( ModelA1Vals, Epochs[2], "Summer 2022 Data" );
    //vector<TGraphErrors*> ModelPlots = Make_A1_NH3_Plots( ModelA1Vals, Epochs[2], "Summer 2022 Data", "Black" );

    //ModelPlots[0]->Draw();
    //ModelPlots[1]->Draw();
    //ModelPlots->SetLineColor(kBlack); ModelPlots->SetMarkerStyle(kFullSquare);
    //ModelPlots->Draw("same");
/*   
    TMultiGraph* mg = new TMultiGraph();
    mg->Add( Plots[0], "p" );
    mg->Add( Plots[1], "l" );
    //mg->Add( ModelPlots[0], "p" );
    mg->SetTitle("A_{1,p} for Summer 2022, Data and Model; X; A_{1,p}");
    mg->Draw("AP");
*/
    // Now make plots of the averaged DF for each of the epochs

/*
    TCanvas* All_DF_Su22_Plots = DF_Q2_Averaged_Plot( Epochs, "Summer", "NH3", "YES" );  
    All_DF_Su22_Plots->Draw();

    TCanvas* All_DF_Fa22_Plots = DF_Q2_Averaged_Plot( Fall22_Epochs, "Fall", "NH3", "YES" );  
    All_DF_Fa22_Plots->Draw();

    TCanvas* All_DF_Sp23_Plots = DF_Q2_Averaged_Plot( Spring23_Epochs, "Spring", "NH3", "YES" );  
    All_DF_Sp23_Plots->Draw();
*/

    //vector<TGraphErrors*> OldData = Old_A1_Data();

}
