/*************************************************************************************
 *
 * Author: Derek Holmberg
 * Date Created: 3/25/25
 * Last Modified: 3/25/25
 *
 * The purpose of this program is to create scaled versions of the CH2 target data
 * for the summer22 dataset. The purpose of the scaling is to create a pseudo-CD2
 * target set for the summer22 dataset. By multiplying the counts in each x, Q2 bin
 * by the ratio of CD2/CH2 counts in the corresponding bin, the pseudo dataset is
 * created.
 *
 * Why use the scaling instead of just plugging in the CD2 counts from the spring23
 * dataset? Using the scaling method accounts for any systematic differences between
 * the run periods, allowing for an "apples to apples" comparison within the sumemr22
 * run period.
 *
***************************************************************************************/

#include "../Input_Functions.h"

using namespace std;

// This is just a dumb class to hold the scaling data
class ScaleBin{

  public:
    double Q2_Min = 0; double Q2_Max = 0;
    double X_Min = 0;  double X_Max = 0;
    double ScaleFactor = 0;
    
    void SetBin( double q2min, double q2max, double xmin, double xmax, double sf){
	Q2_Min = q2min; Q2_Max = q2max; X_Min = xmin; X_Max = xmax; ScaleFactor = sf;
    }

};

void Make_Scaled_Data(){

    DataSet Spring23_Data;

    string textFilePath;
    cout <<"Enter path to CD2 and CH2 directory: ";
    cin >> textFilePath; cout << endl;

    // Read in the CD2 and CH2 data from spring23
    for(int i=17514; i<17677; i++){
	stringstream s1; s1 << "../"<< textFilePath <<"/CD2_"<< i <<"_DF_Data.txt";
	stringstream s2; s2 << "../"<< textFilePath <<"/CH2_"<< i <<"_DF_Data.txt";
	Spring23_Data.AddToBins( s1.str(), "CD2" );
	Spring23_Data.AddToBins( s2.str(), "CH2" );
    }

    //Spring23_Data.Print();
   
    // Output file holding the scaling factor for each of the bins
    ofstream fout("SummerAndFall22_CH2_Scaling_Factors.txt");
 
    for(size_t i=0; i<Q2_Bin_Bounds.size()-1; i++){ // Loop over Q2 bins
	double qmin = Q2_Bin_Bounds[i];
	double qmax = Q2_Bin_Bounds[i+1];
	double qmid = (qmin+qmax)/2.0;
	for(size_t j=0; j<X_Bin_Bounds.size()-1; j++){ // Loop over X bins
	    double xmin = X_Bin_Bounds[j];
	    double xmax = X_Bin_Bounds[j+1];
	    double xmid = (xmin+xmax)/2.0;
	    Bin thisBin = Spring23_Data.getThisBin(qmid, xmid);
	    // Get the total counts and total FC charges for all target types
	    double Nt_CH2 = thisBin.getNt("CH2"); double FCt_CH2 = thisBin.getFCt("CH2");
	    double Nt_CD2 = thisBin.getNt("CD2"); double FCt_CD2 = thisBin.getFCt("CD2");
	    if( FCt_CH2 > 0 && FCt_CD2 > 0 && Nt_CH2 > 0 && Nt_CD2 > 0 ){
		cout << "Scale factor for bin X = " << xmid <<", Q2 = "<< qmid <<" is: "<< (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) << endl;
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<<  (Nt_CD2/FCt_CD2)/(Nt_CH2/FCt_CH2) << endl;

	    }
	    else
		fout << qmin <<"	"<< qmax <<"	"<< xmin <<"	"<< xmax <<"	"<< 0 << endl;
	}
    }
    fout.close();

    // Now go through the process of scaling the Summer22 CH2 files by the new scaling factors
    for( int i=16298; i<16304; i++){
	stringstream sout; sout << "Scaled_CH2_"<<i<<"_DF_Data.txt";
	stringstream sscale; sscale << "SummerAndFall22_CH2_Scaling_Factors.txt";
	stringstream sch2; sch2 << "../" << textFilePath << "/CH2_"<<i<<"_DF_Data.txt";
	ofstream foutscale( sout.str().c_str() );
	ifstream finscale( sscale.str().c_str() );
	ifstream finch2( sch2.str().c_str() );
	// Both input files should have the same number of rows, so this should work lol
	string line1, line2;
	getline( finch2, line2 ); // Throw away header row
	foutscale << line2; // Rewrite the header row into new file
	while( getline( finch2, line2 ) ){
	    getline( finscale, line1 );
	    stringstream sin1(line1);
	    double qmin1, qmax1, xmin1, xmax1, scale;
	    sin1 >> qmin1 >> qmax1 >> xmin1 >> xmax1 >> scale;
	    stringstream sin2(line2);
	    double qmin2, qmax2, xmin2, xmax2, nm, np, fcm, fcp;
	    sin2 >> qmin2 >> qmax2 >> xmin2 >> xmax2 >> nm >> np >> fcm >> fcp;
	    
	    //cout << qmin1<<"   "<< qmax1<<"   "<< xmin1<<"   "<< xmax1<<"   "<< endl;
	    //cout << qmin2<<"   "<< qmax2<<"   "<< xmin2<<"   "<< xmax2<<"   "<< nm<<"   "<< np<<"   "<< fcm<<"   "<< fcp << endl;

	    foutscale << qmin2 <<"	"<< qmax2 <<"	"<< xmin2 <<"	"<< xmax2 <<"	"<< scale*nm <<"	"<< scale*np <<"	"<< fcm <<"	"<< fcp << endl;
	    if( qmin1 != qmin2 || qmax1 != qmax2 || xmin1 != xmin2 || xmax1 != xmax2 ) cout << "ERROR: Files in run "<<i<<" not the same size.\n";

	}
	foutscale.close(); finscale.close(); finch2.close();
    }

    // Now go through the process of scaling the Fall22 CH2 files by the new scaling factors
    for( int i=17118; i<17408; i++){
	stringstream sout; sout << "Scaled_CH2_"<<i<<"_DF_Data.txt";
	stringstream sscale; sscale << "SummerAndFall22_CH2_Scaling_Factors.txt";
	stringstream sch2; sch2 << "../" << textFilePath << "/CH2_"<<i<<"_DF_Data.txt";
	ifstream finch2( sch2.str().c_str() );
	// Both input files should have the same number of rows, so this should work lol
	if( !finch2.fail() ){
	  ofstream foutscale( sout.str().c_str() );
 	  ifstream finscale( sscale.str().c_str() );
	  string line1, line2;
	  getline( finch2, line2 ); // Throw away header row
	  foutscale << line2; // Rewrite the header row into new file
	  while( getline( finch2, line2 ) ){
	    getline( finscale, line1 );
	    stringstream sin1(line1);
	    double qmin1, qmax1, xmin1, xmax1, scale;
	    sin1 >> qmin1 >> qmax1 >> xmin1 >> xmax1 >> scale;
	    stringstream sin2(line2);
	    double qmin2, qmax2, xmin2, xmax2, nm, np, fcm, fcp;
	    sin2 >> qmin2 >> qmax2 >> xmin2 >> xmax2 >> nm >> np >> fcm >> fcp;
	    
	    //cout << qmin1<<"   "<< qmax1<<"   "<< xmin1<<"   "<< xmax1<<"   "<< endl;
	    //cout << qmin2<<"   "<< qmax2<<"   "<< xmin2<<"   "<< xmax2<<"   "<< nm<<"   "<< np<<"   "<< fcm<<"   "<< fcp << endl;

	    foutscale << qmin2 <<"	"<< qmax2 <<"	"<< xmin2 <<"	"<< xmax2 <<"	"<< scale*nm <<"	"<< scale*np <<"	"<< fcm <<"	"<< fcp << endl;
	    if( qmin1 != qmin2 || qmax1 != qmax2 || xmin1 != xmin2 || xmax1 != xmax2 ) cout << "ERROR: Files in run "<<i<<" not the same size.\n";

	  }
	  foutscale.close(); finscale.close();
	}
	finch2.close();
    }


    cout << "All done!\n";

}
