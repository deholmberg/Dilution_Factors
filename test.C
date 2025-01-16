#include "Input_Functions.h"


void test(){

	DataSet AllTest;

	cout << AllTest.getAllTheoryNH3( 1.4, 0.09 ) << endl;
        cout << AllTest.getAllTheoryNH3( 8.0, 0.61 ) << endl;	

	PbPt Test;
	Test.ReadInThisNH3Run(16211);
	Test.Print();
}
