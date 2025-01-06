
using namespace std;

bool isInList( vector<int>& List, int num){

	bool inList = false;
	for(size_t i=0; i<List.size(); i++){
		if( num == List[i]){
			inList = true;
			break;
		}
	}
	return inList;
}

void Make_List(){

	vector<int> Good_Runs; // Good runs from the RGC wiki page

	ifstream fin("Good_Runs.txt");
	string line;
	while( getline(fin, line) ){
		stringstream sin(line);
		int run; sin >> run;
		Good_Runs.push_back(run);
	}
	fin.close();

	// Now reading in from the RCDB list
	ifstream fin2("RCDB_INFO_RGC.txt");
	ofstream fout("RGC_Run_Lists.csv");
	string delimiter = ",";
	string runToCheck;
	while( getline(fin2, line) ){
		runToCheck = line.substr(0, line.find(delimiter));
		//cout << runToCheck << endl;
		bool foundRun = isInList( Good_Runs, stoi(runToCheck) );
		if( foundRun ){
			fout << line << endl;
			cout << runToCheck << endl;
		}
	}
	fin2.close();
	fout.close();

}
