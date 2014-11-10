#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream> // read files
#include <sstream>
#include <iomanip> // setprecision
#include <cmath>   // log

#include "readFiles.h"

using namespace std;


int main(int argc, char *argv[]) {
    // read in data:
	vector<vector<double> > sites;
    readSites( "admixSampleData.sites", sites );
    print2Dvec( sites );

    vector<int> locs;
    readLocs( "admixSampleData.locs", locs);
    print1Dvec( locs );

    vector<vector<string> > hapInfo;
    readHapInfo( "admixSampleData.hapnames", hapInfo );
    print2DvecString( hapInfo );

}



