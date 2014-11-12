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
#include "hmm.h"

using namespace std;

int main(int argc, char *argv[]) {
    // read in data:
    vector<vector<int> > sites;
    readSites( "admixSampleData.sites", sites );
    cout << "Read " << sites.size() << " haplotypes with " << sites[0].size() << " sites:" << endl;
    print2Dvec( sites );

    vector<int> locs;
    readLocs( "admixSampleData.locs", locs);
    cout << endl << "Read " << locs.size() << " physical positions:" << endl;
    print1Dvec( locs );
    cout << endl;

    vector<vector<string> > hapInfo;
    readHapInfo( "admixSampleData.hapnames", hapInfo );
    cout << endl << "Read " << hapInfo.size() << " haplotype definitions:" << endl;
    print2DvecString( hapInfo );

    hmmStates st;
    generateStates( hapInfo, st );
    cout << endl << "Generated " << st.states.size() << " states" << endl << endl;

    // set/get parameters:
    parameters param;
    param.n1 = 5; // get from sites
    param.n2 = 5; // get from sites
    param.S = 20; // get from sites
    param.T = 7;
    param.u1 = 0.8;
    param.rho = 1.0/1000000;
    param.gam = 1.0/10000;
    param.lam = 1.0/500;
    param.theta = 1.0/1000;
    // emissions probabilities:
    emissions emit;
    emit.match = (2.0 * (param.n1 + param.n2) + param.theta ) / ( 2.0 * ( (param.n1 + param.n2) + param.theta) );
    emit.mismatch = param.theta / ( 2.0 * ( ( param.n1 + param.n2 ) + param.theta ) );

    vector<int> obs;
    for(int i=0; i < sites.size(); i++) {
        if( hapInfo[i][0] == "obs" ) {
            for(int j=0; j < param.S; j++) {
                // cout << sites[i][j] << endl;
                obs.push_back( sites[i][j] );
            }
        }
    }

    vector<vector<double> > fwd(st.states.size(), vector<double>(param.S, 0));
    forward( sites, locs, param, emit, st, obs, fwd );

}



