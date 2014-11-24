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
    cout << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
    // print2Dvec( sites );

    vector<int> locs;
    readLocs( "admixSampleData.locs", locs);
    cout << "Read " << locs.size() << " physical positions" << endl;
    // print1Dvec( locs ); cout << endl;

    vector<vector<string> > hapInfo;
    readHapInfo( "admixSampleData.hapnames", hapInfo );
    cout << "Read " << hapInfo.size() << " haplotype definitions" << endl;
    // print2DvecString( hapInfo );

    hmmStates st;
    generateStates( hapInfo, st );
    cout << "Generated " << st.states.size() << " states" << endl << endl;

    // set/get parameters:
    parameters param;
    param.n1 = 0;
    param.n2 = 0;
    for(int i=0; i < hapInfo.size(); i++ ) {
        if( hapInfo[i][1] == "p1" )
            param.n1 += 1;
        if( hapInfo[i][1] == "p2" )
            param.n2 += 1;
    }
    // cout << "n1= " << param.n1 << " | n2= " << param.n2 << endl;
    param.S = locs.size() ;
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
    cout << "n1=" << param.n1 << " n2=" << param.n2 << " S=" << param.S << endl;
    emit.match = log( emit.match );
    emit.mismatch = log( emit.mismatch );

    vector<int> obs;
    for(int i=0; i < hapInfo.size(); i++) {
        if( hapInfo[i][0] == "obs" ) {
            for(int j=0; j < sites.size(); j++) {
                //cout << sites[j][i] ;
                obs.push_back( sites[j][i] );
            }
        }
    }

    cout << "Starting forward algorithm...";
    vector<double> sprob( st.states.size(), 0 );
    getsprob( sites[0], param, emit, st, obs, sprob );

    vector<vector<double> > fwd(param.S, vector<double>(st.states.size(), 0.0));
    forward( sites, locs, param, emit, st, obs, sprob, fwd );
    cout << "finished" << endl;
    // printMat( fwd );

    cout << "Starting backward algorithm...";
    vector<vector<double> > bwd(param.S, vector<double>(st.states.size(), 0.0));
    backward( sites, locs, param, emit, st, obs, bwd );
    cout << "finished" << endl;
    // printMat( bwd );

    cout << "Starting posterior decoding..." << endl;
    vector<vector<double> > pprob(param.S, vector<double>(st.states.size(), 0.0));
    postDecode( fwd, bwd, pprob);
    cout << "finished" << endl;
    // printMat( pprob );

    cout << "Starting Viterbi algorithm...";
    vector<vector<double> > vit(param.S, vector<double>(st.states.size(), 0.0));
    vector<string> vpath( sites.size() );
    vector<double> vprob( sites.size() , 0.0);
    viterbi( sites, locs, param, emit, st, obs, pprob, sprob, vit, vpath, vprob );
    // printMat( vit );
    cout << "finished" << endl;

    // output Viterbi path and probabilites:
    vector<double> gcprob( sites.size(), 0.0 );
    double csum;
    cout << "site\tppos\tpath\tpathProb\tGCProb" << endl;
    for(int j=0; j < pprob.size(); j++ ) {
        csum =0.0;
        for(int i=0; i < st.states.size(); i++ ) {
            if( st.Ghap[i] == "0" ) 
                continue;
            csum += pprob[j][i];
        }
        gcprob[j] = 1.0 - csum;
        // cout << j << "\t" << 1-gcprob[j] << " " << endl;
        cout << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << vprob[j] << "\t" << 1-gcprob[j] << endl;
    }


}



