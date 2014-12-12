#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream> // read files
#include <sstream>
#include <iomanip> // setprecision
#include <cmath>   // log
#include <algorithm>

#include "readFiles.h"
#include "hmm.h"

using namespace std;

int main(int argc, char *argv[]) {
    if( argc < 2 ) {
        //std::cerr << "Usage: " << argv[0] << " sites locs hapnames" << std::endl;
        std::cerr << "Usage: " << argv[0] << " [file prefix, with *.sites, *.locs *.hapnames endings]" << std::endl;
        return 1;
    }
    //cout << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << endl;
    string fname = argv[1];
    string logf = fname + ".log";
    ofstream logfile ( logf.c_str() );
    string matf = fname + ".mat";
    ofstream matfile ( matf.c_str() );
    logfile << argv[0] << " " << argv[1] << endl;
    // read in data:
    vector<vector<int> > sites;
    readSites( fname + ".sites", sites );
    logfile << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
    // print2Dvec( sites );

    vector<int> locs;
    // readLocs( "admixSampleData.locs", locs);
    readLocs( fname + ".locs", locs);
    logfile << "Read " << locs.size() << " physical positions" << endl;
    // print1Dvec( locs ); cout << endl;

    vector<vector<string> > tmp;
    class hapDef hapInfo;
    readHapInfo( fname + ".hapnames", tmp );
    for(int i=0; i < tmp.size(); i++ ) {
        hapInfo.hapName.push_back( tmp[i][0] );
        hapInfo.hapPop.push_back( tmp[i][1] );
        hapInfo.hN.push_back(i+1);
        if( tmp[i][1] == "p1" )
            hapInfo.hP.push_back(1);
        if( tmp[i][1] == "p2" )
            hapInfo.hP.push_back(2);
        if( tmp[i][1] == "p3" )
            hapInfo.hP.push_back(3);
        // cout << hapInfo.hapName[i] << " " << hapInfo.hapPop[i] << " " << hapInfo.hN[i] << " " << hapInfo.hP[i] << endl;
    }
    logfile << "Read " << hapInfo.hN.size() << " haplotype definitions" << endl;
    hmmStates st;
    generateStates( hapInfo, st );
    logfile << "Generated " << st.states.size() << " states" << endl;

    // set/get parameters:
    parameters param;
    param.n1 = 0;
    param.n2 = 0;
    for(int i=0; i < hapInfo.hapPop.size(); i++ ) {
        if( hapInfo.hapPop[i] == "p1" )
            param.n1 += 1;
        if( hapInfo.hapPop[i] == "p2" )
            param.n2 += 1;
    }
    // cout << "n1= " << param.n1 << " | n2= " << param.n2 << endl;
    param.S = locs.size() ;
    param.T = 7;
    param.u1 = 0.8;
    param.rho = 1.0/1000000;
    param.gam = 1.0/100000;
    //param.gam = 1.0/10000;
    param.lam = 1.0/500;
    param.theta = 1.0/1000;
    logfile << "Parameters set:" << endl;
    logfile << "n1 = " << param.n1 << endl;
    logfile << "n2 = " << param.n2 << endl;
    logfile << "nSites = " << param.S << endl;
    logfile << "u1 = " << param.u1 << endl;
    logfile << "rho = " << param.rho << endl;
    logfile << "gamma = " << param.gam << endl;
    logfile << "lambda = " << param.lam << endl;
    logfile << "theta = " << param.theta << endl;
    // emissions probabilities:
    emissions emit;
    emit.match = (2.0 * (param.n1 + param.n2) + param.theta ) / ( 2.0 * ( (param.n1 + param.n2) + param.theta) );
    emit.mismatch = param.theta / ( 2.0 * ( ( param.n1 + param.n2 ) + param.theta ) );
    emit.match = log( emit.match );
    emit.mismatch = log( emit.mismatch );

    vector<int> obs( param.S, 0 );
    for(int i=0; i < hapInfo.hN.size(); i++) {
        if( hapInfo.hapName[i] == "obs" ) {
            for(int j=0; j < sites.size(); j++) {
                obs[j] = sites[j][i];
            }
        }
    }

    logfile << "Starting forward algorithm..." << endl;;
    vector<double> sprob( st.states.size(), 0 );
    getsprob( sites[0], param, emit, st, obs, sprob );
    vector<vector<double> > fwd(param.S, vector<double>(st.states.size(), 0.0));
    forward( sites, locs, param, emit, st, obs, sprob, fwd );
    logfile << "finished" << endl;
    //printMat( fwd );

    logfile << "Starting backward algorithm..." << endl;
    vector<vector<double> > bwd(param.S, vector<double>(st.states.size(), 0.0));
    backward( sites, locs, param, emit, st, obs, bwd );
    logfile << "finished" << endl;
    // printMat( bwd );

    logfile << "Starting posterior decoding..." << endl;
    vector<vector<double> > pprob(param.S, vector<double>(st.states.size(), 0.0));
    vector<string> pppath( sites.size() );
    vector<double> ppprob( sites.size() , 0.0);
    postDecode( fwd, bwd, st, pprob, pppath, ppprob, logfile);
    logfile << "finished" << endl;
    // printMat( pprob );
    writeMat( pprob, matfile );

    logfile << "Starting Viterbi algorithm..." << endl;
    vector<vector<double> > vit(param.S, vector<double>(st.states.size(), 0.0));
    vector<string> vpath( sites.size() );
    vector<double> vprob( sites.size() , 0.0);
    // viterbi( sites, locs, param, emit, st, obs, pprob, sprob, vit, vpath, vprob );
    viterbi( sites, locs, param, emit, st, obs, sprob, vit, vpath, vprob );
    // printMat( vit );
    logfile << "finished" << endl;

    // output Viterbi path and probabilites:
    vector<double> gcprob( sites.size(), 0.0 );
    double csum;
    logfile << "site\tppos\tVpath\tVpathProb\tPpath\tPpathProb\tGCProb" << endl;
    for(int j=0; j < pprob.size(); j++ ) {
        csum =0.0;
        for(int i=0; i < st.states.size(); i++ ) {
            if( st.Ghap[i] == 0 ) 
                continue;
            csum += pprob[j][i];
        }
        gcprob[j] = 1.0 - csum;
        // cout << j << "\t" << 1-gcprob[j] << " " << endl;
        logfile << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << exp(vprob[j]) << "\t" << pppath[j] << "\t" << ppprob[j] << "\t" << 1-gcprob[j] << endl;
    }

    logfile.close();
}



