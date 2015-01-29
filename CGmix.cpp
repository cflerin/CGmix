#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream> // read files
#include <sstream>
#include <iomanip> // setprecision
#include <cmath>   // log
#include <algorithm>

#include "readFiles.h"
#include "hmm.h"

using namespace std;

int gMode;

int main(int argc, char *argv[]) {
    std::streamsize ss = std::cout.precision();
    if( argc < 3 ) {
        std::cerr << "Usage: " << argv[0] << " [file prefix, with *.sites, *.locs *.hapnames endings] [0: limit to crossover analysis; 1: full model; 2: Two-pass model]" << std::endl;
        return 1;
    }
    if( strcmp( argv[2] , "1" ) == 0 ) // run full model
        gMode = 1; //atoi( argv[2] );
    else if( strcmp( argv[2] , "0" ) == 0 ) // run in haplotype-only mode
        gMode = 0;
    else if( strcmp( argv[2] , "2" ) == 0 ) // run first pass, then second
        gMode = 2;
    //
    string fname = argv[1];
    string logf = fname + ".log" + std::to_string(gMode);
    ofstream logfile ( logf.c_str() );
    string pathf = fname + ".path" + std::to_string(gMode);
    ofstream pathfile ( pathf.c_str() );
    string matf = fname + ".mat" + std::to_string(gMode);
    ofstream matfile ( matf.c_str() );
    logfile << argv[0] << " " << argv[1] << endl;

    // read in data:
    vector<vector<int> > sites;
    readSites( fname + ".sites", sites );
    logfile << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
    // print2Dvec( sites );

    vector<int> locs;
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
    param.S = locs.size() ;
    param.T = 7;
    param.u1 = 0.5;
    param.rho = 1.0/100000;
    param.gam = 1.0/10000;
    param.lam = 1.0/500;
    param.theta = 1.0/1000;

    hmmStates st;
    if( ( gMode == 0 ) || ( gMode == 2 ) ) {
        generateXstates( hapInfo, st );
        param.theta = numeric_limits<double>::epsilon(); //1.0/100000;
    } else if( gMode == 1 ) {
        generateStates( hapInfo, st );
    }
    logfile << "Generated " << st.states.size() << " states" << endl;
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

    if( gMode == 0 )
        logfile << "Model: haplotype-only" << endl;
    else if( gMode == 1 )
        logfile << "Model: full haplotype and gene conversion" << endl;
    else if( gMode == 2 )
        logfile << "Model: Two pass approach" << endl;

    logfile << "Starting forward algorithm..." << endl;;
    vector<double> sprob( st.states.size(), 0 );
    if( ( gMode == 0 ) || ( gMode == 2 ) ) {
        getsprobX( sites[0], param, emit, st, obs, sprob );
    } else if( gMode == 1 ) {
        getsprob( sites[0], param, emit, st, obs, sprob );
    }
    vector<vector<double> > fwd(param.S, vector<double>(st.states.size(), 0.0));
    forward( sites, locs, param, emit, st, obs, sprob, fwd );
    // printMat( fwd );

    double Pxa, Pxb;
    logfile << "Starting backward algorithm..." << endl;
    vector<vector<double> > bwd(param.S, vector<double>(st.states.size(), 0.0));
    //backward( sites, locs, param, emit, st, obs, bwd );
    backward( sites, locs, param, emit, st, obs, sprob, bwd, Pxb );
    //
    //logSumExp( fwd[ fwd.size()-1 ], Pxa );
    Pxa = 0.0;
    for(int i=0; i<fwd[ fwd.size()-1 ].size(); i++ ) {
        Pxa += exp( fwd[ fwd.size()-1 ][i] );
    }
    Pxa = log( Pxa );
    Pxb = log( Pxb );
    logfile << setprecision(20) << "Pxa= " << Pxa << endl;
    logfile << setprecision(20) << "Pxb= " << Pxb << endl;
    logfile << "diff= " << Pxa - Pxb << endl;
    // printMat( bwd );

    logfile << "Starting posterior decoding..." << endl;
    vector<vector<double> > pprob(param.S, vector<double>(st.states.size(), 0.0));
    vector<string> pppath( sites.size() );
    vector<double> ppprob( sites.size() , 0.0);
    vector<int> pswitch( sites.size(), 0 );
    postDecode( fwd, bwd, st, Pxa, pprob, pppath, ppprob, pswitch, logfile);
    // printMat( pprob );
    // writeMat( pprob, matfile );

    logfile << "Starting Viterbi algorithm..." << endl;
    vector<vector<double> > vit(param.S, vector<double>(st.states.size(), 0.0));
    vector<string> vpath( sites.size() );
    vector<double> vprob( sites.size() , 0.0);
    viterbi( sites, locs, param, emit, st, obs, sprob, vit, vpath, vprob );
    // printMat( vit );

    if( ( gMode == 0 ) || ( gMode == 2 ) ) {
        // find haplotype switches:
        vector<int> swIx;
        for(int j=0; j < pswitch.size(); j++ ) {
            if( pswitch[j] == 1 )
                swIx.push_back( j );
        }
        //for(int j=0; j<pswitch.size(); j++)
        //    cout << pswitch[j] << " ";
        //cout << endl;
        int width = 10; 
        // expand marks around switches:
        vector<int> repl(2);
        for(int j=0; j < swIx.size(); j++) {
            repl[0] = swIx[j] - width;
            repl[1] = swIx[j] + width;
            if( repl[0] < 0 )
                repl[0] = 0;
            if( repl[1] > (param.S-1) )
                repl[1] = (param.S-1);
            for(int k=repl[0]; k <= repl[1]; k++ )
                pswitch[k] = 1;
        }
    }

    if( ( gMode == 0 ) || ( gMode == 1 ) ) {
        // output Viterbi path and probabilites:
        vector<double> gcprob( sites.size(), 0.0 );
        double csum;
        pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\tGCprob" << endl;
        for(int j=0; j < pprob.size(); j++ ) {
            csum =0.0;
            if( gMode == 1 ) {
                for(int i=0; i < st.states.size(); i++ ) {
                    if( st.Ghap[i] == 0 )
                        continue;
                    csum += pprob[j][i];
                }
            }
            gcprob[j] = 1.0 - csum;
            pathfile << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << exp(vprob[j]) << "\t" << pppath[j] << "\t" << ppprob[j] << "\t" << pswitch[j] << "\t" << 1-gcprob[j] << endl;
        }

        matfile << "forward" << endl;
        writeTmat( fwd, matfile );
        matfile << "backward" << endl;
        writeTmat( bwd, matfile );
        matfile << "posterior" << endl;
        writeTmat( pprob, matfile );
    }

    ////////////////
    //std::fill( pswitch.begin(), pswitch.end(), 0 );
    //for(int j=2; j < 3; j++) {
    //    pswitch[j] = 1;
    //}
    // for(int j=0; j<pswitch.size(); j++ ) {
    //     cout << pswitch[j] << endl;
    // }
    ////////////////

    // restart hmm:
    if( gMode == 2 ) {
        logfile.precision (ss);
        logfile << "Restart hmm" << endl;
        param.theta = 1.0/1000;
        hmmStates st2;
        generateStates( hapInfo, st2 );
        logfile << "Generated " << st2.states.size() << " states" << endl;
        logfile << "Parameters set:" << endl;
        logfile << "n1 = " << param.n1 << endl;
        logfile << "n2 = " << param.n2 << endl;
        logfile << "nSites = " << param.S << endl;
        logfile << "u1 = " << param.u1 << endl;
        logfile << "rho = " << param.rho << endl;
        logfile << "gamma = " << param.gam << endl;
        logfile << "lambda = " << param.lam << endl;
        logfile << "theta = " << param.theta << endl;

        //logfile << "new emissions" << endl;
        // emissions probabilities:
        emissions emit;
        emit.match = (2.0 * (param.n1 + param.n2) + param.theta ) / ( 2.0 * ( (param.n1 + param.n2) + param.theta) );
        emit.mismatch = param.theta / ( 2.0 * ( ( param.n1 + param.n2 ) + param.theta ) );
        emit.match = log( emit.match );
        emit.mismatch = log( emit.mismatch );

        //logfile << "new sprob" << endl;
        if( pswitch[0] == 1 ) {
            sprob.resize( st2.states.size(), 0.0 );
            std::fill( sprob.begin(), sprob.end(), 0.0 );
            getsprob( sites[0], param, emit, st2, obs, sprob );
        } else {
            sprob.resize( st.states.size(), 0.0 );
            std::fill( sprob.begin(), sprob.end(), 0.0 );
            getsprobX( sites[0], param, emit, st, obs, sprob );
        }

        for(int j=0; j < param.S; j++ ) {
            if( pswitch[j] == 0 ) {
                fwd[j].resize( st.states.size(), 0.0 );
                std::fill( fwd[j].begin(), fwd[j].end(), 0.0 );
                bwd[j].resize( st.states.size(), 0.0 );
                std::fill( bwd[j].begin(), bwd[j].end(), 0.0 );
                pprob[j].resize( st.states.size(), 0.0 );
                std::fill( pprob[j].begin(), pprob[j].end(), 0.0 );
                vit[j].resize( st.states.size(), 0.0 );
                std::fill( vit[j].begin(), vit[j].end(), 0.0 );
            } else {
                fwd[j].resize( st2.states.size(), 0.0 );
                std::fill( fwd[j].begin(), fwd[j].end(), 0.0 );
                bwd[j].resize( st2.states.size(), 0.0 );
                std::fill( bwd[j].begin(), bwd[j].end(), 0.0 );
                pprob[j].resize( st2.states.size(), 0.0 );
                std::fill( pprob[j].begin(), pprob[j].end(), 0.0 );
                vit[j].resize( st2.states.size(), 0.0 );
                std::fill( vit[j].begin(), vit[j].end(), 0.0 );
            }
            // cout << "j = " << j << "\t" << fwd[j].size() << endl;
        }

        logfile << "Starting forward algorithm..." << endl;;
        forward2( sites, locs, param, emit, st, st2, obs, sprob, pswitch, fwd );
        matfile << "forward" << endl;
        writeTmat( fwd, matfile );


        logfile << "Starting backward algorithm..." << endl;
        backward2( sites, locs, param, emit, st, st2, obs, sprob, pswitch, bwd, Pxb );
        matfile << "backward" << endl;
        writeTmat( bwd, matfile );

        Pxa = 0.0;
        for(int i=0; i<fwd[ fwd.size()-1 ].size(); i++ ) {
            Pxa += exp( fwd[ fwd.size()-1 ][i] );
        }
        Pxa = log( Pxa );
        Pxb = log( Pxb );
        logfile << setprecision(20) << "Pxa= " << Pxa << endl;
        logfile << setprecision(20) << "Pxb= " << Pxb << endl;
        logfile << "diff= " << Pxa - Pxb << endl;

        logfile << "Starting posterior decoding..." << endl;
        vector<string> pppath2( sites.size() );
        vector<double> ppprob2( sites.size() , 0.0);
        vector<int> pswitch2 = pswitch;
        postDecode( fwd, bwd, st2, Pxa, pprob, pppath2, ppprob2, pswitch2, logfile);
        matfile << "posterior" << endl;
        writeTmat( pprob, matfile );

        logfile << "Starting Viterbi algorithm..." << endl;
        vector<string> vpath2( sites.size() );
        vector<double> vprob2( sites.size() , 0.0);
        viterbi2( sites, locs, param, emit, st2, obs, sprob, pswitch, vit, vpath2, vprob2 );

        // output Viterbi path and probabilites:
        vector<double> gcprob2( sites.size(), 0.0 );
        double csum;
        pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\tVpath2\tVpathProb2\tPpath2\tPpathProb2\tGCprob" << endl;
        for(int j=0; j < pprob.size(); j++ ) {
            csum =0.0;
            for(int i=0; i < pprob[j].size(); i++ ) {
                if( st2.Ghap[i] == 0 )
                    continue;
                csum += pprob[j][i];
            }
            gcprob2[j] = 1.0 - csum;
            pathfile << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << exp(vprob[j]) << "\t" << pppath[j] << "\t" << ppprob[j] << "\t" << pswitch[j] << "\t";
            pathfile << vpath2[j] << "\t" << exp(vprob2[j]) << "\t" << pppath2[j] << "\t" << ppprob2[j] << "\t" << 1-gcprob2[j] << endl;
        }
    } // end of 2nd pass

    logfile << "Done!" << endl;
    logfile.close();
    pathfile.close();
    matfile.close();

}



