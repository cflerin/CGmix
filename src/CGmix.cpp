/*
 * CGmix.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: ccampbell
 */

#include "CGmix.h"

int main(int argc, char *argv[]) {
    std::streamsize ss = std::cout.precision();
    time_t start,end;
    time(&start);
    parameters param(argc, argv);
//    params.print_help();
    param.read_parameters();

    if( param.fname == "empty" ) {
        cout << "Input file not set" << endl;
        exit(1);
    }
    if( param.mode == 99 ) {
        cout << "Mode must be 0, 1 or 2" << endl;
        exit(1);
    }

    ofstream logfile ( param.logf.c_str() );
    ofstream pathfile ( param.pathf.c_str() );
    //ofstream matfile;

    param.print_params(logfile, 0);
    logfile << endl;

    // read in data:
    vector<vector<int> > sites;
    readSites( param.fname + ".sites", sites );
    logfile << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
    // print2Dvec( sites );

    vector<double> locs;
    readLocs( param.fname + ".locs", locs);
    logfile << "Read " << locs.size() << " physical positions" << endl;
    // print1Dvec( locs ); cout << endl;

    vector<vector<string> > tmp;
    hapDef hapInfo;
    readHapInfo( param.fname + ".hapnames", tmp );
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

    geneticMap gMap;
    readGenMap( param.gmfile, gMap );
    logfile << "Read " << gMap.chr.size() << " lines from genetic map" << endl;

    //interpolate values from locs:
    positions pos;
    interpGenMap( gMap, locs, pos );
    logfile << "Interpolated " << pos.cM.size() << " positions from genetic map" << endl;
    
    // set kb throughout:
    for(int i=0; i<pos.pos.size(); i++) {
        //cout << setprecision(20) << pos.pos[i] << "\t" << pos.cM[i] << endl;
        pos.pos[i] = pos.pos[i] / 1000.0;   // kb
        //pos.rate[i] = pos.rate[i] / 1000.0; // cM/kb
        pos.cM[i] = pos.cM[i] / 100.0;     // Morgans
        //cout << setprecision(20) << pos.pos[i] << "\t" << pos.cM[i] << endl;
    }

    // set/get/update parameters:
    for(int i=0; i < hapInfo.hapPop.size(); i++ ) {
        if( hapInfo.hapPop[i] == "p1" )
            param.n1 += 1;
        if( hapInfo.hapPop[i] == "p2" )
            param.n2 += 1;
    }
    param.S = locs.size() ;
    param.u2 = 1.0 - param.u1;

    hmmStates st;
    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        generateXstates( hapInfo, st );
        // set theta to ~0 for reduced first pass:
        param.theta1 = numeric_limits<double>::epsilon();
        param.theta2 = numeric_limits<double>::epsilon();
    } else if( param.mode == 1 ) {
        generateStates( hapInfo, st );
        // set theta to a reasonable number for full model:
        param.theta1 = 0.2 / ( 0.2 + param.n1 );
        param.theta2 = 0.2 / ( 0.2 + param.n2 );
    }
    logfile << "Generated " << st.states.size() << " states" << endl;

    // emissions probabilities:
    param.theta1_match = log( 1.0 - param.theta1 );
    param.theta1_mismatch = log( param.theta1 );
    param.theta2_match = log( 1.0 - param.theta2 );
    param.theta2_mismatch = log( param.theta2 );
    //cout << "theta1 " << exp( param.theta1_match ) << "\t" << exp( param.theta1_mismatch ) << endl;
    //cout << "theta2 " << exp( param.theta2_match ) << "\t" << exp( param.theta2_mismatch ) << endl;

    param.print_params(logfile, 1);

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
    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        getsprobX( sites[0], param, st, obs, sprob );
        param.passAcc = 1; // force LSE sort on first pass
    } else if( param.mode == 1 ) {
        getsprob( sites[0], param, st, obs, sprob );
        param.passAcc = param.highAccuracy; // set to user-selected state
    }
    vector<vector<double> > fwd(param.S, vector<double>(st.states.size(), 0.0));
    forward( sites, pos, param, st, obs, sprob, fwd );

    double Pxa, Pxb;
    logfile << "Starting backward algorithm..." << endl;
    vector<vector<double> > bwd(param.S, vector<double>(st.states.size(), 0.0));
    backward( sites, pos, param, st, obs, sprob, bwd, Pxb );
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
    //if( ( Pxa - Pxb ) > ( 4 * numeric_limits<double>::epsilon() ) ) {
        cout << "P(x) diff=" << Pxa - Pxb << endl;
    //}

    pathVec pvec( param );

    logfile << "Starting posterior decoding..." << endl;
    vector<vector<double> > pprob(param.S, vector<double>(st.states.size(), 0.0));
    //postDecode( fwd, bwd, st, Pxa, pprob, pvec, logfile);
    postDecode( fwd, bwd, st, Pxa, pprob, pvec.pppath, pvec.ppprob, logfile);

    vector<vector<double> > vit;
    if( param.viterbi == 1 ) {
        logfile << "Starting Viterbi algorithm..." << endl;
        vector<vector<double> > vit(param.S, vector<double>(st.states.size(), 0.0));
        viterbi( sites, pos, param, st, obs, sprob, vit, pvec );
        // printMat( vit );
    } else {
        logfile << "Skipping Viterbi algorithm..." << endl;
    }

    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        // find haplotype switches:
        for(int j=1; j < param.S; j++ ) {
            if( pvec.pppath[j] != pvec.pppath[j-1] )
                pvec.pswitch[j] = 1;
        }
        vector<int> swIx;
        for(int j=0; j < param.S; j++ ) {
            if( pvec.pswitch[j] == 1 )
                swIx.push_back( j );
        }
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
                pvec.pswitch[k] = 1;
        }
    }

    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        logfile << "Starting path output" << endl;
        pathOutput( pvec, st, pos, pprob, param, pathfile );
    }

    ////////////////
    if( param.fixPswitch > 0 ) { // testing only
        std::fill( pvec.pswitch.begin(), pvec.pswitch.end(), 0 );
        for(int j=0; j < param.fixPswitch; j++) {
            pvec.pswitch[j] = 1;
        }
    }
    ////////////////

    // restart hmm:
    if( param.mode == 2 ) {
        logfile.precision (ss);
        logfile << endl << "Restart hmm" << endl;
        // set theta to for full model:
        param.theta1 = 0.2 / ( 0.2 + param.n1 );
        param.theta2 = 0.2 / ( 0.2 + param.n2 );

        hmmStates st2;
        generateStates( hapInfo, st2 );
        logfile << "Generated " << st2.states.size() << " states" << endl;
        param.print_params(logfile, 1);

        param.passAcc = param.highAccuracy; // set back to user selected state

        if( pvec.pswitch[0] == 1 ) {
            sprob.resize( st2.states.size(), 0.0 );
            std::fill( sprob.begin(), sprob.end(), 0.0 );
            getsprob( sites[0], param, st2, obs, sprob );
        } else {
            sprob.resize( st.states.size(), 0.0 );
            std::fill( sprob.begin(), sprob.end(), 0.0 );
            getsprobX( sites[0], param, st, obs, sprob );
        }

        for(int j=0; j < param.S; j++ ) {
            if( pvec.pswitch[j] == 0 ) {
                fwd[j].resize( st.states.size(), 0.0 );
                std::fill( fwd[j].begin(), fwd[j].end(), 0.0 );
                bwd[j].resize( st.states.size(), 0.0 );
                std::fill( bwd[j].begin(), bwd[j].end(), 0.0 );
                pprob[j].resize( st.states.size(), 0.0 );
                std::fill( pprob[j].begin(), pprob[j].end(), 0.0 );
                if( param.viterbi == 1 ) {
                    vit[j].resize( st.states.size(), 0.0 );
                    std::fill( vit[j].begin(), vit[j].end(), 0.0 );
                }
            } else {
                fwd[j].resize( st2.states.size(), 0.0 );
                std::fill( fwd[j].begin(), fwd[j].end(), 0.0 );
                bwd[j].resize( st2.states.size(), 0.0 );
                std::fill( bwd[j].begin(), bwd[j].end(), 0.0 );
                pprob[j].resize( st2.states.size(), 0.0 );
                std::fill( pprob[j].begin(), pprob[j].end(), 0.0 );
                if( param.viterbi == 1 ) {
                    vit[j].resize( st2.states.size(), 0.0 );
                    std::fill( vit[j].begin(), vit[j].end(), 0.0 );
                }
            }
            // cout << "j = " << j << "\t" << fwd[j].size() << endl;
        }

        logfile << "Starting forward algorithm..." << endl;;
        forward2( sites, pos, param, st, st2, obs, sprob, pvec.pswitch, fwd );

        logfile << "Starting backward algorithm..." << endl;
        backward2( sites, pos, param, st, st2, obs, sprob, pvec.pswitch, bwd, Pxb );

        Pxa = 0.0;
        for(int i=0; i<fwd[ fwd.size()-1 ].size(); i++ ) {
            Pxa += exp( fwd[ fwd.size()-1 ][i] );
        }
        Pxa = log( Pxa );
        Pxb = log( Pxb );
        logfile << setprecision(20) << "Pxa= " << Pxa << endl;
        logfile << setprecision(20) << "Pxb= " << Pxb << endl;
        logfile << "diff= " << Pxa - Pxb << endl;

        //if( ( Pxa - Pxb ) > ( 4 * numeric_limits<double>::epsilon() ) ) {
            cout << "P(x) diff=" << Pxa - Pxb << endl;
        //}

        logfile << "Starting posterior decoding..." << endl;
        pvec.pppath2.resize( sites.size() );
        pvec.ppprob2.resize( sites.size(), 0.0 );
        pvec.pswitch2 = pvec.pswitch;
        //postDecode( fwd, bwd, st2, Pxa, pprob, pvec, logfile);
        postDecode( fwd, bwd, st2, Pxa, pprob, pvec.pppath2, pvec.ppprob2, logfile);

        if( param.viterbi == 1 ) {
            logfile << "Starting Viterbi algorithm..." << endl;
            pvec.vpath2.resize( sites.size() );
            pvec.vprob2.resize( sites.size(), 0.0 );
            viterbi2( sites, pos, param, st2, obs, sprob, vit, pvec );
        } else {
            logfile << "Skipping Viterbi algorithm..." << endl;
        }

        logfile << "Starting path output" << endl;
        pvec.gcprob.resize( sites.size(), 0.0 );
        std::fill( pvec.gcprob.begin(), pvec.gcprob.end(), 0.0 );
        pvec.gcprobPop1.resize( sites.size(), 0.0 );
        std::fill( pvec.gcprobPop1.begin(), pvec.gcprobPop1.end(), 0.0 );
        pvec.gcprobPop2.resize( sites.size(), 0.0 );
        std::fill( pvec.gcprobPop2.begin(), pvec.gcprobPop2.end(), 0.0 );
        std::fill( pvec.transPGC.begin(), pvec.transPGC.end(), 0.0 );
        pathOutput( pvec, st2, pos, pprob, param, pathfile );

    } // end of 2nd pass

    time(&end);
    double running_time = difftime(end,start);
    //logfile("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
    logfile << "Run Time = " << running_time << endl;

    logfile.close();
    pathfile.close();
    if( param.matrixOutput == 1 ) {
        ofstream matfile ( param.matf.c_str() );
        matfile << "forward" << endl;
        writeTmat( fwd, matfile );
        matfile << "backward" << endl;
        writeTmat( bwd, matfile );
        matfile << "posterior" << endl;
        writeTmat( pprob, matfile );
        matfile.close();
    }
}



