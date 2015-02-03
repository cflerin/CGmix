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


    ofstream logfile ( param.logf.c_str() );
    ofstream pathfile ( param.pathf.c_str() );
    ofstream matfile ( param.matf.c_str() );

    param.print_params(logfile, 0);
    logfile << endl;


    // read in data:
    vector<vector<int> > sites;
    readSites( param.fname + ".sites", sites );
    logfile << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
    // print2Dvec( sites );

    vector<int> locs;
    readLocs( param.fname + ".locs", locs);
    logfile << "Read " << locs.size() << " physical positions" << endl;
    // print1Dvec( locs ); cout << endl;

    vector<vector<string> > tmp;
    class hapDef hapInfo;
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

    // set/get/update parameters:
    for(int i=0; i < hapInfo.hapPop.size(); i++ ) {
        if( hapInfo.hapPop[i] == "p1" )
            param.n1 += 1;
        if( hapInfo.hapPop[i] == "p2" )
            param.n2 += 1;
    }
    param.S = locs.size() ;
    param.theta1 = 0.2 / ( 0.2 / param.n1 );
    param.theta2 = 0.2 / ( 0.2 / param.n2 );
    param.rho1 = param.rho1 / param.n1;
    param.rho2 = param.rho2 / param.n2;
    param.gam1 = param.gam1 / param.n1;
    param.gam2 = param.gam2 / param.n2;
    //////
    param.rho = 1.0/100000;
    param.gam = 1.0/10000;
    param.lam = 1.0/500;
    param.theta = 1.0/1000;
    //////

    hmmStates st;
    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        generateXstates( hapInfo, st );
        param.theta = numeric_limits<double>::epsilon(); //1.0/100000;
    } else if( param.mode == 1 ) {
        generateStates( hapInfo, st );
    }
    logfile << "Generated " << st.states.size() << " states" << endl;

    param.print_params(logfile, 1);

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
    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        getsprobX( sites[0], param, emit, st, obs, sprob );
    } else if( param.mode == 1 ) {
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

    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        // find haplotype switches:
        vector<int> swIx;
        for(int j=0; j < pswitch.size(); j++ ) {
            if( pswitch[j] == 1 )
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
                pswitch[k] = 1;
        }
    }

    logfile << "Starting path output" << endl;
    vector<double> gcprob( sites.size(), 0.0 );
    vector<double> gcprobPop1( sites.size(), 0.0 );
    vector<double> gcprobPop2( sites.size(), 0.0 );
    if( ( param.mode == 0 ) || ( param.mode == 1 ) ) {
        // output Viterbi path and probabilites:
        double csum;
        pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\tGCprob\tGCprobP1\tGCprobP2" << endl;
        for(int j=0; j < pprob.size(); j++ ) {
            for(int i=0; i < pprob[j].size(); i++ ) {
                if( st.Ghap[i] != 0 ) { // gene conversion state
                    gcprob[j] += pprob[j][i];
                    if( st.Gpop[i] == 1 ) {
                        gcprobPop1[j] += pprob[j][i];
                    } else if( st.Gpop[i] == 2 ) {
                        gcprobPop2[j] += pprob[j][i];
                    }
                }
            }
            pathfile << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << exp(vprob[j]) << "\t" << pppath[j] << "\t" << ppprob[j] << "\t" << pswitch[j] << "\t" << gcprob[j] << "\t" << gcprobPop1[j] << "\t" << gcprobPop2[j] << endl;
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
    if( param.mode == 2 ) {
        logfile.precision (ss);
        logfile << endl << "Restart hmm" << endl;
        param.theta = 1.0/1000;
        hmmStates st2;
        generateStates( hapInfo, st2 );
        logfile << "Generated " << st2.states.size() << " states" << endl;
        param.print_params(logfile, 1);

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
        gcprob.resize( sites.size(), 0.0 );
        std::fill( gcprob.begin(), gcprob.end(), 0.0 );
        gcprobPop1.resize( sites.size(), 0.0 );
        std::fill( gcprobPop1.begin(), gcprobPop1.end(), 0.0 );
        gcprobPop2.resize( sites.size(), 0.0 );
        std::fill( gcprobPop2.begin(), gcprobPop2.end(), 0.0 );
        pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\tVpath2\tVpathProb2\tPpath2\tPpathProb2\tGCprob\tGCprobP1\tGCprobP2" << endl;
        for(int j=0; j < pprob.size(); j++ ) {
            for(int i=0; i < pprob[j].size(); i++ ) {
                if( st2.Ghap[i] != 0 ) { // gene conversion state
                    gcprob[j] += pprob[j][i];
                    if( st2.Gpop[i] == 1 ) {
                        gcprobPop1[j] += pprob[j][i];
                    } else if( st2.Gpop[i] == 2 ) {
                        gcprobPop2[j] += pprob[j][i];
                    }
                }
            }
            pathfile << j << "\t" << locs[j] << "\t" << vpath[j] << "\t" << exp(vprob[j]) << "\t" << pppath[j] << "\t" << ppprob[j] << "\t" << pswitch[j] << "\t";
            pathfile << vpath2[j] << "\t" << exp(vprob2[j]) << "\t" << pppath2[j] << "\t" << ppprob2[j] << "\t" << gcprob[j] << "\t" << gcprobPop1[j] << "\t" << gcprobPop2[j] << endl;
        }
    } // end of 2nd pass

    time(&end);
    double running_time = difftime(end,start);
    //logfile("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
    logfile << "Run Time = " << running_time << endl;

    logfile.close();
    pathfile.close();
    matfile.close();
}



