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
    param.print_help();
    param.read_parameters();

    if( param.fname == "empty"  ) {
        cerr << "Input file not set" << endl;
        exit(1);
    }
    if( param.mode == 99 ) {
        cerr << "Mode must be 0, 1, 2, or 3" << endl;
        exit(1);
    }

    ofstream logfile ( param.logf.c_str() );
    ofstream pathfile ( param.pathf.c_str() );

    param.print_params(logfile, 0);
    logfile << endl;

    ////////////////////////////////////
    // read in data:
    vector<vector<int> > sites;
    vector<double> locs;
    hapDef hapInfo;
    vector<int> obs;
    if( param.fname != "unset" ) {
        readSites( param.fname + ".sites", sites );
        logfile << "Read " << sites[0].size() << " haplotypes with " << sites.size() << " sites" << endl;
        // print2Dvec( sites );

        readLocs( param.fname + ".locs", locs);
        logfile << "Read " << locs.size() << " physical positions" << endl;

        vector<vector<string> > tmp;
        readHapInfo( param.fname + ".hapnames", tmp );
        logfile << "Read " << tmp.size() << " haplotype definitions" << endl;
        vector<int> obsIndx;
        for(int i=0; i < tmp.size(); i++ ) {
            if( tmp[i][1] == "p3" ) {
                //hapInfo.hP.push_back(3);
                obsIndx.push_back(i);
                continue;
            }
            hapInfo.hapName.push_back( tmp[i][0] );
            hapInfo.hapPop.push_back( tmp[i][1] );
            hapInfo.hN.push_back(i+1);
            if( tmp[i][1] == "p1" )
                hapInfo.hP.push_back(1);
            if( tmp[i][1] == "p2" )
                hapInfo.hP.push_back(2);
            // cout << hapInfo.hapName[i] << " " << hapInfo.hapPop[i] << " " << hapInfo.hN[i] << " " << hapInfo.hP[i] << endl;
        }
        //logfile << "Read " << hapInfo.hN.size() << " haplotype definitions" << endl;

        if( obsIndx.size() > 1 ) {
            cerr << "Too many p3 haplotypes! Only one allowed per run." << endl;
            exit(1);
        }
        obs.resize( sites.size(), 0.0 );
        for(int j=0; j < sites.size(); j++) {
            obs[j] = sites[j][ obsIndx[0] ];
            sites[j].erase(sites[j].begin()+ obsIndx[0] );
        }
    } else {
        // reference input:
        readTSites( param.ref + ".sites", sites );
        logfile << "Read " << sites[0].size() << " reference haplotypes with " << sites.size() << " sites" << endl;
        // positions:
        readLocs( param.ref + ".locs", locs);
        logfile << "Read " << locs.size() << " reference physical positions" << endl;
        // hapnames
        vector<vector<string> > tmp;
        readHapInfo( param.ref + ".hapnames", tmp );
        for(int i=0; i < tmp.size(); i++ ) {
            hapInfo.hapName.push_back( tmp[i][0] );
            hapInfo.hapPop.push_back( tmp[i][1] );
            hapInfo.hN.push_back(i+1);
            if( tmp[i][1] == "p1" )
                hapInfo.hP.push_back(1);
            if( tmp[i][1] == "p2" )
                hapInfo.hP.push_back(2);
            // cout << hapInfo.hapName[i] << " " << hapInfo.hapPop[i] << " " << hapInfo.hN[i] << " " << hapInfo.hP[i] << endl;
        }
        logfile << "Read " << hapInfo.hN.size() << " haplotype definitions" << endl;

        // read admixed haplotype:
        readAdmixSites( param.admix + ".sites", obs );
        logfile << "Read " << obs.size() << " admixed sites" << endl;
    }

    positions pos;
    if( param.fixedMapRate == -1 ) {
        geneticMap gMap;
        readGenMap( param.mapFile, gMap );
        logfile << "Read " << gMap.chr.size() << " lines from genetic map" << endl;

        //interpolate values from locs:
        interpGenMap( gMap, locs, pos );
        logfile << "Interpolated " << pos.cM.size() << " positions from genetic map" << endl;
        // test for nans in interpolated positions:
        int nancnt = 0;
        vector<int> nanix;
        for(int i=0; i < pos.cM.size(); i++ ) {
            if( std::isnan( pos.cM[i] ) ) {
                nancnt++;
                nanix.push_back( i );
            }
        }
        if( nancnt > 0 ) { // remove offending sites
            for(int i=nanix.size()-1; i>=0; i-- ) {
                sites.erase( sites.begin() + nanix[i] );
                pos.pos.erase( pos.pos.begin() + nanix[i] );
                pos.cM.erase( pos.cM.begin() + nanix[i] );
                obs.erase( obs.begin() + nanix[i] );
            }
            logfile << "Removed " << nancnt << " sites due to non-overlap with genetic map" << endl;
            logfile << "There are now " << sites.size() << " sites with " << pos.pos.size() << " corresponding genetic positions" << endl;
        }
    } else { // use fixed recombination rates
        pos.pos = locs;
        pos.cM.resize( locs.size() );
        double d;
        double cumcM = 0.0;
        for(int i=1; i < pos.pos.size(); i++ ) {
            d = ( pos.pos[i] - pos.pos[i-1] ) / 1000000; // dist in Mb
            cumcM += ( param.fixedMapRate * d ); // cumulative genetic distance, cM
            pos.cM[i] = cumcM;
        }
    }
    
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
    param.S = pos.pos.size() ;
    param.u2 = 1.0 - param.u1;

    hmmStates st;
    if( ( param.mode == 0 ) || ( param.mode == 2 ) || ( param.mode == 3 ) ) {
        generateXstates( hapInfo, st );
        // set theta to ~0 for reduced first pass:
        param.theta1 = numeric_limits<double>::epsilon();
        param.theta2 = numeric_limits<double>::epsilon();
    } else if( param.mode == 1 ) {
        generateStates( hapInfo, st );
        // set theta using watterson's estimator:
        param.theta1 = param.theta2 = 0.0;
        for(int m=1; m < param.n1; m++) { param.theta1 += 1.0/m; }
        for(int m=1; m < param.n2; m++) { param.theta2 += 1.0/m; }
        param.theta1 = 1.0 / param.theta1;
        param.theta2 = 1.0 / param.theta2;
    }
    logfile << "Generated " << st.states.size() << " states" << endl;

    // emissions probabilities:
    param.emit1_match = log( ( 2.0 * param.n1 + param.theta1 )/( 2.0 * ( param.n1 + param.theta1 ) ) );
    param.emit1_mismatch = log( param.theta1 / ( 2.0 * ( param.n1 + param.theta1 ) ) );
    param.emit2_match = log( ( 2.0 * param.n2 + param.theta2 )/( 2.0 * ( param.n2 + param.theta2 ) ) );
    param.emit2_mismatch = log( param.theta2 / ( 2.0 * ( param.n2 + param.theta2 ) ) );

    pathVec pvec( param );
    param.print_params(logfile, 1);

    ////////////////////////////////////
    // calculate alt. allele frequency for each ref population:
    vector<double> aleFrqP1( param.S );
    vector<double> aleFrqP2( param.S );
    for(int j=0; j < param.S; j++) {
        double acp1 = 0.0, acp2 = 0.0, ac0 = 0.0, ac1 = 0.0;
        for(int i=0; i < sites[j].size(); i++) {
            if( sites[j][i] == 1 ) { 
                if( hapInfo.hP[i] == 1 ) {
                    acp1 += 1;
                } else if( hapInfo.hP[i] == 2 ) {
                    acp2 += 1;
                }
                ac1 += 1;
            } else if( sites[j][i] == 0 ) {
                ac0 += 1;
            }
        }
        aleFrqP1[j] = acp1 / ( ac0 + ac1 );
        aleFrqP2[j] = acp2 / ( ac0 + ac1 );
    }
    
    ////////////////////////////////////
    // declare all hmm variables:
    vector<double> sprob( st.states.size(), 0.0 );
    double Pxa, Pxb;
    vector<vector<double> > fwd(param.S, vector<double>(st.states.size(), 0.0));
    vector<double> fwdS;
    vector<vector<double> > bwd(param.S, vector<double>(st.states.size(), 0.0));
    vector<vector<double> > pprob(param.S, vector<double>(st.states.size(), 0.0));
    vector<vector<double> > vit;

    //////////////////////////
    if( param.mode != 3 ) { 
    //////////////////////////

    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        getsprobX( sites[0], param, st, obs, sprob );
        param.passAcc = 1; // force LSE sort on first pass
    } else if( param.mode == 1 ) {
        getsprob( sites[0], param, st, obs, sprob );
        param.passAcc = param.highAccuracy; // set to user-selected state
        std::fill( pvec.pswitch.begin(), pvec.pswitch.end(), 1 );
    }

    logfile << "Starting forward algorithm..." << endl;;
    forward( sites, pos, param, st, obs, sprob, pvec.pswitch, fwd );
    fwdS = fwd[ fwd.size()-1 ];
    logSumExp( fwdS, Pxa, param.passAcc );

    logfile << "Starting backward algorithm..." << endl;
    backward( sites, pos, param, st, obs, sprob, pvec.pswitch, bwd, Pxb );
    //
    logfile << setprecision(20) << "Pxa= " << Pxa << endl;
    logfile << setprecision(20) << "Pxb= " << Pxb << endl;
    logfile << "diff= " << Pxa - Pxb << endl;
    // printMat( bwd );
    // cout << "P(x) diff=" << Pxa - Pxb << endl;

    logfile << "Starting posterior decoding..." << endl;
    //postDecode( fwd, bwd, st, Pxa, pprob, pvec, logfile);
    postDecode( fwd, bwd, st, Pxa, pprob, pvec.pppath, pvec.ppprob, logfile);

    if( ( param.mode == 0 ) || (param.mode == 2 ) || ( param.viterbi == 1 ) ) {
        logfile << "Starting Viterbi algorithm..." << endl;
        vector<vector<double> > vit(param.S, vector<double>(st.states.size(), 0.0));
        viterbi( sites, pos, param, st, obs, sprob, vit, pvec );
        //printMat( vit );
    } else {
        logfile << "Skipping Viterbi algorithm..." << endl;
    }

    if( ( param.mode == 0 ) || ( param.mode == 2 ) ) {
        // find stretches/switches in viterbi path:
        string prevHap = pvec.vpath[0];
        int prevCnt = 1;
        vector<string> elem = split( prevHap, '-' );
        string prevPop = elem[0];
        for(int j=1; j < param.S; j++ ) {
            if( pvec.vpath[j] != pvec.vpath[j-1] ) {
                pvec.rleHap.push_back( prevHap );
                pvec.rleCnt.push_back( prevCnt );
                pvec.rlePrevPop.push_back( prevPop );
                //cout << "rleHap=" << prevHap << "\trleCnt=" << prevCnt << "\trlePop=" << elem[0] << endl;
                prevHap = pvec.vpath[j];
                elem = split( prevHap, '-' );
                prevPop = elem[0];
                prevCnt = 1;
            } else {
                prevCnt += 1;
            }
        }
        pvec.rleHap.push_back( prevHap );
        pvec.rleCnt.push_back( prevCnt );
        pvec.rlePrevPop.push_back( prevPop );
        //cout << "rleHap=" << prevHap << "\trleCnt=" << prevCnt << "\trlePop=" << elem[0] << endl;
        // find cross popultation "islands":
        vector<int> island;
        if( pvec.rlePrevPop[0] != pvec.rlePrevPop[1] ) {
            island.push_back( 1 );
        } else {
            island.push_back( 0 );
        }
        for(int i=1; i < pvec.rleHap.size() - 1; i++) {
            if( ( pvec.rlePrevPop[i-1] != pvec.rlePrevPop[i] ) || ( pvec.rlePrevPop[i+1] != pvec.rlePrevPop[i] ) ) {
                island.push_back( 1 );
            } else {
                island.push_back( 0 );
            }
        }
        int l = pvec.rleHap.size() - 1;
        if( l > 0 ) { // skip if there is only one viterbi stretch
            if( pvec.rlePrevPop[ l ] != pvec.rlePrevPop[ l-1 ] ) {
                island.push_back( 1 );
            } else {
                island.push_back( 0 );
            }
            // find stretches less than X sites, mark for follow-up
            int stepix = 0;
            for(int i=0; i < pvec.rleHap.size(); i++ ) {
                //cout << "rleHap=" << pvec.rleHap[i] << "\trleCnt=" << pvec.rleCnt[i] << "\trlePop=" << pvec.rlePrevPop[i] << "\tisland=" << island[i] << endl;
                if( ( pvec.rleCnt[i] <= param.gcsens ) && ( island[i] == 1 ) ) {
                    for(int j=0; j < pvec.rleCnt[i]; j++ ) {
                        //cout << "pswitch ix: " << stepix+j << endl;
                        pvec.pswitch[ stepix+j ] = 1;
                    }
                }
                stepix += pvec.rleCnt[i];
            }
            vector<int> swIx;
            for(int j=0; j < param.S; j++ ) {
                if( pvec.pswitch[j] == 1 )
                    swIx.push_back( j );
            }
            //int width = 5; // how much to expand on each side of switch marks
            // expand marks around switches:
            vector<int> repl(2);
            for(int j=0; j < swIx.size(); j++) {
                repl[0] = swIx[j] - param.width;
                repl[1] = swIx[j] + param.width;
                if( repl[0] < 0 )
                    repl[0] = 0;
                if( repl[1] > (param.S-1) )
                    repl[1] = (param.S-1);
                for(int k=repl[0]; k <= repl[1]; k++ )
                    pvec.pswitch[k] = 1;
            }
            // for(int i=0; i<pvec.pswitch.size(); i++) { if( pvec.pswitch[i] == 1 ) cout << i << " "; } cout << endl;
        }
    }

    if( ( param.mode == 0 ) || ( param.mode == 1 ) ) {
        logfile << "Starting path output" << endl;
        pathOutput( pvec, st, pos, pprob, param, aleFrqP1, aleFrqP2, pathfile );
    }

    //////////////////////////
    } else { // end skip block for mode 3
    //////////////////////////
        // read previous path file
        // string pathf0 = param.outfname + ".path0";
        readPath( param.pathfile, pvec );
        logfile << "Read first-pass results from " << param.pathfile << endl;
    }

    ////////////////
    if( param.fixPswitch >= 0 ) { // testing only
        std::fill( pvec.pswitch.begin(), pvec.pswitch.end(), 0 );
        for(int j=0; j < param.fixPswitch; j++) {
            pvec.pswitch[j] = 1;
        }
    }
    ////////////////

    // restart hmm:
    if( ( param.mode == 2 ) || ( param.mode == 3 ) ) {
        logfile.precision (ss);
        logfile << endl << "Restart hmm" << endl;

        // set theta using watterson's estimator:
        param.theta1 = param.theta2 = 0.0;
        for(int m=1; m < param.n1; m++) { param.theta1 += 1/m; }
        for(int m=1; m < param.n2; m++) { param.theta2 += 1/m; }
        param.theta1 = 1.0 / param.theta1;
        param.theta2 = 1.0 / param.theta2;
        // emissions probabilities:
        param.emit1_match = log( ( 2.0 * param.n1 + param.theta1 )/( 2.0 * ( param.n1 + param.theta1 ) ) );
        param.emit1_mismatch = log( param.theta1 / ( 2.0 * ( param.n1 + param.theta1 ) ) );
        param.emit2_match = log( ( 2.0 * param.n2 + param.theta2 )/( 2.0 * ( param.n2 + param.theta2 ) ) );
        param.emit2_mismatch = log( param.theta2 / ( 2.0 * ( param.n2 + param.theta2 ) ) );

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
        forward( sites, pos, param, st2, obs, sprob, pvec.pswitch, fwd );
        fwdS.resize( fwd[fwd.size()-1].size() , 0.0 );
        fwdS = fwd[ fwd.size()-1 ];
        logSumExp( fwdS, Pxa, param.passAcc );

        logfile << "Starting backward algorithm..." << endl;
        backward( sites, pos, param, st2, obs, sprob, pvec.pswitch, bwd, Pxb );

        logfile << setprecision(20) << "Pxa= " << Pxa << endl;
        logfile << setprecision(20) << "Pxb= " << Pxb << endl;
        logfile << "diff= " << Pxa - Pxb << endl;
        // cout << "P(x) diff=" << Pxa - Pxb << endl;

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
        std::fill( pvec.transPGC.begin(), pvec.transPGC.end(), 0.0 );
        pathOutput( pvec, st2, pos, pprob, param, aleFrqP1, aleFrqP2, pathfile );

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



