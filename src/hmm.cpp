/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

void generateStates(const hapDef &hapInfo, hmmStates &st ) {
    int nrow = hapInfo.hN.size() ; // - 1 ;
    // G null states:
    vector<string> hapNamesG, hapNamesX, hapPopG, hapPopX;
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        st.Xhap.push_back( hapInfo.hN[j] );
        st.Xpop.push_back( hapInfo.hP[j] );
        st.Xindx.push_back( j ); // check if this is right
        st.Ghap.push_back( 0 );
        st.Gpop.push_back( 0 );
        st.Gindx.push_back( j ); // check if this is right
        st.Epop.push_back( hapInfo.hP[j] ); // X pop
        // keep original naming scheme:
        hapNamesG.push_back( "0" );
        hapPopG.push_back( "0" );
        hapNamesX.push_back( hapInfo.hapName[j] );
        hapPopX.push_back( hapInfo.hapPop[j] );
    }
    // All combinations of X and G:
    for(int i=0; i < nrow; i++) {
        for (int j=0; j < nrow; j++) {
            // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
            st.Xhap.push_back( hapInfo.hN[j] );
            st.Xpop.push_back( hapInfo.hP[j] );
            st.Xindx.push_back( j );
            st.Ghap.push_back( hapInfo.hN[i] );
            st.Gpop.push_back( hapInfo.hP[i] );
            st.Gindx.push_back( i );
            st.Epop.push_back( hapInfo.hP[i] ); // G pop
            // keep original naming scheme:
            hapNamesG.push_back( hapInfo.hapName[i] );
            hapPopG.push_back( hapInfo.hapPop[i] );
            hapNamesX.push_back( hapInfo.hapName[j] );
            hapPopX.push_back( hapInfo.hapPop[j] );
        }
    }
    // state concatenation:
    nrow = st.Xhap.size();
    //cout << "i\tXhap\tXpop\tXindx\tGhap\tGpop\tGindx\tstates\tEmitPop" << endl;
    for(int i=0; i < nrow; i++) {
        string tmp;
        // tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]) + "-" + std::to_string(st.Gpop[i]) + "-" + std::to_string(st.Ghap[i]);
        tmp = hapPopX[i] + "-" + hapNamesX[i] + "-" + hapPopG[i] + "-" + hapNamesG[i];
        //tmp = st.Xpop[i] + "-" + st.Xhap[i] + "-" + st.Gpop[i] + "-" + st.Ghap[i];
        st.states.push_back( tmp );
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" <<
        //     st.Ghap[i] << "\t" << st.Gpop[i] << "\t" << st.Gindx[i] << "\t" << st.states[i] << "\t" << st.Epop[i] << endl;
    }
}

void generateXstates(const hapDef &hapInfo, hmmStates &st ) {
    int nrow = hapInfo.hN.size() ; // - 1 ;
    // Haplotype states only:
    vector<string> hapNames, hapPop;
    for (int j=0; j < nrow; j++) {
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
        st.Xhap.push_back( hapInfo.hN[j] );
        st.Xpop.push_back( hapInfo.hP[j] );
        st.Xindx.push_back( j );
        st.Ghap.push_back( 0 ); // fill this vector with null GC states
        st.Epop.push_back( hapInfo.hP[j] ); // X pop
        // keep original naming scheme:
        hapNames.push_back( hapInfo.hapName[j] );
        hapPop.push_back( hapInfo.hapPop[j] );

    }
    // state concatenation:
    nrow = st.Xhap.size();
    //cout << "i\tXhap\tXpop\tXindx\tGhap\tstates\tEmitPop" << endl;
    string tmp;
    for(int i=0; i < nrow; i++) {
        // tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]);
        tmp = hapPop[i] + "-" + hapNames[i];
        st.states.push_back( tmp );
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" << st.Ghap[i] << "\t" << st.states[i] << "\t" << st.Epop[i] << endl;
    }
}

void getsprob( 
        const vector<int> &sites0, 
        const parameters &p, 
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob ) {
    double e, sprobG_i;
    double g = 1.1 / 100000.0 * p.f; // 1.1cM/Mb * f = 1Morgan/kb
    vector<double> sprobX; // [ pop1, pop2 ]
    sprobX.push_back( log( p.u1 / p.n1 ) ); // pop1
    sprobX.push_back( log( ( 1.0 - p.u1 ) / p.n2 ) ); // pop2
    double sprob_G0 = log( ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * ( p.n1 + p.n2 ) + g * p.T ) );
    vector<double> sprob_G1; // [ pop1, pop2 ]
    sprob_G1.push_back( log( ( g * p.T ) / ( p.lam * ( p.n1 + p.n2 ) + g * p.T ) * ( p.u1 / p.n1 ) ) ); // pop1
    sprob_G1.push_back( log( ( g * p.T ) / ( p.lam * ( p.n1 + p.n2 ) + g * p.T ) * ( ( 1.0 - p.u1 ) / p.n2 ) ) ); // pop2
    //double tmpsum = 0.0;
    for(int i=0; i < st.states.size(); i++) {
        if( st.Ghap[i] == 0 ) { // match on Xpop if G==0
            sprobG_i = sprob_G0; // use G0
        } else { // match on Gpop if G!=0
            sprobG_i = sprob_G1[ st.Gpop[i] - 1 ]; // use pop-specific G1
        }
        if( sites0[ st.Gindx[i] ] == obs[0] ) {
            e = p.emitMatch[ st.Epop[i] - 1 ];
        } else {
            e = p.emitMismatch[ st.Epop[i] - 1 ];
        }
        sprob[i] = sprobX[ st.Xpop[i] - 1 ] + sprobG_i + e;
        //tmpsum += exp( sprobX[ st.Xpop[i] - 1 ] + sprobG_i );
        //cout << "Xpop=" << st.Xpop[i] << " Gpop=" << st.Gpop[i] << " whichGpop=" << whichGpop << " Xprob=" << exp(sprobX[st.Xpop[i]-1]) << " Gprob= " << exp(sprobG_i) << endl;
    }
    //cout << "sprob sum = " << tmpsum << endl;
}

void getsprobX(
        const vector<int> &sites0,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob ) {
    double e;
    vector<double> sprobX; // [ pop1, pop2 ]
    sprobX.push_back( log( p.u1 / p.n1 ) ); // pop1
    sprobX.push_back( log( ( 1.0 - p.u1 ) / p.n2 ) ); // pop2
    //double tmpsum = 0.0;
    for(int i=0; i < st.states.size(); i++) {
        if( sites0[ st.Xindx[i] ] == obs[0] ) {
            e = p.emitMatch[ st.Xpop[i] - 1 ];
        } else {
            e = p.emitMismatch[ st.Xpop[i] - 1 ];
        }
        sprob[i] = sprobX[ st.Xpop[i] - 1 ] + e;
        //tmpsum += exp( sprobX[ st.Xpop[i] - 1 ] );
    }
    //cout << "sprob sum = " << tmpsum << endl;
}

double lookupXtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trXbin ) {
    // first pass at this site: calculate all reusable values:
    if( trXbin[6] == 99 ) {
        // Here, the input r is genetic distance in Morgans, d in kb
        // convert to rate:
        double rrate = r / d;
        double trX1, trX2, trX3, trX4, trX5, trX6;
        // pop1:
        double rho1 = 4.0 * p.Ne1 * rrate;
        trX1 = ( 1.0 - exp( -rrate * d * p.T ) ) * p.u1 / p.n1;
        trX2 = exp( -rrate * d * p.T ) * ( 1.0 - exp( -rho1 * d / p.n1 ) ) / p.n1 + 
               ( 1.0 - exp( -rrate * d * p.T )) * p.u1 / p.n1 ;
        trX3 = exp( -rrate * d * p.T ) * exp( -rho1 * d / p.n1 ) + 
               exp( -rrate * d * p.T ) * ( 1.0 - exp( -rho1 * d / p.n1 ) ) / p.n1 + 
               ( 1.0 - exp( -rrate * d * p.T )) * p.u1 / p.n1 ;
        // pop2:
        double rho2 = 4.0 * p.Ne2 * rrate;
        trX4 = ( 1.0 - exp( -rrate * d * p.T ) ) * ( 1.0 - p.u1 ) / p.n2;
        trX5 = exp( -rrate * d * p.T ) * ( 1.0 - exp( -rho2 * d / p.n2 ) ) / p.n2 + 
               ( 1.0 - exp( -rrate * d * p.T )) * ( 1.0 - p.u1 ) / p.n2;
        trX6 = exp( -rrate * d * p.T ) * exp( -rho2 * d / p.n2 ) + 
               exp( -rrate * d * p.T ) * ( 1.0 - exp( -rho2 * d / p.n2 ) ) / p.n2 + 
               ( 1.0 - exp( -rrate * d * p.T )) * ( 1.0 - p.u1 ) / p.n2 ;
        trXbin[0] = log( trX1 );
        trXbin[1] = log( trX2 );
        trXbin[2] = log( trX3 );
        trXbin[3] = log( trX4 );
        trXbin[4] = log( trX5 );
        trXbin[5] = log( trX6 );
        trXbin[6] = 1.0;
        // prob sums check:
        // cout << "trX sum from p1=" << exp(trXbin[2]) + exp(trXbin[1])*(p.n1-1) + exp(trXbin[3])*p.n2;
        // cout << "\ttrX sum from p2=" << exp(trXbin[5]) + exp(trXbin[4])*(p.n2-1) + exp(trXbin[0])*p.n1 << endl;
    }
    int type;
    double trX;
    if( st.Xpop[from] != st.Xpop[to] ) {
        type = 0;
    } else {
        if( st.Xhap[from] != st.Xhap[to] ) {
            type = 1;
        } else { // if( st.Xhap[from] == st.Xhap[to] )
            type = 2;
        }
    }
    if( st.Xpop[to] == 2 ) // shift type index by 3 if transitioning to a haplotype in pop #2.
        type += 3;
    //cout << " type" << type << " ";
    return( trXbin[ type ] );
}

inline double lookupGtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trGbin ) {
    // first pass at this site: calculate all reusable values:
    if( trGbin[10] == 99 ) {
        // Here, the input r is genetic distance in Morgans, d in kb
        // convert to rate:
        double rrate = r / d;
        double int1, int2p1, int2p2, tmp;
        double g = p.f * rrate;
        double expdlgt = exp( -d * ( p.lam + g * p.T ) );
        double expdl = exp( -d * p.lam );
        int1 = ( ( 1.0 - expdlgt ) * p.lam ) / ( p.lam + g * p.T );
        int2p1 = ( p.u1 - expdl * p.u1 + 
               ( ( -1.0 + expdlgt ) * p.lam * p.u1 ) / ( p.lam + g * p.T )
               ) / p.n1;
        int2p2 = ( ( 1.0 - p.u1 ) - expdl * ( 1.0 - p.u1 ) + 
               ( ( -1.0 + expdlgt ) * p.lam * ( 1.0 - p.u1 ) ) / ( p.lam + g * p.T )
               ) / p.n2;
        // 0 //
        tmp = log( expdl * expdlgt + int1 );
        trGbin[0] = tmp;
        trGbin[5] = tmp;
        // 1 //
        trGbin[1] = log( expdl * ( 1.0 - expdlgt ) * p.u1 / p.n1 + int2p1 );
        // 2 //
        tmp = log( int1 );
        trGbin[2] = tmp;
        trGbin[7] = tmp;
        // 3 //
        trGbin[3] = log( expdl + int2p1 );
        // 4 //
        trGbin[4] = log( int2p1 );
        // 5 //
        // 6 //
        trGbin[6] = log( expdl * ( 1.0 - expdlgt ) * ( 1.0 - p.u1 ) / p.n2 + int2p2 );
        // 7 //
        // 8 //
        trGbin[8] = log( expdl + int2p2 );
        // 9 //
        trGbin[9] = log( int2p2 );
        /////// flag:
        trGbin[10] = 1.0;
        // prob sums check:
        // for(int z=0; z<trGbin.size(); z++) { cout << exp(trGbin[z]) << endl; }
        // cout << "trGbinsum0=" << exp(trGbin[0]) + exp(trGbin[1])*p.n1 + exp(trGbin[6])*p.n2;
        // cout << "\ttrGbinsum1=" << exp(trGbin[2]) + exp(trGbin[3]) + exp(trGbin[4])*(p.n1-1) + exp(trGbin[9])*(p.n2);
        // cout << "\ttrGbinsum2=" << exp(trGbin[2]) + exp(trGbin[8]) + exp(trGbin[9])*(p.n2-1) + exp(trGbin[4])*(p.n1) << endl;
    }
    // return a pre-calculated value:
    int type = 4;
	if ( st.Ghap[from] == 0 ) {
		if ( st.Ghap[from] == st.Ghap[to] )
			type = 0;
		else if ( st.Ghap[to] != 0 )
			type = 1;
	} else {
		if ( st.Ghap[to] == 0 )
			type = 2;
		else if ( st.Ghap[from] == st.Ghap[to] )
			type = 3;
	}
    if( st.Gpop[to] == 2 ) // shift type index by 5 if transitioning to a haplotype in pop #2.
        type += 5;
    //cout << " from" << st.Ghap[from] << " to" << st.Ghap[to];
    //cout << " type" << type << endl;
    // lookup and return a value already seen:
    return( trGbin[ type ] );
}

void forward(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        const vector<int> &pswitch,
        vector<vector<double> > &fwd ) {
    // starting prob:
    fwd[0] = sprob;
    double lsum, trX, trG, e, d, r;
    vector<double> trXbin(7,99); // 6 slots for transitions, 1 for flag
    vector<double> trGbin(11,99); // 10 slots for transitions, 1 for flag
    vector<int> siteIndx;
    vector<double> tmp;
    //cout << "pswitch"<<pswitch.size() << endl;
    for(int j=1; j < sites.size() ; j++ ) {
        //cout << "j = " << j << "\t" << "pswitch=" << pswitch[j] << " Nstates=" << fwd[j].size() << endl;
        d = pos.pos[j] - pos.pos[j-1]; // physical distance in kb
        r = pos.cM[j] - pos.cM[j-1]; // genetic distance in Morgans
        if( pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Xindx;
        } else {
            siteIndx = st.Gindx;
        }
        //cout << "fwd[j].size() = " << fwd[j].size() << endl;
        for(int t=0; t < fwd[j].size(); t++) {
            tmp.resize( fwd[j-1].size(), 0.0 );
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int f=0; f < fwd[j-1].size(); f++) {
                trX = lookupXtrans( t, f, d, r, st, p, trXbin );
                if( ( pswitch[j] == 1 ) || ( pswitch[j-1] == 1 ) ) {
                    //cout << "\ttrx + trG" << endl;
                    trG = lookupGtrans( t, f, d, r, st, p, trGbin );
                    tmp[f] = fwd[j-1][f] + trX + trG;
                } else {
                    tmp[f] = fwd[j-1][f] + trX;
                }
            } // end 'from' loop
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                e = p.emitMatch[ st.Epop[t] - 1 ];
            } else {
                e = p.emitMismatch[ st.Epop[t] - 1 ];
            }
            // cout << "j=" << j << endl;
            logSumExp( tmp, lsum, p.passAcc );
            fwd[j][t] = e + lsum;
        } // end 'to' loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // end j site loop
}

void backward(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        const vector<int> &pswitch,
        vector<vector<double> > &bwd,
        double &Pxb ) {
    double e, lsum, trX, trG, d, r;
    double negInf = - std::numeric_limits<double>::infinity();
    double log_eps = log(numeric_limits<double>::epsilon());
    vector<double> trXbin(7,99); // 6 slots for transitions, 1 for flag
    vector<double> trGbin(11,99); // 10 slots for transitions, 1 for flag
    vector<int> siteIndx;
    vector<double> tmp;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        //cout << "j = " << j << "\t" << "pswitch=" << pswitch[j] << " Nstates=" << bwd[j].size() << endl;
        //d = dvec[j+1] - dvec[j];
        d = pos.pos[j+1] - pos.pos[j];
        r = pos.cM[j+1] - pos.cM[j];
        if( pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Xindx;
        } else {
            siteIndx = st.Gindx;
        }
        for(int f=0; f < bwd[j].size(); f++ ) {
            tmp.resize( bwd[j+1].size(), 0.0 );
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int t=0; t < bwd[j+1].size(); t++ ) {
                trX = lookupXtrans( t, f, d, r, st, p, trXbin);
                if( sites[j+1][ siteIndx[t] ] == obs[j+1] ) {
                    e = p.emitMatch[ st.Epop[t] - 1 ];
                } else {
                    e = p.emitMismatch[ st.Epop[t] - 1 ];
                }
                if( ( pswitch[j] == 1 ) || ( pswitch[j+1] == 1 ) ) {
                    //cout << "\ttrx + trG" << endl;
                    trG = lookupGtrans( t, f, d, r, st, p, trGbin);
                    tmp[t] = bwd[j+1][t] + trX + trG + e ;
                } else {
                    //cout << "\ttrx" << endl;
                    tmp[t] = bwd[j+1][t] + trX + e ;
                }
            } // t loop
            logSumExp( tmp, lsum, p.passAcc );
            bwd[j][f] = lsum;
        } // f loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // j loop
    // termination:
    tmp.resize( bwd[0].size(), 0.0 );
    std::fill( tmp.begin(), tmp.end(), 0.0 );
    for(int f=0; f < bwd[0].size(); f++ ) {
        tmp[f] = bwd[0][f] + sprob[f];
    }
    logSumExp( tmp, Pxb, p.passAcc );
}

void printMat( const vector<vector<double> > &mat ) {
    for(int i=0; i < mat[0].size(); i++ ) {
        for(int j=0; j < mat.size(); j++ ) {
            cout << setprecision(15) << mat[j][i] << " ";
        }
        cout << endl;
    }
}

void writeMat( const vector<vector<double> > &mat, ofstream &matfile ) {
    for(int i=0; i < mat[0].size(); i++ ) {
        for(int j=0; j < mat.size(); j++ ) {
            matfile << setprecision(15) << mat[j][i] << " ";
        }
        matfile << endl;
    }
}

void writeTmat( const vector<vector<double> > &mat, ofstream &matfile ) {
    for(int j=0; j < mat.size(); j++ ) {
        for(int i=0; i < mat[j].size(); i++ ) {
            matfile << setprecision(15) << mat[j][i] << " ";
        }
        matfile << endl;
    }
}

void logSumExp( vector<double> &vec, double &lse, const int &hp ) {
    if( hp == 0 ) {
        double max = *max_element(vec.begin(), vec.end());
        double sum = 0.0;
        double a;
        for(int i=0; i < vec.size(); i++ ) {
            a = vec[i] - max;
            if( a > -37 )
                sum += exp( a );
        }
        lse = max + log(sum) ;
    } else if( hp == 1 ) {
        double sum = 0.0;
        sort(vec.begin(), vec.end());
        double max = vec.back();
        double a;
        for(int i=0; i < vec.size(); i++ )
        {
            sum += exp( vec[i] - max );
        }
        lse = max + log(sum) ;
    } else if( hp == 2 ) {
        double sum = 0.0;
        sort(vec.begin(), vec.end());
        double max = vec.back();
        double a;
        for(int i=0; i < vec.size(); i++ )
        {
            a = vec[i] - max;
            if (a > -37) {
                //cout << "break" << endl;
                break;
            } //else { cout << "nobreak" << endl; }
            sum += exp( a );
        }
        lse = max + log(sum) ;
    }
}

void postDecode(
        const vector<vector<double> > &fwd,
        const vector<vector<double> > &bwd,
        const hmmStates &st,
        const double &Px,
        vector<vector<double> > &pprob,
        vector<string> &pppath,
        vector<double> &ppprob,
        //pathVec &pvec,
        ofstream &logfile
        ) {
    double negInf = - std::numeric_limits<double>::infinity();
    // decoding:
    for(int j=0; j < pprob.size(); j++ ) {
        for(int i=0; i < pprob[j].size(); i++ ) {
            pprob[j][i] = fwd[j][i] + bwd[j][i] - Px;
        }
    }
    // normalize:
    double cmax, csum;
    for(int j=0; j < pprob.size(); j++ ) {
        cmax = negInf;
        // determine max value:
        for(int i=0; i < pprob[j].size(); i++ ) {
            if( pprob[j][i] > cmax ) 
                cmax = pprob[j][i];
        }
        // subtract max and determine sum:
        csum = 0.0;
        for(int i=0; i < pprob[j].size(); i++ ) {
            pprob[j][i] = pprob[j][i] - cmax;
            csum += exp( pprob[j][i] );
        }
        // divide by sum:
        for(int i=0; i < pprob[j].size(); i++ ) {
            pprob[j][i] = exp( pprob[j][i] ) / csum;
        }
    }
    // find most probable path:
    int maxix = -1;
    double pmax;
    for(int j=0; j < pprob.size(); j++ ) {
        pmax = negInf;
        for(int i=0; i < pprob[j].size(); i++ ) {
            if( pprob[j][i] > pmax ) {
                pmax = pprob[j][i];
                maxix = i;
            }
        }
        ppprob[j] = pprob[j][maxix];
        pppath[j] = st.states[maxix];
    }
}

void viterbi(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &vit,
        pathVec &pvec ) {
    // starting prob:
    vit[0] = sprob;
    // recursion:
    double trX, trG, vmax, tmp, e, d, r;
    double negInf = - std::numeric_limits<double>::infinity();
    vector<double> trXbin(7,99); // 6 slots for transitions, 1 for flag
    vector<double> trGbin(11,99); // 10 slots for transitions, 1 for flag
    vector<int> siteIndx;
    if( p.mode == 1 ) {
        siteIndx = st.Gindx;
    } else {
        siteIndx = st.Xindx;
    }
    for(int j=1; j < sites.size() ; j++ ) {
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        for(int t=0; t < st.states.size(); t++ ) {
            vmax = negInf;
            for(int f=0; f < st.states.size(); f++ ) {
                trX = lookupXtrans( t, f, d, r, st, p, trXbin);
                if( p.mode == 1 ) {
                    trG = lookupGtrans( t, f, d, r, st, p, trGbin);
                    tmp = vit[j-1][f] + trX + trG;
                } else {
                    tmp = vit[j-1][f] + trX;
                }
                if( tmp > vmax )
                    vmax = tmp;
            } // end from loop
            // emission prob:
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                e = p.emitMatch[ st.Epop[t] - 1 ];
            } else {
                e = p.emitMismatch[ st.Epop[t] - 1 ];
            }
            vit[j][t] = e + vmax;
        } // end to loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // end j loop
    // termination:
    int maxix = -1;
    // find last max state:
    vmax = negInf;
    maxix = -1;
    int j = sites.size()-1;
    for(int i=0; i < vit[0].size(); i++ ) {
        if( vit[j][i] > vmax ) {
            vmax = vit[j][i];
            maxix = i;
        }
    }
    pvec.vpath[j] = st.states[maxix];
    pvec.vprob[j] = exp( vit[j][maxix] );
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < st.states.size(); f++) {
            trX = lookupXtrans( t, f, d, r, st, p, trXbin);
            if( p.mode == 1 ) {
                trG = lookupGtrans( t, f, d, r, st, p, trGbin);
                tmp = vit[j][f] + trX + trG;
            } else {
                tmp = vit[j][f] + trX;
            }
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        pvec.vpath[j] = st.states[ maxix ];
        pvec.vprob[j] = exp( vit[j][maxix] );
        // cout << vpath[j] << " " << vmax << " " << maxix << endl;
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    }
}

void viterbi2(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &vit,
        pathVec &pvec ) {
    // starting prob, same as fwd:
    vit[0] = sprob;
    // recursion:
    double trX, trG, vmax, tmp, e, d, r;
    double negInf = - std::numeric_limits<double>::infinity();
    //double log_eps = log(numeric_limits<double>::epsilon());
    vector<double> trXbin(7,99); // 6 slots for transitions, 1 for flag
    vector<double> trGbin(11,99); // 10 slots for transitions, 1 for flag
    vector<int> siteIndx;
    unsigned int whichPop;
    for(int j=1; j < sites.size() ; j++ ) {
        //d = dvec[j] - dvec[j-1];
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        if( pvec.pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Xindx;
        } else {
            siteIndx = st.Gindx;
        }
        for(int t=0; t < vit[j].size(); t++ ) {
            vmax = negInf;
            for(int f=0; f < vit[j-1].size(); f++ ) {
                trX = lookupXtrans( t, f, d, r, st, p, trXbin);
                if( ( pvec.pswitch[j] == 1 ) | ( pvec.pswitch[j-1] == 1 ) ) {
                    trG = lookupGtrans( t, f, d, r, st, p, trGbin);
                    tmp = vit[j-1][f] + trX + trG;
                } else {
                    tmp = vit[j-1][f] + trX;
                }
                if( tmp > vmax )
                    vmax = tmp;
            } // end from loop
            // determine population of 'to' haplotype:
            if( pvec.pswitch[j] == 1 ) { // if GC chain...
                if( st.Gpop[t] == 0 ) { // ... and G==0
                    whichPop = st.Xpop[t]; // use Xpop when G==0
                } else { // G!=0
                    whichPop = st.Gpop[t]; // use Gpop when G!=0
                }
            } else { // no GC chain
                whichPop = st.Xpop[t]; // use Xpop when G chain does not exist
            }
            // emission prob:
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                e = p.emitMatch[ st.Epop[t] - 1 ];
            } else {
                e = p.emitMismatch[ st.Epop[t] - 1 ];
            }
            vit[j][t] = e + vmax;
        } // end to loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // end j loop
    // termination:
    int maxix = -1;
    // find last max state:
    vmax = negInf;
    maxix = -1;
    int j = sites.size()-1;
    for(int i=0; i < vit[0].size(); i++ ) {
        if( vit[j][i] > vmax ) {
            vmax = vit[j][i];
            maxix = i;
        }
    } 
    pvec.vpath[j] = st.states[maxix];
    pvec.vprob[j] = exp( vit[j][maxix] );
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < vit[j].size(); f++) {
            trX = lookupXtrans( t, f, d, r, st, p, trXbin);
            //if( pswitch[j] == 1 ) {
            if( ( pvec.pswitch[j] == 1 ) || ( pvec.pswitch[j+1] == 1 ) ) {
                trG = lookupGtrans( t, f, d, r, st, p, trGbin);
                tmp = vit[j][f] + trX + trG;
            } else {
                tmp = vit[j][f] + trX;
            }
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        pvec.vpath[j] = st.states[ maxix ];
        pvec.vprob[j] = exp( vit[j][maxix] );
        // cout << vpath[j] << " " << vmax << " " << maxix << endl;
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    }
}

void which_max( const vector<double> &vec, int maxindx ) {
    double tmp = - std::numeric_limits<double>::infinity();
    double maxelement;
    for(int i=0; i < vec.size(); i++ ) {
        tmp = vec[i];
        if( tmp > maxelement )
            maxelement = tmp;
    }
}

void max( const vector<double> &vec, double maxelement ) {
    double tmp = - std::numeric_limits<double>::infinity();
    for(int i=0; i < vec.size(); i++ ) {
        tmp = vec[i];
        if( tmp > maxelement )
            maxelement = tmp;
    }
}

pathVec::pathVec( parameters param ) {
    pppath.resize( param.S );
    ppprob.resize( param.S, 0.0 );
    pswitch.resize( param.S, 0 );
//
    vpath.resize( param.S, "0" );
    vprob.resize( param.S, 0.0 );
//
    gcprob.resize( param.S, 0.0 );
    transPGC.resize( param.S, 0.0 );
//
    if( ( param.mode == 2 ) || ( param.mode == 3 ) ) {
        pppath2.resize( param.S );
        ppprob2.resize( param.S, 0.0 );
        vpath2.resize( param.S );
        vprob2.resize( param.S, 0.0 );
    }
}

void pathOutput( pathVec &pvec, hmmStates &st, positions &pos, vector<vector<double> > &pprob, parameters &param, vector<double> &aleFrqP1, vector<double> &aleFrqP2, ofstream &pathfile ) {
    if( ( param.mode == 1 ) || ( param.mode == 2 ) || ( param.mode == 3 ) ) {
        // calculate total probability of a cross-population gene conversion:
        for(int j=0; j < pprob.size(); j++ ) {
            for(int i=0; i < pprob[j].size(); i++ ) {
                if( st.Ghap[i] != 0 ) { // gene conversion state
                    pvec.gcprob[j] += pprob[j][i];
                }
                if( st.Xpop[i] == 1 ) {
                    for(int k=0; k < pprob[j].size(); k++ ) {
                        if( st.Gpop[k] == 2 ) {
                            pvec.transPGC[j] += pprob[j][i] * pprob[j][k];
                        }
                    }
                } else if( st.Xpop[i] == 2 ) {
                    for(int k=0; k < pprob[j].size(); k++ ) {
                        if( st.Gpop[k] == 1 ) {
                            pvec.transPGC[j] += pprob[j][i] * pprob[j][k];
                        }
                    }
                }
            }
        }
    } // endif
    // output path and probabilites:
    pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\tAAFp1\tAAFp2\t";
    if( ( param.mode == 2 ) || ( param.mode == 3 ) ) {
        if( param.viterbi == 1 ) 
            pathfile << "Vpath2\tVpathProb2\t";
        pathfile << "Ppath2\tPpathProb2\t";
    }
    pathfile << "GCprob\tGCprobXPop" << endl;
    int intpos;
    for(int j=0; j < pprob.size(); j++ ) {
        intpos = static_cast<int>( pos.pos[j]*1000+0.5 );
        pathfile << j << "\t" << intpos << "\t" << pvec.vpath[j] << "\t" << pvec.vprob[j] << "\t" << pvec.pppath[j] << "\t" << pvec.ppprob[j] << "\t" << pvec.pswitch[j] << "\t" ;
        pathfile << aleFrqP1[j] << "\t" << aleFrqP2[j] << "\t";
        if( ( param.mode == 2 ) || ( param.mode == 3 ) ) {
            if( param.viterbi == 1 ) 
                pathfile << pvec.vpath2[j] << "\t" << pvec.vprob2[j] << "\t";
            pathfile << pvec.pppath2[j] << "\t" << pvec.ppprob2[j] << "\t";
        }
        pathfile << pvec.gcprob[j] << "\t" << pvec.transPGC[j] << endl;
    }
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

