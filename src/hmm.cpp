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
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        st.Xhap.push_back( hapInfo.hN[j] );
        st.Xpop.push_back( hapInfo.hP[j] );
        st.Xindx.push_back( j ); // check if this is right
        st.Ghap.push_back( 0 );
        st.Gpop.push_back( 0 );
        st.Gindx.push_back( j ); // check if this is right
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
        }
    }
    // state concatenation:
    nrow = st.Xhap.size();
    //cout << "i\tXhap\tXpop\tXindx\tGhap\tGpop\tGindx\tstates" << endl;
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]) + "-" + std::to_string(st.Gpop[i]) + "-" + std::to_string(st.Ghap[i]);
        //tmp = st.Xpop[i] + "-" + st.Xhap[i] + "-" + st.Gpop[i] + "-" + st.Ghap[i];
        st.states.push_back( tmp );
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" <<
        //     st.Ghap[i] << "\t" << st.Gpop[i] << "\t" << st.Gindx[i] << "\t" << st.states[i] << endl;
    }
}

void generateXstates(const hapDef &hapInfo, hmmStates &st ) {
    int nrow = hapInfo.hN.size() ; // - 1 ;
    // Haplotype states only:
    for (int j=0; j < nrow; j++) {
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
        st.Xhap.push_back( hapInfo.hN[j] );
        st.Xpop.push_back( hapInfo.hP[j] );
        st.Xindx.push_back( j );
        st.Ghap.push_back( 0 ); // fill this vector with null GC states
    }
    // state concatenation:
    nrow = st.Xhap.size();
    //cout << "i\tXhap\tXpop\tXindx\tGhap\tstates" << endl;
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]);
        st.states.push_back( tmp );
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" << st.states[i] << endl;
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" << st.Ghap[i] << "\t" << st.states[i] << endl;
    }
}

void getXtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, double &trX ) {
    double u, rho;
    int nm;
    if( st.Xpop[to] == 1 ) {
        u = p.u1;
        nm = p.n1;
        //rho = 4.0 * p.Ne1 * (r/d);
        rho = ( 4.0 * p.Ne1 * (r/d) ) / nm;
    }
    if( st.Xpop[to] == 2 ) {
        u = 1.0 - p.u1;
        nm = p.n2;
        rho = ( 4.0 * p.Ne2 * (r/d) ) / nm;
    }
    // cout << "From " << st.Xhap[from] << " To " << st.Xhap[to] ;
    // cout << " | " << u << " " << nm << " ";
    if( st.Xpop[from] != st.Xpop[to] ) { // ancestry and haplotype switch:
        // cout << " anc sw! ";
        // trX = ( 1.0 - exp( -d * p.rho * p.T ) ) * u/nm;
        trX = ( 1.0 - exp( -r * p.T ) ) * u/nm;
    } else { // no ancestry switch
        if( st.Xhap[from] != st.Xhap[to] ) { // hap switch
            // cout << " hap sw! ";
            // trX = exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
            //       ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
            trX = exp( -r * p.T ) * ( 1.0 - exp( -r * rho ) ) / nm + 
                  ( 1.0 - exp( -r * p.T )) * u / nm ;
        }
        if( st.Xhap[from] == st.Xhap[to] ) { // no hap switch
            // cout << " cont! ";
            // trX = exp( -d * p.rho * p.T ) * exp( -d * p.rho ) + 
            //       exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
            //       ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
            trX = exp( -r * p.T ) * exp( -r * rho ) + 
                  exp( -r * p.T ) * ( 1.0 - exp( -r * rho ) ) / nm + 
                  ( 1.0 - exp( -r * p.T )) * u / nm ;
        }
    }
    trX = log(trX);
}

void getsprob( 
        const vector<int> &sites0, 
        const parameters &p, 
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob ) {
    double e;
    double g = 1/100000.0 * p.f; // 1cM/Mb * f = 1Morgan/10000kb
    double gam1 = 4 * p.Ne1 * g;
    double gam2 = 4 * p.Ne2 * g;
    double sprob_G0 = log( 1.0 / ( p.n1 + p.n2 ) * 
            ( p.lam * (p.n1+p.n2) ) / ( p.lam * (p.n1+p.n1) + (gam1/p.n1+gam2/p.n2)) );
    double sprob_G1 = log( 1.0 / ( p.n1 + p.n2 ) * 
            (gam1/p.n1+gam2/p.n2) / ( (p.n1+p.n2) * ( p.lam * (p.n1+p.n2) + (gam1/p.n1+gam2/p.n2) )) );
    for(int i=0; i < st.states.size(); i++) {
        if( sites0[ st.Gindx[i] ] == obs[0] ) {
            if( st.Xpop[i] == 1 ) { e = p.theta1_match; }
            else if( st.Xpop[i] == 2 ) { e = p.theta2_match; }
        } else {
            if( st.Xpop[i] == 1 ) { e = p.theta1_mismatch; }
            else if( st.Xpop[i] == 2 ) { e = p.theta2_mismatch; }
        }
        if( st.Ghap[i] == 0 ) {
            sprob[i] =  sprob_G0 + e;
        } else {
            sprob[i] =  sprob_G1 + e;
        }
    }
}

void getsprobX(
        const vector<int> &sites0,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob ) {
    double e;
    double sprobX = log( 1.0 / ( p.n1 + p.n2 ) );
    for(int i=0; i < st.states.size(); i++) {
        if( sites0[ st.Xindx[i] ] == obs[0] ) {
            if( st.Xpop[i] == 1 ) { e = p.theta1_match; }
            else if( st.Xpop[i] == 2 ) { e = p.theta2_match; }
        } else {
            if( st.Xpop[i] == 1 ) { e = p.theta1_mismatch; }
            else if( st.Xpop[i] == 2 ) { e = p.theta2_mismatch; }
        }
        sprob[i] = sprobX + e;
    }
}

double lookupXtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trXbin ) {
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
    if( trXbin[ type ] == 99 ) { // this transition hasn't been set yet:
        getXtrans( to, from, d, r, st, p, trX);
        trXbin[ type ] = trX;
        return( trX );
    } else { // lookup and return a value already seen:
        trX = trXbin[ type ];
        return( trX );
    }
}

inline double lookupGtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trGbin ) {
    // first pass at this site: calculate all reusable values:
    if( trGbin[10] == 99 ) {
        double int1, int2p1, int2p2, gam1, gam2, tmp;
        double pfrd4 = p.f * 4 * r / d;
        gam1 = pfrd4 * p.Ne1;
        gam2 = pfrd4 * p.Ne2;
        //cout << "r=" << r << "\td=" << d << "\tgam1=" << gam1 << "\tgam2=" << gam2 << endl;
        //gam1 = p.f * ( 4.0 * p.Ne1 * (r/d) );
        //gam2 = p.f * ( 4.0 * p.Ne2 * (r/d) );
        //
        double rpT = r*p.T;
        double plampn1pn2 = p.lam*p.n1*p.n2;
        double tmp1 = exp( -d * (p.lam+(gam1/p.n1 + gam2/p.n2)*rpT));
        double tmp2 = exp( -p.lam * d );
        double tmp3 = exp( -rpT * d * (gam1/p.n1 + gam2/p.n2));
        double gam2pn1gam1pn2 = (gam2*p.n1+gam1*p.n2);
        int1 = ( ( 1.0 - tmp1) * plampn1pn2 ) / ( plampn1pn2 + gam2pn1gam1pn2*rpT );
        int2p1 = ( ( tmp1 - 1.0) * p.lam*p.n2*p.u1 ) /
                 ( plampn1pn2 + gam2pn1gam1pn2*rpT ) +
                 ( p.u1 - tmp2 * p.u1 ) / p.n1;
        int2p2 = ( ( tmp1 - 1.0) * p.lam*p.n1*p.u2 ) /
                 ( plampn1pn2 + gam2pn1gam1pn2*rpT ) +
                 ( p.u2 - tmp2 * p.u2 ) / p.n2;
        // 0 //
        tmp = log( exp( -p.lam * d ) * exp( -r * d * p.T * (gam1/p.n1 + gam2/p.n2) ) + int1 );
        trGbin[0] = tmp;
        trGbin[5] = tmp;
        // 1 //
        trGbin[1] = log( tmp2 * ( 1.0 - tmp3) * p.u1/p.n1 + int2p1 );
        // 2 //
        tmp = log( int1 );
        trGbin[2] = tmp;
        trGbin[7] = tmp;
        // 3 //
        trGbin[3] = log( tmp2 + int2p1 );
        // 4 //
        trGbin[4] = log( int2p1 );
        /////// population 2:
        // 5 //
        // 6 //
        trGbin[6] = log( tmp2 * ( 1.0 - tmp3) * p.u2/p.n2 + int2p2 );
        // 7 //
        // 8 //
        trGbin[8] = log( tmp2 + int2p2 );
        // 9 //
        trGbin[9] = log( int2p2 );
        /////// flag:
        trGbin[10] = 1.0;
        // prob sums check:
        // for(int z=0; z<trGbin.size(); z++) { cout << exp(trGbin[z]) << endl; }
        // cout << "trGbinsum0=" << exp(trGbin[0]) + exp(trGbin[1])*p.n1 + exp(trGbin[6])*p.n2 << endl;
        // cout << "trGbinsum1=" << exp(trGbin[2]) + exp(trGbin[3]) + exp(trGbin[4])*(p.n1-1) + exp(trGbin[9])*(p.n2) << endl;
        // cout << "trGbinsum2=" << exp(trGbin[2]) + exp(trGbin[8]) + exp(trGbin[9])*(p.n2-1) + exp(trGbin[4])*(p.n1) << endl;
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
    vector<double> trXbin(6,99);
    vector<double> trGbin(11,99);
    vector<int> siteIndx;
    vector<double> tmp;
    //cout << "pswitch"<<pswitch.size() << endl;
    for(int j=1; j < sites.size() ; j++ ) {
        //cout << "j = " << j << "\t" << "pswitch=" << pswitch[j] << " Nstates=" << fwd[j].size() << endl;
        //d = dvec[j] - dvec[j-1];
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
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
                    //cout << "\ttrx" << endl;
                    tmp[f] = fwd[j-1][f] + trX;
                }
            } // end 'from' loop
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    vector<double> trXbin(6,99);
    vector<double> trGbin(11,99);
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
                    if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                    else if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
                } else {
                    if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                    else if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    vector<double> trXbin(6,99);
    vector<double> trGbin(11,99);
    vector<int> siteIndx;
    if( p.mode == 1 ) {
        siteIndx = st.Gindx;
    } else {
        siteIndx = st.Xindx;
    }
    for(int j=1; j < sites.size() ; j++ ) {
        //d = dvec[j] - dvec[j-1];
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
                if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    pvec.vprob[j] = vit[j][maxix];
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
        pvec.vprob[j] = vit[j][maxix];
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
    vector<double> trXbin(6,99);
    vector<double> trGbin(11,99);
    vector<int> siteIndx;
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
            // emission prob:
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                else if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    pvec.vprob[j] = vit[j][maxix];
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
        pvec.vprob[j] = vit[j][maxix];
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
    gcprobPop1.resize( param.S, 0.0 );
    gcprobPop2.resize( param.S, 0.0 );
    gcprobXPop.resize( param.S, 0.0 );
    transPGC.resize( param.S, 0.0 );
//
    if( param.mode == 2 ) {
        pppath2.resize( param.S );
        ppprob2.resize( param.S, 0.0 );
        vpath2.resize( param.S );
        vprob2.resize( param.S, 0.0 );
    }
}

void pathOutput( pathVec &pvec, hmmStates &st, positions &pos, vector<vector<double> > &pprob, parameters &param, ofstream &pathfile ) {
    if( ( param.mode == 1 ) || ( param.mode == 2 ) ) {
        double pXi, pGk;
        // calculate total probability of a cross-population gene conversion:
        for(int j=0; j < pprob.size(); j++ ) {
            for(int i=0; i < pprob[j].size(); i++ ) {
                //if( st.Ghap[i] == 0 ) { continue; }
                if( st.Xpop[i] == 1 ) {
                    pXi = pprob[j][i];
                    for(int k=0; k < pprob[j].size(); k++ ) {
                        if( st.Gpop[k] == 2 ) {
                            pGk = pprob[j][k];
                            pvec.transPGC[j] += pXi * pGk;
                        }
                    }
                }
            }
            for(int i=0; i < pprob[j].size(); i++ ) {
                //if( st.Ghap[i] == 0 ) { continue; }
                if( st.Xpop[i] == 2 ) {
                    pXi = pprob[j][i];
                    for(int k=0; k < pprob[j].size(); k++ ) {
                        if( st.Gpop[k] == 1 ) {
                            pGk = pprob[j][k];
                            pvec.transPGC[j] += pXi * pGk;
                        }
                    }
                }
            }
        }
        for(int j=0; j < pprob.size(); j++ ) {
            for(int i=0; i < pprob[j].size(); i++ ) {
                if( st.Ghap[i] != 0 ) { // gene conversion state
                    pvec.gcprob[j] += pprob[j][i];
                    if( st.Gpop[i] == 1 ) {
                        pvec.gcprobPop1[j] += pprob[j][i];
                    } else if( st.Gpop[i] == 2 ) {
                        pvec.gcprobPop2[j] += pprob[j][i];
                    }
                    if( ( st.Xpop[i]==1 && st.Gpop[i]==2 ) || ( st.Xpop[i]==2 && st.Gpop[i]==1 ) ) {
                        // cross-population gene conversion
                        pvec.gcprobXPop[j] += pprob[j][i];
                    }
                }
            }
        }
    } // endif
    // output path and probabilites:
    pathfile << "site\tpos\tVpath\tVpathProb\tPpath\tPpathProb\tPswitch\t";
    if( param.mode == 2 ) 
        pathfile << "Vpath2\tVpathProb2\tPpath2\tPpathProb2\t";
    pathfile << "GCprob\tGCprobP1\tGCprobP2\tGCprobXPop0\tGCprobXPop" << endl;
    int intpos;
    for(int j=0; j < pprob.size(); j++ ) {
        intpos = static_cast<int>( pos.pos[j]*1000+0.5 );
        pathfile << j << "\t" << intpos << "\t" << pvec.vpath[j] << "\t" << exp(pvec.vprob[j]) << "\t" << pvec.pppath[j] << "\t" << pvec.ppprob[j] << "\t" << pvec.pswitch[j] << "\t";
        if( param.mode == 2 ) 
            pathfile << pvec.vpath2[j] << "\t" << exp(pvec.vprob2[j]) << "\t" << pvec.pppath2[j] << "\t" << pvec.ppprob2[j] << "\t";
        pathfile << pvec.gcprob[j] << "\t" << pvec.gcprobPop1[j] << "\t" << pvec.gcprobPop2[j] << "\t" << pvec.gcprobXPop[j] << "\t" << pvec.transPGC[j] << endl;
    }
}

