/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

void generateStates(const hapDef &hapInfo, hmmStates &st ) {
    int nrow = hapInfo.hN.size() - 1 ;
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
    int nrow = hapInfo.hN.size() - 1 ;
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
    //cout << "i\tXhap\tXpop\tXindx\tstates" << endl;
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]);
        st.states.push_back( tmp );
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" << st.states[i] << endl;
        //cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" << st.Ghap[i] << st.states[i] << endl;
    }
}

void getXtrans(const int &to, const int &from, const double &r, const hmmStates &st, const parameters &p, double &trX ) {
    double u, rho;
    int nm;
    if( st.Xpop[to] == 1 ) {
        u = p.u1;
        nm = p.n1;
        rho = p.rho1;
    }
    if( st.Xpop[to] == 2 ) {
        u = 1.0 - p.u1;
        nm = p.n2;
        rho = p.rho2;
    }
    // cout << "From " << st.Xhap[from] << " To " << st.Xhap[to] ;
    // cout << " | " << u << " " << nm << " ";
    if( st.Xpop[from] != st.Xpop[to] ) { // ancestry and haplotype switch:
        // cout << " anc sw! ";
        //trX = ( 1.0 - exp( -d * p.rho * p.T ) ) * u/nm;
        trX = ( 1.0 - exp( -r * p.T ) ) * u/nm;
    } else { // no ancestry switch
        if( st.Xhap[from] != st.Xhap[to] ) { // hap switch
            // cout << " hap sw! ";
            //trX = exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
            //      ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
            trX = exp( -r * p.T ) * ( 1.0 - exp( -r * rho ) ) / nm + 
                  ( 1.0 - exp( -r * p.T )) * u / nm ;
        }
        if( st.Xhap[from] == st.Xhap[to] ) { // no hap switch
            // cout << " cont! ";
            //trX = exp( -d * p.rho * p.T ) * exp( -d * p.rho ) + 
            //      exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
            //      ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
            trX = exp( -r * p.T ) * exp( -r * rho ) + 
                exp( -r * p.T ) * ( 1.0 - exp( -r * rho ) ) / nm + 
                ( 1.0 - exp( -r * p.T )) * u / nm ;

        }
    }
    trX = log(trX);
}

void getGtrans(const int &to, const int &from, const int &d, const hmmStates &st, const parameters &p, double &trG ) {
    double u, a1, a2, a3, a4, a5, gam;
    int nm;
    if( st.Gpop[to] == 1 ) {
        u = p.u1;
        nm = p.n1;
        gam = p.gam1;
    }
    if( st.Gpop[to] == 2 ) {
        u = 1.0 - p.u1;
        nm = p.n2;
        gam = p.gam2;
    }
    // 4 //
    if( ( st.Ghap[from] == st.Ghap[to] ) && ( st.Ghap[from] != 0 ) ) {
        // 5: Pr(G[j+1]=g' | G[j] = g):
        a5 = ( p.lam * u * (p.n1 + p.n2) * ( exp( d*(-gam * p.T / (p.n1 + p.n2) - p.lam)) - 1.0)) / ((p.n1 + p.n2) * nm * p.lam + gam * p.T * nm) + (u - u * exp(-d * p.lam)) / nm ;
        // 4: Pr(G[j+1]=g | G[j] = g):
        a4 = exp( -p.lam * d ) + a5 ;
        trG = a4;     
    }
    // 1 //
    if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] == 0 ) ) {
        // 3: Pr(G[j+1]=0 | G[j] = g):
        a3 = p.lam * (p.n1 + p.n2) * ( 1.0 - exp( -d * (gam * p.T + p.lam * (p.n1 + p.n2) ) / (p.n1 + p.n2) ) ) / (gam * p.T + p.lam * (p.n1 + p.n2) ) ;
        // 1: Pr(G[j+1]=0 | G[j] = 0):
        a1 = exp( -p.lam * d ) * exp( -gam * p.T * d / (p.n1 + p.n2)) + a3 ;
        trG = a1;
    }
    // 2 //
    if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] != 0 ) ) {
        // 5: Pr(G[j+1]=g' | G[j] = g):
        a5 = ( p.lam * u * (p.n1 + p.n2) * ( exp( d*(-gam * p.T / (p.n1 + p.n2) - p.lam)) - 1)) / ((p.n1 + p.n2) * nm * p.lam + gam * p.T * nm) + (u - u * exp(-d * p.lam)) / nm ;
        // 2: Pr(G[j+1]=g | G[j] = 0):
        a2 = exp( -p.lam * d ) * ( 1.0 - exp( -gam * p.T * d / (p.n1 + p.n2))) * u / nm + a5 ;
        trG = a2;
    }
    // 3 //
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] == 0 ) ) {
        // 3: Pr(G[j+1]=0 | G[j] = g):
        a3 = p.lam * (p.n1 + p.n2) * ( 1.0 - exp( -d * (gam * p.T + p.lam * (p.n1 + p.n2) ) / (p.n1 + p.n2) ) ) / (gam * p.T + p.lam * (p.n1 + p.n2) ) ;
        trG = a3;
    }
    // 5 //
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] != 0 ) && ( st.Ghap[from] != st.Ghap[to] ) ) {
        // 5: Pr(G[j+1]=g' | G[j] = g):
        a5 = ( p.lam * u * (p.n1 + p.n2) * ( exp( d*(-gam * p.T / (p.n1 + p.n2) - p.lam)) - 1.0)) / ((p.n1 + p.n2) * nm * p.lam + gam * p.T * nm) + (u - u * exp(-d * p.lam)) / nm ;
        trG = a5;
    }
    trG = log(trG);
}

void getsprob( 
        const vector<int> &sites0, 
        const parameters &p, 
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob ) {
    double e, gam, sprob_G0, sprob_G1;
    // double sprob_G0 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * (p.n1+p.n1) + p.gam * p.T) );
    // double sprob_G1 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.gam * p.T ) / ( ( p.n1 + p.n2 ) * ( p.lam * (p.n1 + p.n2 ) + p.gam * p.T ) ) );
    for(int i=0; i < st.states.size(); i++) {
        // population-specific gamma:
        if( st.Gpop[i] == 1 ) { gam = p.gam1; }
        if( st.Gpop[i] == 2 ) { gam = p.gam2; }
        // emissions:
        if( sites0[ st.Gindx[i] ] == obs[0] ) {
            if( st.Xpop[i] == 1 ) { e = p.theta1_match; }
            if( st.Xpop[i] == 2 ) { e = p.theta2_match; }
        } else {
            if( st.Xpop[i] == 1 ) { e = p.theta1_mismatch; }
            if( st.Xpop[i] == 2 ) { e = p.theta2_mismatch; }
        }
        //
        if( st.Ghap[i] == 0 ) {
            sprob_G0 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * ( p.n1 + p.n1 ) + gam * p.T ) );
            sprob[i] =  sprob_G0 + e;
        } else {
            sprob_G1 = log( 1.0 / ( p.n1 + p.n2 ) * ( gam * p.T ) / ( ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) + gam * p.T ) ) );
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
            if( st.Xpop[i] == 2 ) { e = p.theta2_match; }
        } else {
            if( st.Xpop[i] == 1 ) { e = p.theta1_mismatch; }
            if( st.Xpop[i] == 2 ) { e = p.theta2_mismatch; }
        }
        sprob[i] = sprobX + e;
    }
}

double lookupXtrans(const int &to, const int &from, const double &r, const hmmStates &st, const parameters &p, vector<double> &trXbin ) {
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
        getXtrans( to, from, r, st, p, trX);
        trXbin[ type ] = trX;
        return( trX );
    } else { // lookup and return a value already seen:
        trX = trXbin[ type ];
        return( trX );
    }
}

double lookupGtrans(const int &to, const int &from, const int &d, const hmmStates &st, const parameters &p, vector<double> &trGbin ) {
    double trG;
    int type;
    if( ( st.Ghap[from] == st.Ghap[to] ) && ( st.Ghap[from] != 0 ) )
        type = 3;
    else if( ( st.Ghap[from] == st.Ghap[to] ) && ( st.Ghap[from] == 0 ) )
        type = 0; // overwrites above
    else if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] != 0 ) )
        type = 1;
    else if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] == 0 ) )
        type = 2;
    else //if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] != 0 ) && ( st.Ghap[from] != st.Ghap[to] ) )
        type = 4;
    if( st.Gpop[to] == 2 ) // shift type index by 3 if transitioning to a haplotype in pop #2.
        type += 5;
    //cout << " from" << st.Ghap[from] << " to" << st.Ghap[to];
    //cout << " type" << type << " ";
    if( trGbin[ type ] == 99 ) { // this transition hasn't been set yet:
        getGtrans( to, from, d, st, p, trG);
        trGbin[ type ] = trG;
        return( trG );
    } else { // lookup and return a value already seen:
        trG = trGbin[ type ];
        return( trG );
    }
}

void forward( 
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &fwd ) {
    int d;
    double lsum, trX, trG, e, r;
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
    vector<int> siteIndx;
    vector<double> tmp( st.states.size() );
    if( p.mode == 1 ) { // full model
        siteIndx = st.Gindx;
    } else { // reduced
        siteIndx = st.Xindx;
    }
    // starting prob: 
    fwd[0] = sprob;
    for(int j=1; j < sites.size() ; j++ ) {
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        for(int t=0; t < st.states.size(); t++) {
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int f=0; f < st.states.size(); f++) {
                trX = lookupXtrans( t, f, r, st, p, trXbin );
                if( p.mode == 1 ) {
                    trG = lookupGtrans( t, f, d, st, p, trGbin );
                    tmp[f] = fwd[j-1][f] + trX + trG;
                } else {
                    tmp[f] = fwd[j-1][f] + trX;
                }
                //cout << "t= " << st.states[t] << " f= " << st.states[f] << " fwd[j-1][f]= " << fwd[j-1][f] << " trX= " << trX << endl;
            } // end 'from' loop
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
            }
            logSumExp( tmp, lsum );
            fwd[j][t] = e + lsum;
            // cout << t << endl;
        } // end 'to' loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // end j site loop
}

void forward2(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const hmmStates &st2,
        const vector<int> &obs,
        const vector<double> &sprob,
        const vector<int> &pswitch,
        vector<vector<double> > &fwd ) {
    // starting prob:
    fwd[0] = sprob;
    int d;
    double lsum, trX, trG, e, r;
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
    vector<int> siteIndx;
    vector<double> tmp;
    //cout << "pswitch"<<pswitch.size() << endl;
    for(int j=1; j < sites.size() ; j++ ) {
        //cout << "j = " << j << "\t" << "pswitch=" << pswitch[j] << " Nstates=" << fwd[j].size() << endl;
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        if( pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Xindx;
        } else {
            siteIndx = st2.Gindx;
        }
        //cout << "fwd[j].size() = " << fwd[j].size() << endl;
        for(int t=0; t < fwd[j].size(); t++) {
            tmp.resize( fwd[j-1].size(), 0.0 );
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int f=0; f < fwd[j-1].size(); f++) {
                //cout << "t=" << t << " " << st2.states[t] << "\tf=" << f << " " << st2.states[f] ;
                trX = lookupXtrans( t, f, r, st2, p, trXbin );
                if( ( pswitch[j] == 1 ) || ( pswitch[j-1] == 1 ) ) {
                    //cout << "\ttrx + trG" << endl;
                    trG = lookupGtrans( t, f, d, st2, p, trGbin );
                    tmp[f] = fwd[j-1][f] + trX + trG;
                } else {
                    //cout << "\ttrx" << endl;
                    tmp[f] = fwd[j-1][f] + trX;
                }
            } // end 'from' loop
            if( sites[j][ siteIndx[t] ] == obs[j] ) {
                if( st2.Xpop[t] == 1 ) { e = p.theta1_match; }
                if( st2.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st2.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                if( st2.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
            }
            logSumExp( tmp, lsum );
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
        vector<vector<double> > &bwd,
        double &Pxb) {
    double e, lsum, trX, trG, r;
    int d;
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
    vector<int> siteIndx;
    vector<double> tmp( st.states.size() );
    if( p.mode == 1 ) { // full model
        siteIndx = st.Gindx;
    } else { // reduced
        siteIndx = st.Xindx;
    }
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        for(int f=0; f < st.states.size(); f++ ) {
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int t=0; t < st.states.size(); t++ ) {
                trX = lookupXtrans( t, f, r, st, p, trXbin);
                //cout << "j= " << j << "\tf= " << f << "\tt= " << t;
                if( sites[j+1][ siteIndx[t] ] == obs[j+1] ) {
                    if( st.Xpop[t] == 1 ) { e = p.theta1_match; }
                    if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
                } else {
                    if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                    if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
                }
                if( p.mode == 1 ) {
                    trG = lookupGtrans( t, f, d, st, p, trGbin);
                    //cout << "\ttrG= " << trG;
                    tmp[t] = bwd[j+1][t] + trX + trG + e;
                } else {
                    tmp[t] = bwd[j+1][t] + trX + e;
                    //cout << "\tbwd[j+1][t]= " << bwd[j+1][t] << "\ttrX= " << trX << "\te= " << e << "\ttmp[t]= " << tmp[t] ;
                }
                //cout << endl;
            } // t loop
            logSumExp( tmp, lsum );
            //cout << "j= " << j << "\tf= " << f << "\tlsum= " << lsum << endl << endl;
            bwd[j][f] = lsum;
        } // f loop
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        std::fill( trGbin.begin(), trGbin.end(), 99 );
    } // j loop
    // termination
    tmp.resize( bwd[0].size(), 0.0 );
    std::fill( tmp.begin(), tmp.end(), 0.0 );
    for(int f=0; f < bwd[0].size(); f++ ) {
        tmp[f] = bwd[0][f] + sprob[f]; // e already included within sprob...
    }
    Pxb = 0.0;
    for(int i=0; i<tmp.size(); i++ ) {
        Pxb += exp( tmp[i] );
    }
    //logSumExp(tmp,Pxb);
}

void backward2(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const hmmStates &st2,
        const vector<int> &obs,
        const vector<double> &sprob,
        const vector<int> &pswitch,
        vector<vector<double> > &bwd,
        double &Pxb ) {
    double e, lsum, trX, trG, r;
    int d;
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
    vector<int> siteIndx;
    vector<double> tmp;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        //cout << "j = " << j << "\t" << "pswitch=" << pswitch[j] << " Nstates=" << bwd[j].size() << endl;
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        if( pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Xindx;
        } else {
            siteIndx = st2.Gindx;
        }
        for(int f=0; f < bwd[j].size(); f++ ) {
            tmp.resize( bwd[j+1].size(), 0.0 );
            std::fill( tmp.begin(), tmp.end(), 0.0 );
            for(int t=0; t < bwd[j+1].size(); t++ ) {
                //cout << "t=" << t << " " << st2.states[t] << "\tf=" << f << " " << st2.states[f] ;
                trX = lookupXtrans( t, f, r, st2, p, trXbin);
                if( sites[j+1][ siteIndx[t] ] == obs[j+1] ) {
                    if( st2.Xpop[t] == 1 ) { e = p.theta1_match; }
                    if( st2.Xpop[t] == 2 ) { e = p.theta2_match; }
                } else {
                    if( st2.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                    if( st2.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
                }
                if( ( pswitch[j] == 1 ) || ( pswitch[j+1] == 1 ) ) {
                    //cout << "\ttrx + trG" << endl;
                    trG = lookupGtrans( t, f, d, st2, p, trGbin);
                    tmp[t] = bwd[j+1][t] + trX + trG + e ;
                } else {
                    //cout << "\ttrx" << endl;
                    tmp[t] = bwd[j+1][t] + trX + e ;
                }
            } // t loop
            logSumExp( tmp, lsum );
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
    Pxb = 0.0;
    for(int i=0; i<tmp.size(); i++ ) {
        Pxb += exp( tmp[i] );
    }
    //logSumExp(tmp,Pxb);
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

void logSumExp( const vector<double> &vec, double &lse ) {
    double max = - std::numeric_limits<double>::max();
    for(int i=0; i < vec.size(); i++ ) {
        if( vec[i] > max )
            max = vec[i];
    }
    //cout << setprecision(64) << "max =" << max << endl;
    double sum = 0.0;
    vector<double> vec2(vec);
    sort(vec2.begin(), vec2.end());
    for(int i=0; i < vec2.size(); i++ ) {
        sum += exp( vec2[i] - max );
        //cout << i << " " << sum << endl;
    }
    lse = max + log(sum) ;
}

void postDecode(
        const vector<vector<double> > &fwd,
        const vector<vector<double> > &bwd,
        const hmmStates &st,
        const double &Px,
        vector<vector<double> > &pprob,
        vector<string> &pppath,
        vector<double> &ppprob,
        vector<int> &pswitch,
        ofstream &logfile
        ) {
    double negInf = - std::numeric_limits<double>::infinity();
    // decoding:
    for(int j=0; j < pprob.size(); j++ ) {
        // cout << "j= " << j << "\tfwd[j]= " << fwd[j].size() << "\tbwd[j]= " << bwd[j].size() << "\tpprob[j]= " << pprob[j].size() << "\tstLen= " << stLen << endl;
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
        if( ( j > 0 ) && ( pppath[j] != pppath[j-1] ) )
            pswitch[j] = 1;
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
        vector<string> &vpath,
        vector<double> &vprob ) {
    // starting prob:
    vit[0] = sprob;
    // recursion:
    int d;
    double trX, trG, vmax, tmp, e, r;
    double negInf = - std::numeric_limits<double>::infinity();
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
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
                trX = lookupXtrans( t, f, r, st, p, trXbin);
                if( p.mode == 1 ) {
                    trG = lookupGtrans( t, f, d, st, p, trGbin);
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
                if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    vpath[j] = st.states[maxix];
    vprob[j] = vit[j][maxix];
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < st.states.size(); f++) {
            trX = lookupXtrans( t, f, r, st, p, trXbin);
            if( p.mode == 1 ) {
                trG = lookupGtrans( t, f, d, st, p, trGbin);
                tmp = vit[j][f] + trX + trG;
            } else {
                tmp = vit[j][f] + trX;
            }
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        vpath[j] = st.states[ maxix ];
        vprob[j] = vit[j][maxix];
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
        const vector<int> &pswitch,
        vector<vector<double> > &vit,
        vector<string> &vpath,
        vector<double> &vprob ) {
    // starting prob, same as fwd:
    vit[0] = sprob;
    // recursion:
    int d;
    double trX, trG, vmax, tmp, e, r;
    double negInf = - std::numeric_limits<double>::infinity();
    vector<double> trXbin(6,99);
    vector<double> trGbin(10,99);
    vector<int> siteIndx;
    for(int j=1; j < sites.size() ; j++ ) {
        d = pos.pos[j] - pos.pos[j-1];
        r = pos.cM[j] - pos.cM[j-1];
        if( pswitch[j] == 0 ) { // choose sites to match for emission
            siteIndx = st.Gindx;
        } else {
            siteIndx = st.Xindx;
        }
        for(int t=0; t < vit[j].size(); t++ ) {
            vmax = negInf;
            for(int f=0; f < vit[j-1].size(); f++ ) {
                trX = lookupXtrans( t, f, r, st, p, trXbin);
                if( ( pswitch[j] == 1 ) | ( pswitch[j-1] == 1 ) ) {
                    trG = lookupGtrans( t, f, d, st, p, trGbin);
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
                if( st.Xpop[t] == 2 ) { e = p.theta2_match; }
            } else {
                if( st.Xpop[t] == 1 ) { e = p.theta1_mismatch; }
                if( st.Xpop[t] == 2 ) { e = p.theta2_mismatch; }
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
    vpath[j] = st.states[maxix];
    vprob[j] = vit[j][maxix];
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < vit[j].size(); f++) {
            trX = lookupXtrans( t, f, r, st, p, trXbin);
            if( pswitch[j] == 1 ) {
                trG = lookupGtrans( t, f, d, st, p, trGbin);
                tmp = vit[j][f] + trX + trG;
            } else {
                tmp = vit[j][f] + trX;
            }
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        vpath[j] = st.states[ maxix ];
        vprob[j] = vit[j][maxix];
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
