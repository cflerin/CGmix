/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

void generateStates(const vector<vector<string> >& hapInfo, class hmmStates& st ) {
    int nrow = hapInfo.size() - 1 ;
    // G null states:
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        st.Xhap.push_back( hapInfo[j][0] );
        st.Xpop.push_back( hapInfo[j][1] );
        st.Xindx.push_back( j ); // check if this is right
        st.Ghap.push_back( "0" );
        st.Gpop.push_back( "0" );
        st.Gindx.push_back( j ); // check if this is right
    }
    // All combinations of X and G:
    for(int i=0; i < nrow; i++) {
        for (int j=0; j < nrow; j++) {
            // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
            st.Xhap.push_back( hapInfo[j][0] );
            st.Xpop.push_back( hapInfo[j][1] );
            st.Xindx.push_back( j );
            st.Ghap.push_back( hapInfo[i][0] );
            st.Gpop.push_back( hapInfo[i][1] );
            st.Gindx.push_back( i );
        }
    }
    // state concatenation:
    nrow = st.Xhap.size();
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = st.Xpop[i] + "-" + st.Xhap[i] + "-" + st.Gpop[i] + "-" + st.Ghap[i];
        st.states.push_back( tmp );
        // cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" <<
        //     st.Ghap[i] << "\t" << st.Gpop[i] << "\t" << st.Gindx[i] << "\t" << st.states[i] << endl;
    }
}

// void getXtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }
// 
// void getGtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }

void forward( 
        const vector<vector<int> >& sites,
        const vector<int>& dvec,
        const struct parameters& param,
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        vector<vector<double> > fwd ) {
    cout << fwd.size() << " " << fwd[0].size() << endl;
    // starting prob: 
    double e;
    for(int i=0; i < st.states.size(); i++) {
        if( sites[ st.Gindx[i] ][0] == obs[0] ) {
            e = emit.match;
        } else {
            e = emit.mismatch;
        }
        // cout << e << endl;
        if( st.Ghap[i] == "0" ) {
            fwd[i][0] = log( 1.0 / ( param.n1 + param.n2 ) * ( param.lam * ( param.n1 + param.n2 ) + param.gam * param.T) * e );
        } else {
            fwd[i][0] = log( 1.0 / ( param.n1 + param.n2 ) * ( param.gam * param.T ) / ( ( param.n1 + param.n2 ) * ( param.lam * (param.n1 + param.n2 ) + param.gam * param.T ) ) * e );
        }
    }
    int d, cnt;
    double lsum, tmp, trX, trG;
    double negInf = - std::numeric_limits<double>::infinity();
    for(int j=1; j < sites[0].size() ; j++ ) {
        d = dvec[j] - dvec[j-1];
        for(int t=0; t < st.states.size(); t++) {
            lsum = negInf;
            for(int f=0; f < st.states.size(); f++) {
                // trX = getXtrans( f, t, st, param );
                trX = 0.05;
                // trG = getGtrans( f, t, st, param );
                trG = 0.06;
                tmp = fwd[f][j-1] + log( trX * trG );
                if( tmp > negInf ) lsum = tmp + log( 1 + exp( lsum - tmp ) );
            } // end 'from' loop
            // emission prob:
            if( sites[ st.Gindx[t] ][j] == obs[j] ) {
                e = emit.match;
            } else {
                e = emit.mismatch;
            }
            fwd[t][j] = log( e ) + lsum;
        } // end 'to' loop
    } // end j site loop
}






