/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

void generateStates(const vector<vector<string> >& hapInfo, struct states& st ) {
    int nrow = hapInfo.size() - 2 ;
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
    }
}

// void getXtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }
// 
// void getGtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }

void forward( 
        const vector<vector<int> >& sites,
        const struct parameters& param,
        const struct emissions& emit,
        const struct states& st,
        const vector<int>& obs,
        vector<vector<double> > fwd ) {
    // starting prob: 
    for(int i=0; i < st.states.size(); i++) {
        double e;
        if( sites[ st.Gindx[i] ][0] == obs[0] ) {
            e = emit.match;
        } else {
            e = emit.mismatch;
        }
        // cout << e << endl;
    }
}






