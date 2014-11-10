/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;


void generateStates(const vector<vector<string> >& hapInfo, vector<vector<string> > &st ) {
    int nrow = hapInfo.size() - 2 ;
    // cout << nrow << endl;
    vector<string> Xpop;
    vector<string> Xhap;
    vector<string> Gpop;
    vector<string> Ghap;
    // G null states:
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        Xhap.push_back( hapInfo[j][0] );
        Xpop.push_back( hapInfo[j][1] );
        Ghap.push_back( "0" );
        Gpop.push_back( "0" );
    }
    // All combinations of X and G:
    for(int i=0; i < nrow; i++) {
        for (int j=0; j < nrow; j++) {
            // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
            Xhap.push_back( hapInfo[j][0] );
            Xpop.push_back( hapInfo[j][1] );
            Ghap.push_back( hapInfo[i][0] );
            Gpop.push_back( hapInfo[i][1] );
        }
    }
    st.push_back( Xpop );
    st.push_back( Xhap );
    st.push_back( Gpop );
    st.push_back( Ghap );
    // state concatenation:
    vector<string> states;
    nrow = st[0].size();
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = Xpop[i] + Xhap[i] + Gpop[i] + Ghap[i];
        states.push_back( tmp );
    }
    st.push_back( states );
    // cout << states.size() << "\t" << Xpop.size() << endl;
}

// void getXtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }
// 
// void getGtrans(int to, int from, vector<vector<string> > st, vector<double> param ) {
// }






