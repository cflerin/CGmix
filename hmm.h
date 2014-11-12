/*
 * hmm.h
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#ifndef HMM_H_
#define HMM_H_

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>

using namespace std;

struct parameters {
    int n1, n2, S;
    double T, 
           u1,
           rho,
           gam,
           lam,
           theta;
};
struct emissions {
    double match, mismatch;
};
class hmmStates 
{
    public:
        hmmStates(){};
        ~hmmStates(){};
        vector<string> Xpop, Xhap, Gpop, Ghap, states;
        vector<int> Xindx, Gindx;
};

void generateStates(const vector<vector<string> >& hapInfo, class hmmStates& st ) ;

void forward( 
        const vector<vector<int> >& sites,
        const vector<int>& dvec,
        const struct parameters& param,
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        vector<vector<double> > fwd );

#endif


