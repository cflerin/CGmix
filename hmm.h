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
#include <iomanip>
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

void generateStates(const vector<vector<string> >& hapInfo, class hmmStates& st );

void getXtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trX );
void getGtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trG );

void forward( 
        const vector<vector<int> >& sites,
        const vector<int>& dvec,
        const struct parameters& p,
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        vector<vector<double> >& fwd );

void backward(
		const vector<vector<int> >& sites,
		const vector<int>& dvec,
		const struct parameters& p,
		const struct emissions& emit,
		const class hmmStates& st,
		const vector<int>& obs,
		vector<vector<double> >& bwd );

void printMat( const vector<vector<double> >& mat );
void logSumExp( const vector<double>& vec, double& lse );
void postDecode(
		const vector<vector<double> >& fwd,
		const vector<vector<double> >& bwd,
		vector<vector<double> >& pprob
		);

void viterbi(
		const vector<vector<int> >& sites,
		const vector<int>& dvec,
		const struct parameters& p,
		const struct emissions& emit,
		const class hmmStates& st,
		const vector<int>& obs,
		vector<vector<double> >& vit,
		vector<string>& vpath );

void max( const vector<double>& vec, double maxelement );


#endif


