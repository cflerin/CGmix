/*
 * hmm.h
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#ifndef HMM_H_
#define HMM_H_

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
// #include <sstream>

#include "readFiles.h"
#include "parameters.h"

using namespace std;
extern ofstream matfile;

class hmmStates 
{
    public:
        hmmStates(){};
        ~hmmStates(){};
        vector<string> states;
        vector<int> Xpop, Xhap, Gpop, Ghap, Xindx, Gindx;
};
class pathVec
{
    public:
        pathVec( parameters param );
        ~pathVec(){};
        vector<int> pswitch;
        vector<string> vpath, pppath, pppath2;
        vector<double> vprob, gcprob, gcprobPop1, gcprobPop2, ppprob, gcprobXPop, transPGC;
        //
        vector<double> ppprob2, vprob2;
        vector<int> pswitch2;
        vector<string> vpath2;
};

void generateStates(const hapDef &hapInfo, hmmStates &st );
void generateXstates(const hapDef &hapInfo, hmmStates &st );

void getXtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, double &trX );

void getsprob( 
        const vector<int> &sites0, 
        const parameters &p, 
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob );
void getsprobX( 
        const vector<int> &sites0, 
        const parameters &p, 
        const hmmStates &st,
        const vector<int> &obs,
        vector<double> &sprob );

double lookupXtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trXbin );
inline double lookupGtrans(const int &to, const int &from, const double &d, const double &r, const hmmStates &st, const parameters &p, vector<double> &trGbin );

void forward( 
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &fwd );
void forward2( 
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const hmmStates &st2,
        const vector<int> &obs,
        const vector<double> &sprob,
        const vector<int> &pswitch,
        vector<vector<double> > &fwd );

void backward(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &bwd,
        double &Pxb );
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
        double &Pxb);

void printMat( const vector<vector<double> > &mat );
void writeMat( const vector<vector<double> > &mat, ofstream &matfile );
void writeTmat( const vector<vector<double> > &mat, ofstream &matfile );
void logSumExp( const vector<double> &vec, double &lse, const int &hp );
void postDecode(
        const vector<vector<double> > &fwd,
        const vector<vector<double> > &bwd,
        const hmmStates &st,
        const double &Px,
        vector<vector<double> > &pprob,
        pathVec &pvec,
        ofstream &logfile);

void viterbi(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &vit,
        pathVec &pvec );
void viterbi2(
        const vector<vector<int> > &sites,
        const positions &pos,
        const parameters &p,
        const hmmStates &st,
        const vector<int> &obs,
        const vector<double> &sprob,
        vector<vector<double> > &vit,
        pathVec &pvec );

void max( const vector<double> &vec, double maxelement );

void pathOutput( pathVec &pvec, hmmStates &st, positions &pos, vector<vector<double> > &pprob, parameters &param, ofstream &pathfile );

#endif


