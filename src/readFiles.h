/*
 * readFiles.h
 *
 *  Created on: Nov 7, 2014
 *      Author: ccampbell
 */

#ifndef READFILES_H_
#define READFILES_H_

#include <vector>
//#include <string>
#include <iostream>
#include <fstream> // read files
#include <sstream>
#include <zlib.h>
//#include <stdio.h>
#include <cstring>
#include <limits>

using namespace std;

class hapDef 
{
    public:
        hapDef(){};
        ~hapDef(){};
        vector<string> hapName, hapPop;
        vector<int> hN, hP;
};
class geneticMap 
{
    public:
        geneticMap(){};
        ~geneticMap(){};
        vector<string> chr;
        vector<double> pos, rate, cM;
};
class positions 
{
    public:
        positions(){};
        ~positions(){};
        vector<double> pos, cM;
};

void readSites(const string &sitefname, vector<vector<int> > &sites );
void readLocs(const string &locfname, vector<double> &locs );
void readHapInfo(const string &fname, vector<vector<string> > &hapInfo );

void print1Dvec(const vector<int> &vec);
void print2Dvec(const vector<vector<int> > &vec);
void print2DvecString(const vector<vector<string> > &vec);
void print1DvecString(const vector<string> &vec);

void readGenMap(const string &fname, geneticMap &gMap );
void interpGenMap(const geneticMap &gMap, const vector<double> &locs, positions &pos );

#endif


