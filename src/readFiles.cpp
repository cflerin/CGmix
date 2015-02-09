/*
 * readFiles.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: ccampbell
 */

#include "readFiles.h"

/*
void print1Dvec(const vector<double> &vec) {
    for (vector<double>::const_iterator iter = vec.begin();
        iter != vec.end(); ++iter) {
        cout << *iter << " ";
    }
}
*/

void print1Dvec(const vector<int> &vec) {
    for (vector<int>::const_iterator iter = vec.begin();
        iter != vec.end(); ++iter) {
        cout << *iter << " ";
    }
}

void print2Dvec(const vector<vector<int> > &vec) {
for(int i=0; i < vec[0].size(); i++) {
      for (int j=0; j < vec.size(); j++)
        cout << vec[j][i] << " "; 
      cout << endl;
   }
}

void print2DvecString(const vector<vector<string> > &vec) {
for(int i=0; i < vec.size(); i++) {
      for (int j=0; j < vec[i].size(); j++)
        cout << vec[i][j] << " "; 
      cout << endl;
   }
}

void print1DvecString(const vector<string> &vec) {
    for (vector<string>::const_iterator iter = vec.begin();
        iter != vec.end(); ++iter) {
        cout << *iter << " ";
    }
}

void readSites(const string &fname, vector<vector<int> > &sites ) {
    ifstream file ( fname.c_str() );
    unsigned int cnt = 0;
    if (file.is_open()) {
        while (file.good()) {
            string line;
            getline(file, line);
            if( line == "" )
                continue;
            stringstream ss(line);
            string field;
            if( cnt == 0 ) { // read in first row:
                vector<int> row;
                while (getline(ss, field, '\t')) {
                    stringstream fs(field);
                    int f = 0;
                    fs >> f;
                    row.push_back(f);
                }
                // fill vector of vectors:
                for(int j=0; j < row.size(); j++ ) {
                    vector<int> tmp;
                    tmp.push_back( row[j] );
                    sites.push_back( tmp );

                }
            } else { // read the rest of the file:
                int j = 0;
                while (getline(ss, field, '\t')) {
                    stringstream fs(field);
                    int f = 0;
                    fs >> f;
                    sites[j].push_back(f);
                    j++;
                }
            }
            cnt++;
        }
        file.close();
    } else
        cout << "Failed to open sites file: " << fname << endl;
}

void readHapInfo(const string &fname, vector<vector<string> > &hapInfo ) {
    ifstream file ( fname.c_str() );
    if (file.is_open()) {
        while (file.good()) {
            vector<string> row;
            string line;
            getline(file, line);
            if( line == "" )
                continue;
            stringstream ss(line);
            string field;
            while( getline(ss, field, '\t')) {
                stringstream fs(field);
                string f = "";
                fs >> f;
                row.push_back(f);
            }
            hapInfo.push_back(row);
        }
        file.close();
    } else
        cout << "Failed to open haplotype info file: " << fname << endl;
}

void readLocs(const string &fname, vector<double> &locs ) {
    ifstream file ( fname.c_str() );
    if (file.is_open()) {
        while (file.good()) {
            string line;
            getline(file, line);
            if( line == "" )
                continue;
            stringstream ss(line);
            string field;
            while (getline(ss, field, '\t')) {
                stringstream fs(field);
                double f = 0.0;
                fs >> f;
                locs.push_back(f);
            }
        }
        file.close();
    } else
        cout << "Failed to open loc file: " << fname << endl;
}

void readGenMap(const string &fname, geneticMap &gMap ) {
    unsigned int gzMAX_LINE_LEN = 1024*1024;
    char *gz_readbuffer = new char[gzMAX_LINE_LEN];
    gzFile gzfile_in = gzopen(fname.c_str(), "rb");
    unsigned int lineCnt = 0;
    while( ! gzeof(gzfile_in) ) {
        // read 1 line:
        bool again = true;
        char * tmp;
        string out = "";
        while (again == true) {
            tmp = gzgets(gzfile_in, gz_readbuffer, gzMAX_LINE_LEN);
            if (tmp == NULL)
                return;
            out.append(gz_readbuffer);
            if ((strlen(gz_readbuffer) != gzMAX_LINE_LEN-1) || (gz_readbuffer[gzMAX_LINE_LEN-2] == '\n'))
                again = false;
        }
        out.erase( out.find_last_not_of(" \t\n\r") + 1); // Trim whitespace at end of line (required in gzipped case!)
        // end 1 line
        if( lineCnt == 0 ) { // skip header line
            lineCnt++;
            continue;
        }
        // fill genMap values:
        string field;
        stringstream ss(out);
        double tmpDbl;
        unsigned int fieldCnt = 0;
        while( std::getline( ss, field, '\t' ) ) {
            if( fieldCnt == 0 ) { 
                gMap.chr.push_back( field ); 
            }
            if( fieldCnt == 1 ) { 
                istringstream( field ) >> tmpDbl;
                gMap.pos.push_back( tmpDbl );
            }
            if( fieldCnt == 2 ) { 
                istringstream( field ) >> tmpDbl;
                gMap.rate.push_back( tmpDbl ); 
            }
            if( fieldCnt == 3 ) { 
                istringstream( field ) >> tmpDbl;
                gMap.cM.push_back( tmpDbl ); 
            }
            fieldCnt++;
        }
        // cout << gMap.chr[ gMap.chr.size()-1 ] << "\t" << gMap.pos[ gMap.pos.size()-1 ] << "\t" << gMap.rate[ gMap.rate.size()-1 ] << "\t" << gMap.cM[ gMap.cM.size()-1 ] << endl;
        lineCnt++;
    }
    gzclose( gzfile_in );
}

void interpGenMap(const geneticMap &gMap, const vector<double> &locs, positions &pos ) {
    pos.pos.resize( locs.size() );
    pos.cM.resize( locs.size() );
    unsigned int jstart = 0;
    double x, tmp;
    vector<int> flank(2), flankix(2);
    for(int i=0; i < locs.size(); i++ ) {
        x = locs[i];
        flank[0] =   std::numeric_limits<int>::max();
        flank[1] = - std::numeric_limits<int>::max();
        std::fill( flankix.begin(), flankix.end(), 0 );
        for(int j=jstart; j < gMap.pos.size(); j++ ) {
            tmp = x - gMap.pos[j];
            if( ( tmp>0 ) && ( tmp < flank[0] ) ) {
                flank[0] = gMap.pos[j];
                flankix[0] = j;
                continue;
            }
            if( ( tmp<0 ) && ( tmp > flank[1] ) ) {
                flank[1] = gMap.pos[j];
                flankix[1] = j;
                if(j>5) jstart = j - 5;
                break;
            }
        }
        pos.pos[i] = locs[i];
        pos.cM[i] =
            gMap.cM[ flankix[0] ] + ( gMap.cM[ flankix[1] ] - gMap.cM[ flankix[0] ] ) * ( x - flank[0] ) / ( flank[1] - flank[0] );
        //cout << pos.pos[i] << "\t" << pos.cM[i] << endl;
    }
}

