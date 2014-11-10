/*
 * readFiles.cpp
 *
 *  Created on: Nov 7, 2014
 *      Author: ccampbell
 */

#include "readFiles.h"

/*
void print1Dvec(const vector<double>& vec) {
    for (vector<double>::const_iterator iter = vec.begin();
        iter != vec.end(); ++iter) {
        cout << *iter << " ";
    }
}
*/

void print1Dvec(const vector<int>& vec) {
    for (vector<int>::const_iterator iter = vec.begin();
        iter != vec.end(); ++iter) {
        cout << *iter << " ";
    }
}

void print2Dvec(const vector<vector<int> >& vec) {
for(int i=0; i < vec.size(); i++) {
      for (int j=0; j < vec[i].size(); j++)
        cout << vec[i][j] << " "; 
      cout << endl;
   }
}

void print2DvecString(const vector<vector<string> >& vec) {
for(int i=0; i < vec.size(); i++) {
      for (int j=0; j < vec[i].size(); j++)
        cout << vec[i][j] << " "; 
      cout << endl;
   }
}



void readSites(const string &fname, vector<vector<int> > &sites ) {
	ifstream file ( fname.c_str() );
	if (file.is_open()) {
		while (file.good()) {
			vector<int> row;
			string line;
			getline(file, line);
			stringstream ss(line);
			string field;
			while (getline(ss, field, '\t')) {
				stringstream fs(field);
				int f = 0;
				fs >> f;
				row.push_back(f);
			}
			// write_vec(row);
			// cout << endl;
			sites.push_back(row);
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
			stringstream ss(line);
			string field;
			while (getline(ss, field, '\t')) {
				stringstream fs(field);
				string f = "";
				fs >> f;
				row.push_back(f);
			}
			//write_vec(row);
			// cout << endl;
			hapInfo.push_back(row);
		}
		file.close();
	} else
		cout << "Failed to open haplotype info file: " << fname << endl;
}

void readLocs(const string &fname, vector<int> &locs ) {
	ifstream file ( fname.c_str() );
	vector<double> row;
	if (file.is_open()) {
		while (file.good()) {
			string line;
			getline(file, line);
			stringstream ss(line);
			string field;
			while (getline(ss, field, '\t')) {
				stringstream fs(field);
				int f = 0;
				fs >> f;
				row.push_back(f);
			}
			// write_vec(row);
			// cout << endl;
		}
		file.close();
	} else
		cout << "Failed to open loc file: " << fname << endl;
}


