/*
 * parameters.h
 *
 *  Created on: Jan 2, 2015
 *      Author: ccampbell
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <vector>
#include <iostream>
#include <stdint.h>
#include <unistd.h>
#include <fstream> // read files

using namespace std;

class parameters
{
    public:
        int mode;
        string fname;
        string outfname; // determined
        string logf;    // determined
        string pathf;   // determined
        string matf;    // determined
        string gmfile;
        int n1;         // determined
        int n2;         // determined
        int S;          // determined
        double T;
        double u1;
        double u2;
        double Ne1;
        double Ne2;
        double f;
        double lam;
        double theta1;  // determined
        double theta2;  // determined
        double theta1_match;    // determined
        double theta1_mismatch; // determined
        double theta2_match;    // determined
        double theta2_mismatch; // determined
        int highAccuracy; // whether to use sort in logSumExp (1) or not (0)
        int passAcc; // if two-pass model, use this for force sort on first pass
        int viterbi; // whether to run viterbi or not
        int fixPswitch; // how many sites (from index 0) to force into full model. Testing purposes only.
        int matrixOutput; // whether to output full forward/backward/posterior matrices

        parameters(int argc, char *argv[]);
        ~parameters(){};
        void read_parameters();
        void print_help();
        void print_params(ofstream &logfile, const int which);

    private:
        //void check_parameters();
        static void error(string err_msg, int code);
        vector<string> argv;
        string get_arg(unsigned int i);
};

#endif

