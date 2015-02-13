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
        //double rho1;
        //double rho2;
        //double gam1;
        //double gam2;
        double lam;
        double theta1;  // determined
        double theta2;  // determined
        double theta1_match;    // determined
        double theta1_mismatch; // determined
        double theta2_match;    // determined
        double theta2_mismatch; // determined

        // double rho;
        // double gam;
        // double theta;

        parameters(int argc, char *argv[]);
        ~parameters(){};

        void read_parameters();
        void print_help();
        //void print_params();
        void print_params(ofstream &logfile, const int which);

    private:
        //void check_parameters();
        static void error(string err_msg, int code);
        vector<string> argv;
        string get_arg(unsigned int i);
};





#endif
