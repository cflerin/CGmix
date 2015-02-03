/*
 * parameters.cpp
 *
 *  Created on: Jan 2, 2015
 *      Author: ccampbell
 */

#include "parameters.h"


parameters::parameters(int argc, char *argv[]) {
	string tmp;
	for (int i=0; i<argc; i++) {
		tmp = argv[i];
		this->argv.push_back(tmp);
	}

    mode = 99;
    fname = "empty" ;
    outfname = "unset" ;
    logf = "empty";
    pathf = "empty";
    matf = "empty";
    gmfile = "/home/ccampbell/resources/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr22.txt.gz";

    n1 = 0; // set later
    n2 = 0; // set later
    S = 0; // set later
    T = 7.0;
    u1 = 0.5;
    rho1 = 600.0;
    rho2 = 900.0;
    gam1 = 0.0;
    gam2 = 0.0;
    lam = 1.0/500;
    theta1 = 0.0;
    theta2 = 0.0;
    theta3 = 0.01;

    rho = 600.0;
    gam = 100.0;
    theta = 1.0/1000;
}

string parameters::get_arg(unsigned int i) {
    if (i>=argv.size())
        error("Requested Missing Argument",76);
    return argv[i];
}

void parameters::error(string err_msg, int code) {
    // LOG.printLOG("\n\nError: " + err_msg + "\n\n");
    cout << "\nError: " << err_msg << "\n" << endl;
    exit(code);
}


void parameters::read_parameters() {
    unsigned int i=1;
    string in_str; 
    while ( i < argv.size() ) { 
        in_str = argv[i];
        if (in_str == "--in") { fname = get_arg(i+1); i++; }
        else if (in_str == "--gmfile") { gmfile = get_arg(i+1); i++; }
        else if (in_str == "--mode") { mode = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--T") { T = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--u1") { u1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--rho1") { rho1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--rho2") { rho2 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--gam1") { gam1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--gam2") { gam2 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--lam") { lam = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta1") { theta1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta2") { theta2 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta3") { theta3 = atof( get_arg(i+1).c_str() ); i++; }
        //////
        else if (in_str == "--rho") { rho = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--gam") { gam = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta") { theta = atof( get_arg(i+1).c_str() ); i++; }
        //////

        else
            error("Unknown option: " + string(in_str), 0);
        i++;
    }
    if( outfname == "unset" ) { outfname = fname; }
    logf = outfname + ".log" + std::to_string(mode);
    pathf = outfname + ".path" + std::to_string(mode);
    matf = outfname + ".mat" + std::to_string(mode);
    //check_parameters();
}

void parameters::print_params(ofstream &logfile, const int which) {
    if( which == 0 ) {
        //logfile << argv << endl;
        logfile << "Mode " << mode;
        if( mode == 0 ) { logfile << ". Haplotype-only model." << endl; }
        if( mode == 1 ) { logfile << ". Full haplotype and gene conversion model." << endl; }
        if( mode == 2 ) { logfile << ". Two-pass model." << endl; }
        logfile << "I/O parameters set:" << endl;
        logfile << "Input file: " << fname << "[(.sites|.locs|.hapnames)]" << endl;
        logfile << "Log file: " << logf << endl;
        logfile << "Path file: " << pathf << endl;
        logfile << "Matrix output file: " << matf << endl;
    }
    if( which == 1 ) {
        //
        logfile << endl << "Model parameters set:" << endl;
        logfile << "n1 = " << n1 << endl;
        logfile << "n2 = " << n2 << endl;
        logfile << "nSites = " << S << endl;
        logfile << "T = " << T << endl;
        logfile << "u1 = " << u1 << endl;
        logfile << "rho1 = " << rho1 << endl;
        logfile << "rho2 = " << rho2 << endl;
        logfile << "gamma1 = " << gam1 << endl;
        logfile << "gamma2 = " << gam2 << endl;
        logfile << "lambda = " << lam << endl;
        logfile << "theta1 = " << theta1 << endl;
        logfile << "theta2 = " << theta2 << endl;
        logfile << "theta3 = " << theta3 << endl;
        //
        logfile << "rho = " << rho << endl;
        logfile << "gamma = " << gam << endl;
        logfile << "theta = " << theta << endl;
    }


}
