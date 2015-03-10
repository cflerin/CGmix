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
    fname = "unset" ;
    ref = "unset" ;
    admix = "unset" ;
    outfname = "unset" ;
    logf = "unset";
    pathf = "unset";
    matf = "unset";
    gmfile = "/home/ccampbell/resources/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr22.txt.gz";
    n1 = 0; // set later
    n2 = 0; // set later
    S = 0; // set later
    T = 7.0;
    u1 = 0.5;
    Ne1 = 10000;
    Ne2 = 18000;
    f = 10.0; // scaling factor for g ( f=gamma/rho, g=f*r, and gamma=4*Ne*g )
    lam = 1.0/500 * 1000; // tract length (kb)
    theta1 = 0.0;
    theta2 = 0.0;
    fixPswitch = -1;
    highAccuracy = 0;
    viterbi = 0;
    matrixOutput = 0;
    gcsens = 7;
    width = 3;
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
        else if (in_str == "--out") { outfname = get_arg(i+1); i++; }
        else if (in_str == "--gmfile") { gmfile = get_arg(i+1); i++; }
        else if (in_str == "--mode") { mode = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--T") { T = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--u1") { u1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--Ne1") { Ne1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--Ne2") { Ne2 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--f") { f = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--lam") { lam = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta1") { theta1 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--theta2") { theta2 = atof( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--fixPswitch") { fixPswitch = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--highAccuracy") { highAccuracy = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--viterbi") { viterbi = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--matrixOutput") { matrixOutput = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--ref") { ref = get_arg(i+1).c_str(); i++; }
        else if (in_str == "--admix") { admix = get_arg(i+1).c_str(); i++; }
        else if (in_str == "--gcsens") { gcsens = atoi( get_arg(i+1).c_str() ); i++; }
        else if (in_str == "--width") { width = atoi( get_arg(i+1).c_str() ); i++; }
        else
            error("Unknown option: " + string(in_str), 0);
        i++;
    }
    if( ( fname != "unset" ) && ( outfname == "unset" ) ) { 
        outfname = fname; 
    //} else if( ( fname == "unset" ) && ( admix == "unset" ) || ( ref == "unset" ) ) { 
    //    cerr << "Both --admix and --ref must be set!" << endl;
    //    exit(1);
    } else if( ( ref == "unset" ) && ( outfname == "unset" ) ) { 
        outfname = admix; 
    }
    logf = outfname + ".log" + std::to_string(mode);
    pathf = outfname + ".path" + std::to_string(mode);
    matf = outfname + ".mat" + std::to_string(mode);
    //check_parameters();
}

void parameters::print_params(ofstream &logfile, const int which) {
    if( which == 0 ) {
        logfile << "CGmix: Crossover and Gene converion detection in admixed populations." << endl;
        logfile << "Mode " << mode;
        if( mode == 0 ) { logfile << ". Haplotype-only model." << endl; }
        if( mode == 1 ) { logfile << ". Full haplotype and gene conversion model." << endl; }
        if( mode == 2 ) { logfile << ". Two-pass model." << endl; }
        logfile << "I/O parameters set:" << endl;
        if( fname != "unset" )
            logfile << "Input file: " << fname << "[(.sites|.locs|.hapnames)]" << endl;
        if( ref != "unset" )
            logfile << "Reference input file: " << ref << "[(.sites|.locs|.hapnames)]" << endl;
        if( admix != "unset" )
            logfile << "Admix input file: " << admix << "[(.sites|.locs|.hapnames)]" << endl;
        logfile << "Log file: " << logf << endl;
        logfile << "Path file: " << pathf << endl;
        logfile << "Matrix output file: ";
        if( matrixOutput == 1 ) { 
            logfile << matf << endl;
        } else {
            logfile << "NA" << endl;
        }
        logfile << "Viterbi decoding: ";
        if( viterbi == 1 ) {
            logfile << "on" << endl;
        } else {
            logfile << "off" << endl;
        }
    }
    if( which == 1 ) {
        //
        logfile << endl << "Model parameters set:" << endl;
        logfile << "n1 = " << n1 << endl;
        logfile << "n2 = " << n2 << endl;
        logfile << "nSites = " << S << endl;
        logfile << "T = " << T << endl;
        logfile << "u1 = " << u1 << endl;
        logfile << "Ne1 = " << Ne1 << endl;
        logfile << "Ne2 = " << Ne2 << endl;
        logfile << "f = " << f << endl;
        logfile << "lambda = " << lam << endl;
        logfile << "theta1 = " << theta1 << endl;
        logfile << "theta2 = " << theta2 << endl;
        logfile << "fixPswitch = " << fixPswitch << endl;
        logfile << "highAccuracy = " << highAccuracy << endl;
        logfile << "GCsensitivity = " << gcsens << endl;
        logfile << "Expansion width = " << width << endl;
    }


}

