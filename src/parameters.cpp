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
    mapFile = "";
    fixedMapRate = -1;
    n1 = 0; // set later
    n2 = 0; // set later
    S = 0; // set later
    T = 7.0;
    u1 = 0.2;
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
        else if (in_str == "--mapFile") { mapFile = get_arg(i+1); i++; }
        else if (in_str == "--fixedMapRate") { fixedMapRate = atof( get_arg(i+1).c_str() ); i++; }
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
        else if (in_str == "--pathfile") { pathfile = get_arg(i+1).c_str(); i++; }
        else
            error("Unknown option: " + string(in_str), 0);
        i++;
    }
    if( ( fname != "unset" ) && ( outfname == "unset" ) ) { 
        outfname = fname; 
    //} else if( ( fname == "unset" ) && ( admix == "unset" ) || ( ref == "unset" ) ) { 
    //    cerr << "Both --admix and --ref must be set!" << endl;
    //    exit(1);
    } else if( ( ( admix != "unset" ) && ( ref != "unset" ) ) && ( outfname == "unset" ) ) { 
        outfname = admix; 
    }
    logf = outfname + ".log" + std::to_string(mode);
    pathf = outfname + ".path" + std::to_string(mode);
    matf = outfname + ".mat" + std::to_string(mode);
    if( fixedMapRate != -1 ) {
        mapFile = "None";
    }
    //check_parameters();
}

void parameters::print_params(ofstream &logfile, const int which) {
    if( which == 0 ) {
        logfile << "CGmix: Crossover and gene conversion detection in admixed populations." << endl;
        logfile << "Input command: ";
        for(int i=0; i < argv.size(); i++ ) {
            logfile << argv[i] << " ";
        }
        logfile << endl;
        logfile << "Mode " << mode;
        if( mode == 0 ) { logfile << ". Haplotype-only model." << endl; }
        if( mode == 1 ) { logfile << ". Full haplotype and gene conversion model." << endl; }
        if( mode == 2 ) { logfile << ". Two-pass model." << endl; }
        if( mode == 3 ) { logfile << ". Second-pass model (requires output from mode 0).." << endl; }
        logfile << "I/O parameters set:" << endl;
        if( fname != "unset" )
            logfile << "Input file: " << fname << "[(.sites|.locs|.hapnames)]" << endl;
        if( ref != "unset" )
            logfile << "Reference input file: " << ref << "[(.sites|.locs|.hapnames)]" << endl;
        if( admix != "unset" )
            logfile << "Admix input file: " << admix << "[(.sites|.locs|.hapnames)]" << endl;
        if( fixedMapRate == -1 ) {
            logfile << "Genetic map: " << mapFile << endl;
        } else {
            logfile << "Using fixed recombination rate of: " << fixedMapRate << " cM/Mb"  << endl;
        }
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

void parameters::print_help() {
    string in_str;
    if( argv.size() <= 1 ) {
        argv.push_back("--?");
    }
    for(int i=0; i < argv.size(); i++) {
        in_str = argv[i];
        if ((in_str == "-h") || (in_str == "-?") || (in_str == "-help") || (in_str == "--?") || (in_str == "--help") || (in_str == "--h")) {
            cout << "CGmix: Crossover and gene conversion detection in admixed populations." << endl;
            cout << endl;
            cout << "Input files:" << endl;
            /////
            /////
            cout << "\t--ref" << "\t\tFile prefix for the reference haplotypes data:" << endl;
            cout << "\t\t\t(.sites) file: Each columns contains a haplotype, each row a variant coded as 0 or 1, whitespace delimited." << endl;
            cout << "\t\t\t(.locs)  file: One column, each line gives a physical position in integer base pairs on the chromosome." << endl;
            cout << "\t\t\t(.hapnames)  file: Two columns, the first gives the haplotype name, the second give the population it belongs to ('p1' or 'p2')." << endl;
            cout << "\t--admix" << "\t\tFile prefix for the admixed haplotype data.  Same 3 files as above, except only one admixed haplotype is allowed and it must be labeled as population 'p3'." << endl;
            cout << "\t--in" << "\t\tAlternate input format from --ref/--admix options above. Same 3 files as above, except the admixed haplotype is included in the same file as the reference haplotypes and the 'sites' and 'locs' files are transposed." << endl;
            /////

            /////
            cout << "Model parameters:" << endl;
            cout << "\t--mode" << "\t\tHMM run mode:" << endl;
            cout << "\t\t\t0 runs a first-pass using only the crossover chain." << endl;
            cout << "\t\t\t1 runs the full haplotype and gene conversion model. Potentially very slow." << endl;
            cout << "\t\t\t2 runs the first pass mode (0) to select interesting sites, then restarts, applying the full model only on those sites." << endl;
            cout << "\t\t\t3 takes input from a previous first pass run (mode 0) to mark interesting sites, and applies the full model only on those sites." << endl;
            cout << "\t--T" << "\t\tNumber of generations since admixture. [Default = " << T << "]." << endl;
            cout << "\t--u1" << "\t\tAncestry contribution from pop 1. [Default, Europeans = " << u1 << "]." << endl;
            cout << "\t--Ne1" << "\t\tEffective population size for pop 1. [Default, Europeans = " << Ne1 << "]." << endl;
            cout << "\t--Ne2" << "\t\tEffective population size for pop 2. [Default, Africans = " << Ne2 << "]." << endl;
            cout << "\t--f" << "\t\tRatio of gene conversion to crossover, g / r. [Default = " << f << "]." << endl;
            cout << "\t--lam" << "\t\tLambda, rate of terminating a gene conversion tract (1/lambda = tract length in kb). [Default = " << lam << "]." << endl;
            cout << "\t--mapFile" << "\tPath to a genetic map file (gzipped). Must follow the same convention as the Hapmap2 genetic map file." << endl;
            cout << "\t--fixedMapRate" << "\tInstead of using a genetic map, use this value as a fixed estimate of the genome-wide recombination rate (usually 1.1 cM/Mb)." << endl;
            /////
            cout << "Output parameters:" << endl;
            cout << "\t--out" << "\t\tOutput file prefix. Optional, defaults to value of '--admix'. CGmix produces two files:" << endl;
            cout << "\t\t\t(.log[mode]) file: Contains all parameters and log info for the run." << endl;
            cout << "\t\t\t(.path[mode]) file: Contains model outputs for each site." << endl;
            cout << "\t--matrixOutput" << "\tOutput forward, backward, and posterior matrices. Can use a lot of disk space." << endl;

            exit(0);
        }
    }
}

