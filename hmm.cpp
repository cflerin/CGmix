/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

void generateStates(const vector<vector<string> >& hapInfo, class hmmStates& st ) {
    int nrow = hapInfo.size() - 1 ;
    // G null states:
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        st.Xhap.push_back( hapInfo[j][0] );
        st.Xpop.push_back( hapInfo[j][1] );
        st.Xindx.push_back( j ); // check if this is right
        st.Ghap.push_back( "0" );
        st.Gpop.push_back( "0" );
        st.Gindx.push_back( j ); // check if this is right
    }
    // All combinations of X and G:
    for(int i=0; i < nrow; i++) {
        for (int j=0; j < nrow; j++) {
            // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
            st.Xhap.push_back( hapInfo[j][0] );
            st.Xpop.push_back( hapInfo[j][1] );
            st.Xindx.push_back( j );
            st.Ghap.push_back( hapInfo[i][0] );
            st.Gpop.push_back( hapInfo[i][1] );
            st.Gindx.push_back( i );
        }
    }
    // state concatenation:
    nrow = st.Xhap.size();
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = st.Xpop[i] + "-" + st.Xhap[i] + "-" + st.Gpop[i] + "-" + st.Ghap[i];
        st.states.push_back( tmp );
        // cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" <<
        //      st.Ghap[i] << "\t" << st.Gpop[i] << "\t" << st.Gindx[i] << "\t" << st.states[i] << endl;
    }
}

void getXtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trX ) {
    double u;
    int nm;
    if( st.Xpop[to] == "p1" ) {
        u = p.u1;
        nm = p.n1;
    }
    if( st.Xpop[to] == "p2" ) {
        u = 1.0 - p.u1;
        nm = p.n2;
    }
    // cout << "From " << st.Xhap[from] << " To " << st.Xhap[to] ;
    // cout << " | " << u << " " << nm << " ";
    if( st.Xpop[from] != st.Xpop[to] ) { // ancestry and haplotype switch:
        // cout << " anc sw! ";
        trX = ( 1.0 - exp( -d * p.rho * p.T ) ) * u/nm;
    } else { // no ancestry switch
        if( st.Xhap[from] != st.Xhap[to] ) { // hap switch
            // cout << " hap sw! ";
            trX = exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
                  ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
        }
        if( st.Xhap[from] == st.Xhap[to] ) { // no hap switch
            // cout << " cont! ";
            trX = exp( -d * p.rho * p.T ) * exp( -d * p.rho ) + 
                  exp( -d * p.rho * p.T ) * ( 1.0 - exp( -d * p.rho ) ) / nm + 
                  ( 1.0 - exp( -d * p.rho * p.T )) * u / nm ;
        }
    }
}

void getGtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trG ) {
    double u, a1, a2, a3, a4, a5;
    int nm;
    if( st.Gpop[to] == "p1" ) {
        u = p.u1;
        nm = p.n1;
    }
    if( st.Gpop[to] == "p2" ) {
        u = 1.0 - p.u1;
        nm = p.n2;
    }
    // 3: Pr(G[j+1]=0 | G[j] = g):
    a3 = p.lam * (p.n1 + p.n2) * ( 1 - exp( -d * (p.gam * p.T + p.lam * (p.n1 + p.n2) ) / (p.n1 + p.n2) ) ) / (p.gam * p.T + p.lam * (p.n1 + p.n2) ) ; 
    // 5: Pr(G[j+1]=g' | G[j] = g):
    a5 = ( p.lam * u * (p.n1 + p.n2) * ( exp( d*(-p.gam * p.T / (p.n1 + p.n2) - p.lam)) - 1)) / ((p.n1 + p.n2) * nm * p.lam + p.gam * p.T * nm) + (u - u * exp(-d * p.lam)) / nm ;
    // 1: Pr(G[j+1]=0 | G[j] = 0):
    a1 = exp( -p.lam * d ) * exp( -p.gam * p.T * d / (p.n1 + p.n2)) + a3 ;
    // 2: Pr(G[j+1]=g | G[j] = 0):
    a2 = exp( -p.lam * d ) * ( 1 - exp( -p.gam * p.T * d / (p.n1 + p.n2))) * u / nm + a5 ;
    // 4: Pr(G[j+1]=g | G[j] = g):
    a4 = exp( -p.lam * d ) + a5 ;
    if( st.Ghap[from] == st.Ghap[to] )              trG = a4; // erroneously includes 0 -> 0; overwritten on next line
    if( ( st.Ghap[from] == "0" ) & ( st.Ghap[to] == "0" ) ) trG = a1;
    if( st.Ghap[from] == "0" & st.Ghap[to] != "0" ) trG = a2;
    if( st.Ghap[from] != "0" & st.Ghap[to] == "0" ) trG = a3;
    if( st.Ghap[from] != "0" & st.Ghap[to] != "0" & st.Ghap[from] != st.Ghap[to] ) trG = a5;
}

void forward( 
        const vector<vector<int> >& sites,
        const vector<int>& dvec,
        const struct parameters& p,
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        vector<vector<double> >& fwd ) {
    // cout << fwd.size() << " " << fwd[0].size() << endl;
    // starting prob: 
    double e;
    double sprob_G0 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * (p.n1+p.n1) + p.gam * p.T) );
    double sprob_G1 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.gam * p.T ) / ( ( p.n1 + p.n2 ) * ( p.lam * (p.n1 + p.n2 ) + p.gam * p.T ) ) );
    for(int i=0; i < st.states.size(); i++) {
        if( sites[0][ st.Gindx[i] ] == obs[0] ) {
            e = emit.match;
        } else {
            e = emit.mismatch;
        }
        if( st.Ghap[i] == "0" ) {
            fwd[0][i] =  sprob_G0 + e;
        } else {
            fwd[0][i] =  sprob_G1 + e;
        }
    }
    int d;
    double lsum, tmp, trX, trG;
    double negInf = - std::numeric_limits<double>::infinity();
    for(int j=1; j < sites.size() ; j++ ) {
        d = dvec[j] - dvec[j-1];
        for(int t=0; t < st.states.size(); t++) {
            lsum = negInf;
            for(int f=0; f < st.states.size(); f++) {
                getXtrans( t, f, d, st, p, trX);
                // cout << "\t" <<trX << "\t" << trG << endl;
                getGtrans( t, f, d, st, p, trG);
                tmp = fwd[j-1][f] + log( trX * trG );
                if( tmp > negInf ) lsum = tmp + log( 1 + exp( lsum - tmp ) );
            } // end 'from' loop
            // emission prob:
            if( sites[j][ st.Gindx[t] ] == obs[j] ) {
                e = emit.match;
            } else {
                e = emit.mismatch;
            }
            fwd[j][t] = e + lsum;
            // cout << fwd[t][j] << "\t";
        } // end 'to' loop
        // cout << endl;
    } // end j site loop
}

void backward(
		const vector<vector<int> >& sites,
		const vector<int>& dvec,
		const struct parameters& p,
		const struct emissions& emit,
		const class hmmStates& st,
		const vector<int>& obs,
		vector<vector<double> >& bwd ) {
	// starting prob:
	for(int i=0; i < st.states.size(); i++) {
		bwd[i][ bwd[i].size()-1 ] = 0.0;
	}
	double e, lsum, tmp, trX, trG;
	double negInf = - std::numeric_limits<double>::infinity();
	int d;
	for(int j = sites[0].size()-2; j >= 0 ; j-- ) {
		d = dvec[j+1] - dvec[j];
		for(int f=0; f < st.states.size(); f++ ) {
			lsum = negInf;
			for(int t=0; t < st.states.size(); t++ ) {
				getXtrans( t, f, d, st, p, trX);
				getGtrans( t, f, d, st, p, trG);
				// emission prob:
				if( sites[ st.Gindx[t] ][j+1] == obs[j+1] ) {
					e = emit.match;
				} else {
					e = emit.mismatch;
				}
				tmp = bwd[t][j+1] + log( trX * trG * e );
				if( tmp > negInf ) lsum = tmp + log( 1 + exp( lsum - tmp ) );
			} // t loop
			bwd[f][j] = lsum;
		} // f loop
	} // j loop
}


void printMat( const vector<vector<double> >& mat ) {
    for(int i=0; i < mat.size(); i++ ) {
        for(int j=0; j < mat[0].size(); j++ ) {
            cout << setprecision(15) << mat[i][j] << " ";
        }
        cout << endl;
    }
}

void logSumExp( const vector<double>& vec, double& lse ) {
	double max = - std::numeric_limits<double>::infinity();
	lse = 0.0;
	for(int i=0; i < vec.size(); i++ ) {
		if( vec[i] > max )
			max = vec[i];
	}
	for(int i=0; i < vec.size(); i++ )
		lse += exp( vec[i] - max );
	lse = max + log( lse) ;
}

void postDecode(
		const vector<vector<double> >& fwd,
		const vector<vector<double> >& bwd,
		vector<vector<double> >& pprob
		) {
	vector<double> sprob;
	for(int i=0; i < fwd.size(); i++) {
		sprob.push_back( fwd[i][ fwd[i].size()-1 ] );
		//cout << fwd[i][ fwd[i].size()-1 ] << " " << endl;
	}
	double Pxa, tmp;
	logSumExp( sprob, Pxa );
	cout << "Pxa= " << Pxa << endl;
    double negInf = - std::numeric_limits<double>::infinity();
	// backward:
	double Pxb = bwd[0][0] + fwd[0][0];
	for(int i=1; i < fwd.size(); i++ ) {
		tmp = bwd[i][0] + fwd[i][0];
		// Pxb += bwd[0][0] + fwd[0][ fwd[0].size()-1 ];
		if( tmp > negInf ) 
            Pxb = tmp + log( 1 + exp( Pxb - tmp ) );
	}
	cout << "Pxb= " << Pxb << endl;
    // decoding:
    for(int j=0; j < fwd[0].size(); j++ ) {
        for(int i=0; i < fwd.size(); i++ ) {
            pprob[i][j] = fwd[i][j] + bwd[i][j] - Pxa;
        }
    }
    // normalize:
    double cmax, csum;
    for(int j=0; j < fwd[0].size(); j++ ) {
        cmax = negInf;
        // determine max value:
        for(int i=0; i < fwd.size(); i++ ) {
            if( pprob[i][j] > cmax ) 
                cmax = pprob[i][j];
        }
        // subtract max and determine sum:
        csum = 0.0;
        for(int i=0; i < fwd.size(); i++ ) {
            pprob[i][j] = pprob[i][j] - cmax;
            csum += exp( pprob[i][j] );
        }
        // divide by sum:
        for(int i=0; i < fwd.size(); i++ ) {
            pprob[i][j] = exp( pprob[i][j] ) / csum;
        }
    }
}

void viterbi(
		const vector<vector<int> >& sites,
		const vector<int>& dvec,
		const struct parameters& p,
		const struct emissions& emit,
		const class hmmStates& st,
		const vector<int>& obs,
        const vector<vector<double> >& pprob,
		vector<vector<double> >& vit,
		vector<string>& vpath,
		vector<double>& vprob ) {
    // starting prob, same as fwd:
    double e;
    for(int i=0; i < st.states.size(); i++) {
        if( sites[ st.Gindx[i] ][0] == obs[0] ) {
            e = emit.match;
        } else {
            e = emit.mismatch;
        }
        // cout << st.Ghap[i] << endl;
        if( st.Ghap[i] == "0" ) {
            vit[i][0] = log( 1.0 / ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * (p.n1+p.n1) + p.gam * p.T) * e );
        } else {
            vit[i][0] = log( 1.0 / ( p.n1 + p.n2 ) * ( p.gam * p.T ) / ( ( p.n1 + p.n2 ) * ( p.lam * (p.n1 + p.n2 ) + p.gam * p.T ) ) * e );
        }
    }
    // recursion:
    int d;
    double trX, trG, vmax, tmp;
    double negInf = - std::numeric_limits<double>::infinity();
    // vector<double> tmp( sites.size()-1 );
    for(int j=1; j < sites[0].size() ; j++ ) {
    	d = dvec[j] - dvec[j-1];
    	for(int t=0; t < st.states.size(); t++ ) {
    	    vmax = negInf;
    	    for(int f=0; f < st.states.size(); f++ ) {
    	    	getXtrans( t, f, d, st, p, trX);
    	    	getGtrans( t, f, d, st, p, trG);
                tmp = vit[f][j-1] + log( trX * trG );
                if( tmp > vmax )
                    vmax = tmp;
    	    } // end from loop
    	    // emission prob:
    	    if( sites[ st.Gindx[t] ][j] == obs[j] ) {
    	    	e = emit.match;
    	    } else {
    	    	e = emit.mismatch;
    	    }
    	    vit[t][j] = log( e ) + vmax;
    	} // end to loop
    } // end j loop
    // termination:
    int maxix = -1;
    // find last max state:
    vmax = negInf;
    maxix = -1;
    int j = sites[0].size()-1;
    for(int i=0; i < vit.size(); i++ ) {
        if( vit[i][j] > vmax ) {
            vmax = vit[i][j];
            maxix = i;
        }
    } 
    vpath[j] = st.states[maxix];
    vprob[j] = pprob[j][maxix];
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites[0].size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < st.states.size(); f++) {
            getXtrans( t, f, d, st, p, trX);
            getGtrans( t, f, d, st, p, trG);
            tmp = vit[f][j] + log( trX * trG );
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        vpath[j] = st.states[ maxix ];
        vprob[j] = pprob[j][maxix];
        // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    }
}

void which_max( const vector<double>& vec, int maxindx ) {
	double tmp = - std::numeric_limits<double>::infinity();
    double maxelement;
	for(int i=0; i < vec.size(); i++ ) {
		tmp = vec[i];
		if( tmp > maxelement )
			maxelement = tmp;
	}
}

void max( const vector<double>& vec, double maxelement ) {
	double tmp = - std::numeric_limits<double>::infinity();
	for(int i=0; i < vec.size(); i++ ) {
		tmp = vec[i];
		if( tmp > maxelement )
			maxelement = tmp;
	}
}
