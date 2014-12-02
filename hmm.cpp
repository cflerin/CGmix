/*
 * hmm.cpp
 *
 *  Created on: Nov 10, 2014
 *      Author: ccampbell
 */

#include "hmm.h"

using namespace std;

// void generateStates(const vector<vector<string> >& hapInfo, class hmmStates& st ) {
void generateStates(const class hapDef& hapInfo, class hmmStates& st ) {
    int nrow = hapInfo.hN.size() - 1 ;
    // G null states:
    for (int j=0; j < nrow; j++) { 
        // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << "0" << "\t" << "0" << "\t" << endl;
        st.Xhap.push_back( hapInfo.hN[j] );
        st.Xpop.push_back( hapInfo.hP[j] );
        st.Xindx.push_back( j ); // check if this is right
        st.Ghap.push_back( 0 );
        st.Gpop.push_back( 0 );
        st.Gindx.push_back( j ); // check if this is right
    }
    // All combinations of X and G:
    for(int i=0; i < nrow; i++) {
        for (int j=0; j < nrow; j++) {
            // cout << hapInfo[j][0] << "\t" << hapInfo[j][1] << "\t" << hapInfo[i][0] << "\t" << hapInfo[i][1] << "\t" << endl;
            st.Xhap.push_back( hapInfo.hN[j] );
            st.Xpop.push_back( hapInfo.hP[j] );
            st.Xindx.push_back( j );
            st.Ghap.push_back( hapInfo.hN[i] );
            st.Gpop.push_back( hapInfo.hP[i] );
            st.Gindx.push_back( i );
        }
    }
    // state concatenation:
    nrow = st.Xhap.size();
    for(int i=0; i < nrow; i++) {
        string tmp;
        tmp = std::to_string(st.Xpop[i]) + "-" + std::to_string(st.Xhap[i]) + "-" + std::to_string(st.Gpop[i]) + "-" + std::to_string(st.Ghap[i]);
        //tmp = st.Xpop[i] + "-" + st.Xhap[i] + "-" + st.Gpop[i] + "-" + st.Ghap[i];
        st.states.push_back( tmp );
        // cout << i << "\t" << st.Xhap[i] << "\t" << st.Xpop[i] << "\t" << st.Xindx[i] << "\t" <<
        //      st.Ghap[i] << "\t" << st.Gpop[i] << "\t" << st.Gindx[i] << "\t" << st.states[i] << endl;
    }
}

void getXtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trX ) {
    double u;
    int nm;
    if( st.Xpop[to] == 1 ) {
        u = p.u1;
        nm = p.n1;
    }
    if( st.Xpop[to] == 2 ) {
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
    trX = log(trX);
}

void getGtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, double& trG ) {
    double u, a1, a2, a3, a4, a5;
    int nm;
    if( st.Gpop[to] == 1 ) {
        u = p.u1;
        nm = p.n1;
    }
    if( st.Gpop[to] == 2 ) {
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
    if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] == 0 ) ) trG = a1;
    if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] != 0 ) ) trG = a2;
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] == 0 ) ) trG = a3;
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] != 0 ) && ( st.Ghap[from] != st.Ghap[to] ) ) trG = a5;
    trG = log(trG);
}

void getsprob( 
        const vector<int>& sites0, 
        const struct parameters& p, 
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        vector<double>& sprob ) {
    double e;
    double sprob_G0 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.lam * ( p.n1 + p.n2 ) ) / ( p.lam * (p.n1+p.n1) + p.gam * p.T) );
    double sprob_G1 = log( 1.0 / ( p.n1 + p.n2 ) * ( p.gam * p.T ) / ( ( p.n1 + p.n2 ) * ( p.lam * (p.n1 + p.n2 ) + p.gam * p.T ) ) );
    for(int i=0; i < st.states.size(); i++) {
        if( sites0[ st.Gindx[i] ] == obs[0] ) {
            e = emit.match;
        } else {
            e = emit.mismatch;
        }
        if( st.Ghap[i] == 0 ) {
            sprob[i] =  sprob_G0 + e;
        } else {
            sprob[i] =  sprob_G1 + e;
        }
    }
}

double lookupXtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, vector<double>& trXbin ) {
    int type;
    double trX;
    if( st.Xpop[from] != st.Xpop[to] ) {
        type = 0;
    } else {
        if( st.Xhap[from] != st.Xhap[to] ) {
            type = 1;
        } else { // if( st.Xhap[from] == st.Xhap[to] )
            type = 2;
        }
    }
    if( st.Xpop[to] == 2 ) // shift type index by 3 if transitioning to a haplotype in pop #2.
        type += 3;
    if( trXbin[ type ] == 99 ) { // this transition hasn't been set yet:
        getXtrans( to, from, d, st, p, trX);
        trXbin[ type ] = trX;
        return( trX );
    } else { // lookup and return a value already seen:
        trX = trXbin[ type ];
        return( trX );
    }
}

double lookupGtrans(const int& to, const int& from, const int& d, const class hmmStates& st, const struct parameters& p, class trBin& trGbin ) {
	double trG;
	int type;
    if( st.Ghap[from] == st.Ghap[to] )
        type = 4;
    if( ( st.Ghap[from] == st.Ghap[to] ) && ( st.Ghap[from] == 0 ) )
        type = 1; // overwrites above
    if( ( st.Ghap[from] == 0 ) && ( st.Ghap[to] == 0 ) )
        type = 2;
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] == 0 ) )
        type = 3;
    if( ( st.Ghap[from] != 0 ) && ( st.Ghap[to] != 0 ) && ( st.Ghap[from] != st.Ghap[to] ) )
        type = 5;
    if( trGbin.tr.size() == 0 ) { // start:
        getGtrans( to, from, d, st, p, trG);
        trGbin.tr.push_back( trG );
        trGbin.type.push_back( type );
        trGbin.toPop.push_back( st.Gpop[to] );
        // cout << endl << "started trG = " <<  trG << " size=" << trGbin.tr.size() << endl;
        return(trG);
    }
    // search for existing transition within trBin:
    for(int i=0; i < trGbin.tr.size(); i++ ) {
        if( ( trGbin.type[i] == type ) && ( trGbin.toPop[i] == st.Gpop[to] ) ) {
            trG = trGbin.tr[i];
            // cout << "reused trG = " <<  trG << " size=" << trGbin.tr.size() << endl;
            return(trG);
        }
    }
    // generate new transition if not found above:
    getGtrans( to, from, d, st, p, trG);
    trGbin.tr.push_back( trG );
    trGbin.type.push_back( type );
    trGbin.toPop.push_back( st.Gpop[to] );
    // cout << "added trG = " <<  trG << " size=" << trGbin.tr.size() << endl;
    return(trG);
}

void forward( 
        const vector<vector<int> >& sites,
        const vector<int>& dvec,
        const struct parameters& p,
        const struct emissions& emit,
        const class hmmStates& st,
        const vector<int>& obs,
        const vector<double>& sprob,
        vector<vector<double> >& fwd ) {
    // cout << fwd.size() << " " << fwd[0].size() << endl;
    // starting prob: 
    fwd[0] = sprob;
    int d;
    double lsum, tmp, trXa, trGa, e;
    double negInf = - std::numeric_limits<double>::infinity();
    //trBin trXbin, trGbin;
    vector<double> trXbin(6,99);
    trBin trGbin;
    int cnt;
    for(int j=1; j < sites.size() ; j++ ) {
        //cout << "XXXXXXX| j=" << j << endl;
        d = dvec[j] - dvec[j-1];
        //double test;
        for(int t=0; t < st.states.size(); t++) {
            lsum = negInf;
            for(int f=0; f < st.states.size(); f++) {
                //getXtrans( t, f, d, st, p, trXa);
                //test = trXa;
                trXa = lookupXtrans( t, f, d, st, p, trXbin );
                //cout << "from: " << st.states[f] << " to: " << st.states[t] << "\t";
                //cout << "realtrX= " << test << " lookuptrX= " << trXa << endl;
                //if( test != trXa ) {
                //    cout << "FAIL" << endl;
                //    exit(1);
                //}
                getGtrans( t, f, d, st, p, trGa);
                //trGa = lookupGtrans( t, f, d, st, p, trGbin );
                tmp = fwd[j-1][f] + trXa + trGa;
                if( tmp > negInf ) lsum = tmp + log( 1 + exp( lsum - tmp ) );
            } // end 'from' loop
            // emission prob:
            if( sites[j][ st.Gindx[t] ] == obs[j] ) {
                e = emit.match;
            } else {
                e = emit.mismatch;
            }
            fwd[j][t] = e + lsum;
            // cout << t << endl;
        } // end 'to' loop
        cnt = 0;
        for(int k=0; k < trXbin.size(); k++) {
            if( trXbin[k] == 99 )
                cout << trXbin[k] << "\t" ;
                cnt++;
        }
        // cout << "trXbin99: " <<  cnt << " trGbin: " << trGbin.tr.size() << endl;
        std::fill( trXbin.begin(), trXbin.end(), 99 );
        // trXbin.Xswitch.clear();
        // trXbin.Pswitch.clear();
        // trXbin.tr.clear();
        // trXbin.toPop.clear();
        trGbin.tr.clear();
        trGbin.type.clear();
        trGbin.toPop.clear();
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
	double e, lsum, tmp, trXa, trGa;
	double negInf = - std::numeric_limits<double>::infinity();
	int d;
	trBin trXbin, trGbin;
	for(int j = sites.size()-2; j >= 0 ; j-- ) {
		d = dvec[j+1] - dvec[j];
		for(int f=0; f < st.states.size(); f++ ) {
			lsum = negInf;
			for(int t=0; t < st.states.size(); t++ ) {
				getXtrans( t, f, d, st, p, trXa);
				getGtrans( t, f, d, st, p, trGa);
				//trXa = lookupXtrans( t, f, d, st, p, trXbin);
				//trGa = lookupGtrans( t, f, d, st, p, trGbin);
				// emission prob:
				if( sites[j+1][ st.Gindx[t] ] == obs[j+1] ) {
					e = emit.match;
				} else {
					e = emit.mismatch;
				}
				tmp = bwd[j+1][t] + trXa + trGa + e ;
				if( tmp > negInf ) lsum = tmp + log( 1 + exp( lsum - tmp ) );
			} // t loop
			bwd[j][f] = lsum;
		} // f loop
        // trXbin.Xswitch.clear();
        // trXbin.Pswitch.clear();
        // trXbin.tr.clear();
        // trXbin.toPop.clear();
        trGbin.tr.clear();
        trGbin.type.clear();
        trGbin.toPop.clear();
	} // j loop
}

void printMat( const vector<vector<double> >& mat ) {
    for(int i=0; i < mat[0].size(); i++ ) {
        for(int j=0; j < mat.size(); j++ ) {
            cout << setprecision(15) << mat[j][i] << " ";
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
	double Pxa, tmp;
	logSumExp( fwd[ fwd.size()-1 ], Pxa );
	cout << "Pxa= " << Pxa << endl;
    double negInf = - std::numeric_limits<double>::infinity();
	// backward:
	double Pxb = bwd[0][0] + fwd[0][0];
	for(int i=1; i < fwd[0].size(); i++ ) {
		tmp = bwd[0][i] + fwd[0][i];
		if( tmp > negInf ) 
            Pxb = tmp + log( 1 + exp( Pxb - tmp ) );
	}
	cout << "Pxb= " << Pxb << endl;
    // decoding:
    for(int j=0; j < fwd.size(); j++ ) {
        for(int i=0; i < fwd[0].size(); i++ ) {
            pprob[j][i] = fwd[j][i] + bwd[j][i] - Pxa;
        }
    }
    // normalize:
    double cmax, csum;
    for(int j=0; j < fwd.size(); j++ ) {
        cmax = negInf;
        // determine max value:
        for(int i=0; i < fwd[0].size(); i++ ) {
            if( pprob[j][i] > cmax ) 
                cmax = pprob[j][i];
        }
        // subtract max and determine sum:
        csum = 0.0;
        for(int i=0; i < fwd[0].size(); i++ ) {
            pprob[j][i] = pprob[j][i] - cmax;
            csum += exp( pprob[j][i] );
        }
        // divide by sum:
        for(int i=0; i < fwd[0].size(); i++ ) {
            pprob[j][i] = exp( pprob[j][i] ) / csum;
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
        const vector<double>& sprob,
		vector<vector<double> >& vit,
		vector<string>& vpath,
		vector<double>& vprob ) {
    // starting prob, same as fwd:
    vit[0] = sprob;
    // recursion:
    int d;
    double trXa, trGa, vmax, tmp, e;
    double negInf = - std::numeric_limits<double>::infinity();
	trBin trXbin, trGbin;
    for(int j=1; j < sites.size() ; j++ ) {
    	d = dvec[j] - dvec[j-1];
    	for(int t=0; t < st.states.size(); t++ ) {
    	    vmax = negInf;
    	    for(int f=0; f < st.states.size(); f++ ) {
    	    	getXtrans( t, f, d, st, p, trXa);
    	    	getGtrans( t, f, d, st, p, trGa);
				//trXa = lookupXtrans( t, f, d, st, p, trXbin);
				//trGa = lookupGtrans( t, f, d, st, p, trGbin);
                tmp = vit[j-1][f] + trXa + trGa;
                if( tmp > vmax )
                    vmax = tmp;
    	    } // end from loop
    	    // emission prob:
    	    if( sites[j][ st.Gindx[t] ] == obs[j] ) {
    	    	e = emit.match;
    	    } else {
    	    	e = emit.mismatch;
    	    }
    	    vit[j][t] = e + vmax;
    	} // end to loop
        // trXbin.Xswitch.clear();
        // trXbin.Pswitch.clear();
        // trXbin.tr.clear();
        // trXbin.toPop.clear();
        trGbin.tr.clear();
        trGbin.type.clear();
        trGbin.toPop.clear();
    } // end j loop
    // termination:
    int maxix = -1;
    // find last max state:
    vmax = negInf;
    maxix = -1;
    int j = sites.size()-1;
    for(int i=0; i < vit[0].size(); i++ ) {
        if( vit[j][i] > vmax ) {
            vmax = vit[j][i];
            maxix = i;
        }
    } 
    vpath[j] = st.states[maxix];
    vprob[j] = pprob[j][maxix];
    // cout << vpath[j] << " " << vmax << " " << maxix << endl;
    // traceback:
    int t;
    for(int j = sites.size()-2; j >= 0 ; j-- ) {
        // find previous max state:
        t = maxix;
        vmax = negInf;
        for(int f=0; f < st.states.size(); f++) {
            getXtrans( t, f, d, st, p, trXa);
            getGtrans( t, f, d, st, p, trGa);
			//trXa = lookupXtrans( t, f, d, st, p, trXbin);
			//trGa = lookupGtrans( t, f, d, st, p, trGbin);
            tmp = vit[j][f] + trXa + trGa;
            if( tmp > vmax ) {
                vmax = tmp;
                maxix = f;
            }
        }
        vpath[j] = st.states[ maxix ];
        vprob[j] = pprob[j][maxix];
        // cout << vpath[j] << " " << vmax << " " << maxix << endl;
        // trXbin.Xswitch.clear();
        // trXbin.Pswitch.clear();
        // trXbin.tr.clear();
        // trXbin.toPop.clear();
        trGbin.tr.clear();
        trGbin.type.clear();
        trGbin.toPop.clear();
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
