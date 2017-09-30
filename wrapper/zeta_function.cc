// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "wrapper.h"
#include "dr.h"
#include "dr_nd.h"
#include "hypersurface.h"
#include "hypersurface_nd.h"
#include "tools.h"
#include "vec_int64.h"
#include <cstdint>
#include <map>
#include <vector>
#include <NTL/vector.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

using namespace std;
using namespace NTL;

void zeta_function(ZZX &zeta, const map< Vec<int64_t>, int64_t, vi64less> &f,const int64_t &p, bool verbose )
{
    
    Vec<int64_t> monomial = f.cbegin()->first;
    int64_t n = monomial.length() - 1;
    int64_t d = 0;
    for(int64_t i = 0; i < n + 1; ++i)
        d += monomial[i];

    Vec<int64_t> N;
    Vec<int64_t> charpoly_prec;
    int64_t precision;
    //populates precision, N, and charpoly_prec
    default_args(precision, N, charpoly_prec, p, n, d);
    {
        zz_pPush push(p);
        ZZ_pPush push2(power_ZZ(p,precision));

        map< Vec<int64_t>, zz_p , vi64less> f_zz_p;
        for(map< Vec<int64_t>, int64_t, vi64less>::const_iterator fit; fit != f.end(); ++fit)
            f_zz_p[ fit->first ] = conv<zz_p>(fit->second);

        if(!isSmooth(f_zz_p))
        {
            cout << "f is not smooth!" <<endl;
            abort();
        }
        // try to find a change of variables
        Mat<zz_p> M = find_change_of_variables(f_zz_p, p*1000 + 1000);
        bool is_ND = !IsZero(M);
        Mat<ZZ_p> Frob;
        if( is_ND )
        {
            if (verbose)
                cout<<"Found a change of variables!"<<endl;
            map< Vec<int64_t>, zz_p, vi64less> f_map = change_of_variables<zz_p>(f_zz_p, M);
            hypersurface_non_degenerate hs_ND(p, precision, f_map, verbose);
            assert(hs_ND.dR->coKernels_J_basis.length() + 1 == charpoly_prec.length() );
            Frob = hs_ND.frob_matrix_ND(N);
        }
        else
        {
            if(verbose)
                cout<<"Wasn't able to find a suitable change of variables to make it non degenerate!"<<endl;
            hypersurface hs(p, precision, f_zz_p, verbose);
            assert(hs.dR->coKernels_J_basis.length() + 1 == charpoly_prec.length() );
            Frob = hs.frob_matrix_J(N);
        }
        
        if (verbose)
            cout << "Frob = "<<Frob<<endl;
        Mat<ZZ> Frob_ZZ;
        Frob_ZZ = conv< Mat<ZZ> >(Frob);
        Vec<ZZ> cp = charpoly_frob(Frob_ZZ, charpoly_prec, p, n - 1);
        cout <<"Characteristic polynomial = "<< cp <<endl;
        zeta = conv<ZZX>(cp);
    }
}


void zeta_function(ZZX &zeta, vector< vector<int64_t> > monomials, vector<int64_t> coefficients, const int64_t &p, bool verbose )
{
    map< Vec<int64_t>, int64_t, vi64less> f_map;
    assert( monomials.size() == coefficients.size() );
    for(size_t i = 0; i < monomials.size(); ++i) {
        Vec<int64_t> monomial;
        monomial.SetLength(monomials[i].size(), 0);
        for(size_t j = 0; j < monomials[i].size(); ++j) 
            monomial[j] = monomials[i][j];
        f_map[ monomial ] = coefficients[i];
    }
    zeta_function(zeta, f_map, p, verbose);
}
