// Copyright 2017-2020 Edgar Costa
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
#include <sstream>
#include <NTL/vector.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>


#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace NTL;

void zeta_function(
    ZZX &zeta, // the zeta function
    Mat<ZZ> &Frob_ZZ, // the Frobenius matrix
    const map< Vec<int64_t>, ZZ, vi64less> &f, //f as a vector
    const int64_t &p, //the prime p
    const bool &verbose, //enable/disable verbose mode
    const int &threads, //number of threads
    const int &abs_precision, // in case we want Frob correct at least mod p^abs_precision,
    const bool &increase_precision_to_deduce_zeta, // in case we want compute Frob with enough precision to deduce the zeta function
    const bool &find_better_model // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    )
{
    #ifdef _OPENMP
    omp_set_num_threads(threads);
    #else
    assert(threads > 0);
    #endif
    Vec<int64_t> monomial = f.cbegin()->first;
    int64_t n = monomial.length() - 1;
    int64_t d = 0;
    for(int64_t i = 0; i < n + 1; ++i)
        d += monomial[i];


    Vec<int64_t> N;
    Vec<int64_t> charpoly_prec;
    int64_t precision;
    bool lift_charpoly = true;
    //populates precision, N, and charpoly_prec
    default_args(precision, N, charpoly_prec, p, n, d);
    if(abs_precision > 0) {
      // lets figure out an upper bound for the relative precision vector
      Vec<int64_t> r_vector;
      r_vector.SetLength(n);
      for(int64_t i = 0; i < n; ++i) {
        r_vector[i] = N[i] + (i + 1) - n;
      }

      // we can now modify r_vector accordingly
      if(increase_precision_to_deduce_zeta) {
        for(int64_t i = 0; i < n; ++i) {
          r_vector[i] = max(r_vector[i], abs_precision - (n - 1) + i);
          assert(r_vector[i] + n - (i + 1) >=  N[i]);
          N[i] = r_vector[i] + n - (i + 1);
          precision = max(precision, r_vector[i] + max(int64_t(0), valuation_of_factorial(p * (i + N[i]) - 1, p) - i));
        }
      } else {
        precision = 1;
        for(int64_t i = 0; i < n; ++i) {
          // the columns corresponding to H^(i,(n-1) - i) have valuation (n-1) - i
          int relative_prec = abs_precision - (n - 1) + i;
          // if we drop the relative prec of a column, we can no longer lift
          lift_charpoly = lift_charpoly && (relative_prec >= r_vector[i]);
          if(relative_prec <= 0) {
            r_vector[i] = N[i] = 0; //these columns will be zero regardless
          } else {
            r_vector[i] = abs_precision - (n - 1) + i;
            N[i] = r_vector[i] + n - (i + 1);
            precision = max(precision, r_vector[i] + max(int64_t(0), valuation_of_factorial(p * (i + N[i]) - 1, p) - i));
          }
        }
      }
      assert(abs_precision <= precision);
      // too messy for small p
      if( p < 2*n + max(r_vector) ) {
        cout <<"p is too small, optional parameters on precision only implemented for p >= 2*n + r = ";
        cout << 2*n + max(r_vector) << endl;
        cout<<"bye bye"<<endl;
        abort();
      }
    }

    {
        zz_pPush push(p);
        ZZ_pPush push2(power_ZZ(p,precision));

        map< Vec<int64_t>, zz_p , vi64less> fp;
        for(map< Vec<int64_t>, ZZ, vi64less>::const_iterator fit = f.begin(); fit != f.end(); ++fit)
            fp[ fit->first ] = conv<zz_p>(fit->second);

        map< Vec<int64_t>, ZZ_p, vi64less> f_ZZp;
        for(map< Vec<int64_t>, ZZ, vi64less>::const_iterator fit = f.begin(); fit != f.end(); ++fit)
            f_ZZp[ fit->first ] = conv<ZZ_p>(fit->second);

        if(!isSmooth(fp))
        {
            cerr << "f is not smooth!" <<endl;
            abort();
        }
        bool is_ND = isND(fp);
        if(!is_ND and find_better_model) {
          Mat<zz_p> M = find_change_of_variables(fp, p*100 + 100);
          if(!IsZero(M)) {
            if (verbose) {
                cout<<"Found a change of variables!"<<endl;
            }
            Mat<ZZ_p> M_ZZ_p = conv< Mat<ZZ_p> >(conv<Mat<ZZ> >(M));
            f_ZZp = change_of_variables<ZZ_p>(f_ZZp, M_ZZ_p);
            is_ND = true;
          } else {
            if(verbose) {
                cout<<"Wasn't able to find a suitable change of variables to make it non degenerate!"<<endl;
            }
          }

        }
        Mat<ZZ_p> Frob;

        if( is_ND )
        {
            hypersurface_non_degenerate hs_ND(p, precision, f_ZZp, verbose);
            assert(hs_ND.dR->coKernels_J_basis.length() + 1 == charpoly_prec.length() );
            Frob = hs_ND.frob_matrix_ND(N);
            // FIXME?
            //delete hs_ND.dR_ND;
        }
        else
        {
            hypersurface hs(p, precision, f_ZZp, verbose);
            assert(hs.dR->coKernels_J_basis.length() + 1 == charpoly_prec.length() );
            Frob = hs.frob_matrix_J(N);
            // FIXME?
            //delete hs.dR;
        }
        Frob_ZZ = conv< Mat<ZZ> >(Frob);
        Vec<ZZ> cp;
        if (verbose)
            cout << "Frob = "<<Frob<<endl;

        if(lift_charpoly) {
          cp = charpoly_frob(Frob_ZZ, charpoly_prec, p, n - 1);
          if(verbose)
            cout <<"Characteristic polynomial = "<< cp <<endl;
        } else {
          cp = charpoly(Frob_ZZ); // not enough precision
          // reduce the output precision accordingly
          ZZ p_abs_precision = power_ZZ(p, abs_precision);
          for(size_t i=0; i < cp.length(); ++i) {
            cp[i] %= p_abs_precision;
          }
        }

        zeta = conv<ZZX>(cp);
    }
}


void zeta_function(
    ZZX &zeta, // the zeta function
    Mat<ZZ> &Frob_ZZ, // the Frobenius matrix
    const vector< vector<int64_t> > &monomials, // monomials of f
    const vector<ZZ> &coefficients, // coefficient of f
    const int64_t &p, //the prime p
    const bool &verbose, //enable/disable verbose mode
    const int &threads, //number of threads
    const int &abs_precision, // in case we want Frob correct mod p^abs_precision,
    const bool &increase_precision_to_deduce_zeta, // in case we want Frob with enough precision to deduce the zeta function
    const bool &find_better_model // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    )
{
    map< Vec<int64_t>, ZZ, vi64less> f_map;
    assert( monomials.size() == coefficients.size() );
    for(uint64_t i = 0; i < monomials.size(); ++i) {
        Vec<int64_t> monomial;
        monomial.SetLength(monomials[i].size(), 0);
        for(uint64_t j = 0; j < monomials[i].size(); ++j)
            monomial[j] = monomials[i][j];
        f_map[ monomial ] = coefficients[i];
    }
    zeta_function(zeta, Frob_ZZ, f_map,
        p, verbose, threads,
        abs_precision, increase_precision_to_deduce_zeta,
        find_better_model);
}

void zeta_function(
    ZZX &zeta, // the zeta function
    Mat<ZZ> &Frob_ZZ, // the Frobenius matrix
    const vector< vector<int64_t> > &monomials,
    const vector<int64_t> &coefficients,
    const int64_t &p, //the prime p
    const bool &verbose, //enable/disable verbose mode
    const int &threads, //number of threads
    const int &abs_precision, // in case we want Frob correct mod p^abs_precision,
    const bool &increase_precision_to_deduce_zeta, // in case we want compute Frob with enough precision to deduce the zeta function
    const bool &find_better_model // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    )
{
  zeta_function(zeta, Frob_ZZ, monomials, vector<ZZ>(coefficients.begin(), coefficients.end()),
      p, verbose, threads,
      abs_precision, increase_precision_to_deduce_zeta,
      find_better_model);
}


/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      f.keys()
 *      f.values()
 */
void zeta_function(
    ZZX &zeta, // the zeta function
    Mat<ZZ> &Frob_ZZ, // the Frobenius matrix
    const char* input,
    const bool &verbose, //enable/disable verbose mode
    const int &threads, //number of threads
    const int &abs_precision, // in case we want Frob correct mod p^abs_precision,
    const bool &increase_precision_to_deduce_zeta, // in case we want compute Frob with enough precision to deduce the zeta function
    const bool &find_better_model // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    )
{
    int64_t p;
    map< Vec<int64_t>, ZZ, vi64less> f;
    std::stringstream buffer;
    buffer << input;
    buffer >> p;
    buffer >> f;
    zeta_function(zeta, Frob_ZZ, f,
        p, verbose, threads,
        abs_precision, increase_precision_to_deduce_zeta,
        find_better_model);
}
